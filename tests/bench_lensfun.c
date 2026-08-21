/*
    LensSerious — is it actually faster, and where?

    Copyright (C) 2026 Aurélien PIERRE.  License: LGPL-3.0-or-later.
*/

/** @file bench_lensfun.c
 *
 * @brief The design claim of this library is that a lens correction is data, and that
 * treating it as data is faster than treating it as a library of callbacks. This measures
 * that claim against liblensfun, stage by stage, rather than as one number -- because the
 * stages differ by wildly different factors and a single figure would flatter it.
 *
 * Everything here is CPU-side on both sides, deliberately. The GPU numbers live in
 * tests/bench_opencl.c and in Ansel; comparing a GPU kernel against a CPU callback would
 * prove nothing about the design.
 *
 * The two per-pixel stages are measured BOTH single-threaded and across OpenMP threads,
 * because a pixel pipeline runs them in parallel and the two libraries need not scale the
 * same way. Rows are independent in both -- lensfun's modifier is const during apply and
 * writes only the caller's buffer, which is what lets Ansel parallelise it too -- so this
 * is a fair test of the same work, not of one library's threading model.
 *
 * Run it with the database built by tools/import_lensfun_xml.c:
 *     bench_lensfun lenses.db [megapixels]
 */

#include "lensserious_db.h"

#include <lensfun.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef _OPENMP
  #include <omp.h>
  #define LS_OMP_FOR _Pragma("omp parallel for schedule(static)")
  #define LS_THREADS omp_get_max_threads()
#else
  #define LS_OMP_FOR
  #define LS_THREADS 1
#endif

#if defined(__GNUC__) || defined(__clang__)
  #define LS_B_RESTRICT __restrict__
#else
  #define LS_B_RESTRICT
#endif

static double now_ms(void)
{
  struct timespec t;
  clock_gettime(CLOCK_MONOTONIC, &t);
  return t.tv_sec * 1e3 + t.tv_nsec / 1e6;
}

/* Each timed stage runs REPEATS times and the BEST is kept. A shared machine only ever
 * makes a measurement slower, never faster, so the minimum is the closest thing to the
 * cost of the code; the mean measures the machine's other tenants. Taking one sample gave
 * lensfun's single-threaded map as 333 ms in one run and 532 ms in the next. */
enum { REPEATS = 5 };

#define TIME_BEST(best, body)                                                             \
  do {                                                                                     \
    (best) = 1e30;                                                                         \
    for(int _r = 0; _r < REPEATS; _r++)                                                    \
    {                                                                                      \
      const double _t0 = now_ms();                                                         \
      body;                                                                                \
      const double _d = now_ms() - _t0;                                                    \
      if(_d < (best)) (best) = _d;                                                          \
    }                                                                                      \
  } while(0)

/* One bilinear tap per channel, each at its own coordinates -- the three channels of a
 * subpixel map do not land on the same pixel, which is the whole reason the map has six
 * floats and not two. Cheap on purpose: this stands in for a resampler so that the MEMORY
 * behaviour of fetching coordinates is what differs between the variants, not the filter. */
static inline void _resample_pixel(const float *LS_B_RESTRICT src, int w, int h,
                                   const float *LS_B_RESTRICT c6,
                                   float *LS_B_RESTRICT out)
{
  for(int ch = 0; ch < 3; ch++)
  {
    float fx = c6[ch * 2], fy = c6[ch * 2 + 1];
    fx = (fx < 0.f) ? 0.f : ((fx > (float)(w - 2)) ? (float)(w - 2) : fx);
    fy = (fy < 0.f) ? 0.f : ((fy > (float)(h - 2)) ? (float)(h - 2) : fy);
    const int ix = (int)fx, iy = (int)fy;
    const float tx = fx - (float)ix, ty = fy - (float)iy;
    const float *p00 = src + ((size_t)iy * w + ix) * 4 + ch;
    const float *p10 = p00 + 4;
    const float *p01 = p00 + (size_t)w * 4;
    const float *p11 = p01 + 4;
    out[ch] = (*p00 * (1.f - tx) + *p10 * tx) * (1.f - ty)
            + (*p01 * (1.f - tx) + *p11 * tx) * ty;
  }
  out[3] = 1.f;
}

static void _resample_row(const float *LS_B_RESTRICT src, int w, int h,
                          const float *LS_B_RESTRICT map, float *LS_B_RESTRICT dst, int count)
{
  for(int x = 0; x < count; x++)
    _resample_pixel(src, w, h, map + (size_t)x * 6, dst + (size_t)x * 4);
}

static void row(const char *what, double lf, double ls, const char *unit)
{
  if(lf > 0.0 && ls > 0.0)
    printf("  %-34s %10.4f %10.4f %s   %8.1fx\n", what, lf, ls, unit, lf / ls);
  else
    printf("  %-34s %10.4f %10.4f %s\n", what, lf, ls, unit);
}

int main(int argc, char **argv)
{
  if(argc < 2)
  {
    fprintf(stderr, "usage: %s <lenses.db> [megapixels]\n", argv[0]);
    return 2;
  }
  const double mpix = (argc > 2) ? atof(argv[2]) : 24.0;
  const int W = (int)(sqrt(mpix * 1e6 * 3.0 / 2.0) + 0.5);
  const int H = (int)(W * 2.0 / 3.0 + 0.5);

  /* A lens with distortion, TCA and vignetting all calibrated, so every stage has work. */
  const char *maker = "Nikon";
  const char *model = "Nikon AF-S Nikkor 16-35mm f/4G ED VR";
  const float focal = 16.f, aperture = 4.f, distance = 10.f, crop = 1.f;

  printf("bench: %s, %dx%d (%.1f Mpx), single-threaded\n", model, W, H, W * (double)H / 1e6);
  printf("  %-34s %10s %10s %s\n", "", "lensfun", "LensSerious", "");

  /* ======================================================================
   * PART ONE: getting the calibration data. Paid once per lens, on the CPU,
   * whatever the pixels are processed on.
   * ====================================================================== */

  printf("\n  getting the calibration data\n");

  /* --- 1. getting a database at all ------------------------------------- */

  double t = now_ms();
  lfDatabase *ldb = lf_db_new();
  if(!ldb || lf_db_load(ldb) != LF_NO_ERROR)
  {
    fprintf(stderr, "no system lensfun database\n");
    return 1;
  }
  const double t_lf_load = now_ms() - t;

  t = now_ms();
  ls_db_t *db = ls_db_open(argv[1]);
  const double t_ls_open = now_ms() - t;
  if(!db)
  {
    fprintf(stderr, "cannot open `%s'\n", argv[1]);
    return 1;
  }
  row("open the database", t_lf_load, t_ls_open, "ms");

  /* --- 2. finding the lens ---------------------------------------------- */

  enum { LOOKUPS = 200 };
  t = now_ms();
  for(int i = 0; i < LOOKUPS; i++)
  {
    const lfLens **r = lf_db_find_lenses_hd(ldb, NULL, maker, model, 0);
    if(r) lf_free(r);
  }
  const double t_lf_find = (now_ms() - t) / LOOKUPS;

  ls_lens_t lens;
  t = now_ms();
  for(int i = 0; i < LOOKUPS; i++)
  {
    ls_db_match_t m[1];
    if(ls_db_match_lens(db, maker, model, 0, m, 1) > 0) ls_db_lens_by_id(db, m[0].lens_id, &lens);
  }
  const double t_ls_find = (now_ms() - t) / LOOKUPS;
  row("find the lens (fuzzy)", t_lf_find, t_ls_find, "ms");

  const lfLens **hits = lf_db_find_lenses_hd(ldb, NULL, maker, model, 0);
  if(!hits || !hits[0])
  {
    fprintf(stderr, "lensfun cannot find `%s'\n", model);
    return 1;
  }
  const lfLens *lf = hits[0];
  {
    ls_db_match_t m[1];
    if(ls_db_match_lens(db, maker, model, 0, m, 1) <= 0
       || ls_db_lens_by_id(db, m[0].lens_id, &lens) != 1)
    {
      fprintf(stderr, "LensSerious cannot find `%s'\n", model);
      return 1;
    }
  }

  /* --- 3. resolving it at one shooting configuration -------------------- */

  enum { INITS = 2000 };
  t = now_ms();
  for(int i = 0; i < INITS; i++)
  {
    lfModifier *m = lf_modifier_new(lf, crop, W, H);
    lf_modifier_initialize(m, lf, LF_PF_F32, focal, aperture, distance, 1.f, LF_RECTILINEAR,
                           LF_MODIFY_DISTORTION | LF_MODIFY_TCA | LF_MODIFY_VIGNETTING, 0);
    lf_modifier_destroy(m);
  }
  const double t_lf_init = (now_ms() - t) / INITS;

  ls_modifier_t mod;
  t = now_ms();
  for(int i = 0; i < INITS; i++)
    ls_modifier_init(&mod, &lens, crop, W, H, focal, aperture, distance, 1.f,
                     LS_LENS_UNKNOWN,
                     LS_ENABLE_DISTORTION | LS_ENABLE_TCA | LS_ENABLE_VIGNETTING, 0);
  const double t_ls_init = (now_ms() - t) / INITS;
  row("resolve at focal/aperture", t_lf_init, t_ls_init, "ms");

  /* Open + find + resolve: everything between "I have a filename and some EXIF" and "I
   * can correct a pixel". This is the number a caller actually waits for, and the only
   * one of the three that is dominated by a term LensSerious does not have at all. */
  const double t_lf_data = t_lf_load + t_lf_find + t_lf_init;
  const double t_ls_data = t_ls_open + t_ls_find + t_ls_init;
  row("  = open + find + resolve", t_lf_data, t_ls_data, "ms");

  /* ======================================================================
   * PART TWO: processing the pixels. Paid per image, and the part a GPU can take.
   * ====================================================================== */

  printf("\n  processing the pixels (%.1f Mpx)\n", W * (double)H / 1e6);

  /* --- 4. the geometry map, the whole frame ----------------------------- */

  float *buf = (float *)malloc((size_t)W * 6 * sizeof(float));
  if(!buf) return 1;

  lfModifier *m = lf_modifier_new(lf, crop, W, H);
  lf_modifier_initialize(m, lf, LF_PF_F32, focal, aperture, distance, 1.f, LF_RECTILINEAR,
                         LF_MODIFY_DISTORTION | LF_MODIFY_TCA, 0);
  double t_lf_map;
  TIME_BEST(t_lf_map,
            for(int y = 0; y < H; y++)
              lf_modifier_apply_subpixel_geometry_distortion(m, 0.f, (float)y, W, 1, buf));

  ls_modifier_init(&mod, &lens, crop, W, H, focal, aperture, distance, 1.f, LS_LENS_UNKNOWN,
                   LS_ENABLE_DISTORTION | LS_ENABLE_TCA, 0);
  double t_ls_map;
  TIME_BEST(t_ls_map,
            for(int y = 0; y < H; y++)
              ls_modifier_apply_subpixel_geometry(&mod, 0.f, (float)y, W, 1, buf));
  row("build the whole geometry map", t_lf_map, t_ls_map, "ms");

  /* The same work across threads. Each thread writes its own row of the shared map, so the
   * only thing shared for writing is the buffer, and no two threads touch the same bytes. */
  double t_lf_map_mt = 0.0, t_ls_map_mt = 0.0;
#ifdef _OPENMP
  {
    float *big = (float *)malloc((size_t)W * H * 6 * sizeof(float));
    if(big)
    {
      TIME_BEST(t_lf_map_mt,
                _Pragma("omp parallel for schedule(static)")
                for(int y = 0; y < H; y++)
                  lf_modifier_apply_subpixel_geometry_distortion(m, 0.f, (float)y, W, 1,
                                                                 big + (size_t)y * W * 6));
      TIME_BEST(t_ls_map_mt,
                _Pragma("omp parallel for schedule(static)")
                for(int y = 0; y < H; y++)
                  ls_modifier_apply_subpixel_geometry(&mod, 0.f, (float)y, W, 1,
                                                      big + (size_t)y * W * 6));
      free(big);
    }
  }
#endif
  lf_modifier_destroy(m);

  /* --- 5. vignetting over the whole frame ------------------------------- */

  float *rgba = (float *)calloc((size_t)W * 4, sizeof(float));
  if(!rgba) return 1;
  for(int i = 0; i < W * 4; i++) rgba[i] = 0.5f;

  m = lf_modifier_new(lf, crop, W, H);
  lf_modifier_initialize(m, lf, LF_PF_F32, focal, aperture, distance, 1.f, LF_RECTILINEAR,
                         LF_MODIFY_VIGNETTING, 0);
  double t_lf_vig;
  TIME_BEST(t_lf_vig,
            for(int y = 0; y < H; y++)
              lf_modifier_apply_color_modification(m, rgba, 0.f, (float)y, W, 1,
                                                   LF_CR_4(RED, GREEN, BLUE, UNKNOWN),
                                                   W * 4 * (int)sizeof(float)));
  lf_modifier_destroy(m);

  ls_modifier_init(&mod, &lens, crop, W, H, focal, aperture, distance, 1.f, LS_LENS_UNKNOWN,
                   LS_ENABLE_VIGNETTING, 0);
  double t_ls_vig;
  TIME_BEST(t_ls_vig,
            for(int y = 0; y < H; y++)
              ls_modifier_apply_vignetting(&mod, 0.f, (float)y, W, 1, rgba,
                                           W * 4 * (int)sizeof(float)));
  row("apply vignetting, whole frame", t_lf_vig, t_ls_vig, "ms");

  double t_lf_vig_mt = 0.0, t_ls_vig_mt = 0.0;
#ifdef _OPENMP
  {
    float *big = (float *)malloc((size_t)W * H * 4 * sizeof(float));
    if(big)
    {
      for(size_t i = 0; i < (size_t)W * H * 4; i++) big[i] = 0.5f;
      lfModifier *vm = lf_modifier_new(lf, crop, W, H);
      lf_modifier_initialize(vm, lf, LF_PF_F32, focal, aperture, distance, 1.f, LF_RECTILINEAR,
                             LF_MODIFY_VIGNETTING, 0);
      TIME_BEST(t_lf_vig_mt,
                _Pragma("omp parallel for schedule(static)")
                for(int y = 0; y < H; y++)
                  lf_modifier_apply_color_modification(vm, big + (size_t)y * W * 4, 0.f,
                                                       (float)y, W, 1,
                                                       LF_CR_4(RED, GREEN, BLUE, UNKNOWN),
                                                       W * 4 * (int)sizeof(float)));
      lf_modifier_destroy(vm);

      ls_modifier_init(&mod, &lens, crop, W, H, focal, aperture, distance, 1.f, LS_LENS_UNKNOWN,
                       LS_ENABLE_VIGNETTING, 0);
      TIME_BEST(t_ls_vig_mt,
                _Pragma("omp parallel for schedule(static)")
                for(int y = 0; y < H; y++)
                  ls_modifier_apply_vignetting(&mod, 0.f, (float)y, W, 1,
                                               big + (size_t)y * W * 4,
                                               W * 4 * (int)sizeof(float)));
      free(big);
    }
  }
#endif

#ifdef _OPENMP
  if(t_lf_map_mt > 0.0)
  {
    printf("\n  the same two stages across %d OpenMP threads\n", omp_get_max_threads());
    row("build the whole geometry map", t_lf_map_mt, t_ls_map_mt, "ms");
    row("apply vignetting, whole frame", t_lf_vig_mt, t_ls_vig_mt, "ms");
    printf("  %-34s %10.2fx %9.2fx\n", "  scaling vs. its own 1 thread",
           t_lf_map / t_lf_map_mt, t_ls_map / t_ls_map_mt);
  }
#else
  (void)t_lf_map_mt; (void)t_ls_map_mt; (void)t_lf_vig_mt; (void)t_ls_vig_mt;
#endif

  /* --- 6. the map, versus never building one ----------------------------- */
  /*
   * The point of the whole library, measured on the CPU this time.
   *
   * A correction is only ever wanted so that something can RESAMPLE with it. lensfun can
   * only deliver it as a buffer: its callbacks fill six floats per pixel, and the consumer
   * then reads them back. LensSerious can be called per pixel, so the consumer can evaluate
   * the coordinates where it needs them and never materialise anything -- which is exactly
   * what the OpenCL kernels do, and there is no reason the CPU cannot do the same.
   *
   * All three produce the same image. The difference is 1.1 GB of memory traffic: 549 MB
   * written by the map pass and 549 MB read back by the resampler.
   */
  {
    const size_t pix = (size_t)W * H;
    float *src = (float *)malloc(pix * 4 * sizeof(float));
    float *dst = (float *)malloc(pix * 4 * sizeof(float));
    float *map = (float *)malloc(pix * 6 * sizeof(float));
    if(src && dst && map)
    {
      for(size_t i = 0; i < pix * 4; i++) src[i] = 0.5f;

      lfModifier *fm = lf_modifier_new(lf, crop, W, H);
      lf_modifier_initialize(fm, lf, LF_PF_F32, focal, aperture, distance, 1.f,
                             LF_RECTILINEAR, LF_MODIFY_DISTORTION | LF_MODIFY_TCA, 0);
      ls_modifier_init(&mod, &lens, crop, W, H, focal, aperture, distance, 1.f,
                       LS_LENS_UNKNOWN, LS_ENABLE_DISTORTION | LS_ENABLE_TCA, 0);
      ls_eval_t ep;
      ls_eval_from_modifier(&mod, &ep);

      double t_lf_two = 0.0, t_ls_two = 0.0, t_ls_fused = 0.0;

      /* lensfun: fill the map, then resample from it. */
      TIME_BEST(t_lf_two, {
        LS_OMP_FOR
        for(int y = 0; y < H; y++)
          lf_modifier_apply_subpixel_geometry_distortion(fm, 0.f, (float)y, W, 1,
                                                         map + (size_t)y * W * 6);
        LS_OMP_FOR
        for(int y = 0; y < H; y++)
          _resample_row(src, W, H, map + (size_t)y * W * 6, dst + (size_t)y * W * 4, W);
      });

      /* LensSerious, the same two passes, for a like-for-like comparison. */
      TIME_BEST(t_ls_two, {
        LS_OMP_FOR
        for(int y = 0; y < H; y++)
          ls_modifier_apply_subpixel_geometry(&mod, 0.f, (float)y, W, 1,
                                              map + (size_t)y * W * 6);
        LS_OMP_FOR
        for(int y = 0; y < H; y++)
          _resample_row(src, W, H, map + (size_t)y * W * 6, dst + (size_t)y * W * 4, W);
      });

      /* LensSerious, fused: the coordinates never leave a register. */
      TIME_BEST(t_ls_fused, {
        LS_OMP_FOR
        for(int y = 0; y < H; y++)
        {
          float *drow = dst + (size_t)y * W * 4;
          for(int x = 0; x < W; x++)
          {
            float c6[6];
            ls_eval_map(&ep, (float)x, (float)y, c6);
            _resample_pixel(src, W, H, c6, drow + (size_t)x * 4);
          }
        }
      });

      printf("\n  correct-and-resample a whole frame, %d thread(s)\n", LS_THREADS);
      printf("  %-34s %10.2f ms\n", "lensfun: map, then resample", t_lf_two);
      printf("  %-34s %10.2f ms\n", "LensSerious: map, then resample", t_ls_two);
      printf("  %-34s %10.2f ms   %8.1fx\n", "LensSerious: fused, no map",
             t_ls_fused, t_lf_two / t_ls_fused);
      printf("  the fused pass never writes the %zu MB map, nor reads it back\n",
             pix * 6 * sizeof(float) / 1048576);

      lf_modifier_destroy(fm);
    }
    free(src);
    free(dst);
    free(map);
  }

  /* --- what a session actually pays ------------------------------------- */

  const double lf_session = t_lf_data + t_lf_map + t_lf_vig;
  const double ls_session = t_ls_data + t_ls_map + t_ls_vig;
  printf("\n");
  row("one image, cold start", lf_session, ls_session, "ms");
  printf("  %-34s %10.4f %10.4f ms   %8.1fx\n", "one more image, warm",
         t_lf_init + t_lf_map + t_lf_vig, t_ls_init + t_ls_map + t_ls_vig,
         (t_lf_init + t_lf_map + t_lf_vig) / (t_ls_init + t_ls_map + t_ls_vig));
  printf("\n  The map row is the one that matters for a pixel pipeline, and it is the one\n"
         "  where a GPU removes the cost entirely rather than reducing it: see Ansel's\n"
         "  lens_distort_*_ls kernels, which evaluate it per work-item from ~80 bytes.\n");

  free(buf);
  free(rgba);
  lf_free(hits);
  ls_db_close(db);
  lf_db_destroy(ldb);
  return 0;
}
