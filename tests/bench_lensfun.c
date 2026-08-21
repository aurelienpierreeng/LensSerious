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
 * Everything here is single-threaded and CPU-side on both sides, deliberately. The GPU
 * numbers live in Ansel, where the kernel that consumes the map can be measured end to
 * end; comparing a GPU kernel against a CPU callback would prove nothing about the design.
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

static double now_ms(void)
{
  struct timespec t;
  clock_gettime(CLOCK_MONOTONIC, &t);
  return t.tv_sec * 1e3 + t.tv_nsec / 1e6;
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
                     LS_ENABLE_DISTORTION | LS_ENABLE_TCA | LS_ENABLE_VIGNETTING);
  const double t_ls_init = (now_ms() - t) / INITS;
  row("resolve at focal/aperture", t_lf_init, t_ls_init, "ms");

  /* --- 4. the geometry map, the whole frame ----------------------------- */

  float *buf = (float *)malloc((size_t)W * 6 * sizeof(float));
  if(!buf) return 1;

  lfModifier *m = lf_modifier_new(lf, crop, W, H);
  lf_modifier_initialize(m, lf, LF_PF_F32, focal, aperture, distance, 1.f, LF_RECTILINEAR,
                         LF_MODIFY_DISTORTION | LF_MODIFY_TCA, 0);
  t = now_ms();
  for(int y = 0; y < H; y++) lf_modifier_apply_subpixel_geometry_distortion(m, 0.f, (float)y, W, 1, buf);
  const double t_lf_map = now_ms() - t;
  lf_modifier_destroy(m);

  ls_modifier_init(&mod, &lens, crop, W, H, focal, aperture, distance, 1.f, LS_LENS_UNKNOWN,
                   LS_ENABLE_DISTORTION | LS_ENABLE_TCA);
  t = now_ms();
  for(int y = 0; y < H; y++) ls_modifier_apply_subpixel_geometry(&mod, 0.f, (float)y, W, 1, buf);
  const double t_ls_map = now_ms() - t;
  row("build the whole geometry map", t_lf_map, t_ls_map, "ms");

  /* --- 5. vignetting over the whole frame ------------------------------- */

  float *rgba = (float *)calloc((size_t)W * 4, sizeof(float));
  if(!rgba) return 1;
  for(int i = 0; i < W * 4; i++) rgba[i] = 0.5f;

  m = lf_modifier_new(lf, crop, W, H);
  lf_modifier_initialize(m, lf, LF_PF_F32, focal, aperture, distance, 1.f, LF_RECTILINEAR,
                         LF_MODIFY_VIGNETTING, 0);
  t = now_ms();
  for(int y = 0; y < H; y++)
    lf_modifier_apply_color_modification(m, rgba, 0.f, (float)y, W, 1,
                                         LF_CR_4(RED, GREEN, BLUE, UNKNOWN), W * 4 * (int)sizeof(float));
  const double t_lf_vig = now_ms() - t;
  lf_modifier_destroy(m);

  ls_modifier_init(&mod, &lens, crop, W, H, focal, aperture, distance, 1.f, LS_LENS_UNKNOWN,
                   LS_ENABLE_VIGNETTING);
  t = now_ms();
  for(int y = 0; y < H; y++)
    ls_modifier_apply_vignetting(&mod, 0.f, (float)y, W, 1, rgba, W * 4 * (int)sizeof(float));
  const double t_ls_vig = now_ms() - t;
  row("apply vignetting, whole frame", t_lf_vig, t_ls_vig, "ms");

  /* --- what a session actually pays ------------------------------------- */

  const double lf_session = t_lf_load + t_lf_find + t_lf_init + t_lf_map + t_lf_vig;
  const double ls_session = t_ls_open + t_ls_find + t_ls_init + t_ls_map + t_ls_vig;
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
