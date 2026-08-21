/*
    LensSerious — implementation.

    Copyright (C) 2026 Aurélien PIERRE.
    Ported from Lensfun 0.3.4, Copyright (C) 2007 Andrew Zabolotny and contributors.
    License: LGPL-3.0-or-later (see LICENSE).

    PORTING RULE, load-bearing: every numeric convention below reproduces lensfun 0.3.4
    behaviour ON PURPOSE, including the ones that look like bugs. The (width−1) sizing,
    the one-sided Hermite tangents, evaluating interpolation weights in the coefficient
    domain — a deviation here is not a cleanup, it silently moves every pixel of every
    corrected image. tests/parity_lensfun.c is the arbiter; change nothing the harness
    cannot re-verify.
*/

#include "lensserious.h"
#include "lensserious_eval.h"

#if defined(__GNUC__) || defined(__clang__)
  #define LS_RESTRICT __restrict__
#else
  #define LS_RESTRICT
#endif

#include <float.h>
#include <math.h>
#include <string.h>

/* ------------------------------------------------------------------------- */
/* Interpolation — ported from auxfun.cpp:_lf_interpolate() and               */
/* lens.cpp:__insert_spline()/InterpolateDistortion()/InterpolateTCA().       */
/* ------------------------------------------------------------------------- */

/* Hermite with one-sided tangents at the ends: y1/y4 == FLT_MAX means "no outer
 * neighbour", and the tangent degrades to the chord. Upstream auxfun.cpp:441. */
static float _interpolate(float y1, float y2, float y3, float y4, float t)
{
  const float t2 = t * t;
  const float t3 = t2 * t;
  const float tg2 = (y1 == FLT_MAX) ? y3 - y2 : (y3 - y1) * 0.5f;
  const float tg3 = (y4 == FLT_MAX) ? y3 - y2 : (y4 - y2) * 0.5f;
  return (2.f * t3 - 3.f * t2 + 1.f) * y2 + (t3 - 2.f * t2 + t) * tg2
       + (-2.f * t3 + 3.f * t2) * y3 + (t3 - t2) * tg3;
}

/* Keep the two nearest calibrations on each side of the target, by signed distance.
 * Upstream lens.cpp:802. Indices: [0] far-below, [1] near-below, [2] near-above,
 * [3] far-above. */
typedef struct spline_t { const void *v[4]; float d[4]; } spline_t;

static void _spline_init(spline_t *s)
{
  memset(s->v, 0, sizeof(s->v));
  s->d[0] = s->d[1] = -FLT_MAX;
  s->d[2] = s->d[3] = FLT_MAX;
}

static void _spline_insert(spline_t *s, float dist, const void *val)
{
  if(dist < 0)
  {
    if(dist > s->d[1]) { s->d[0] = s->d[1]; s->d[1] = dist; s->v[0] = s->v[1]; s->v[1] = val; }
    else if(dist > s->d[0]) { s->d[0] = dist; s->v[0] = val; }
  }
  else
  {
    if(dist < s->d[2]) { s->d[3] = s->d[2]; s->d[2] = dist; s->v[3] = s->v[2]; s->v[2] = val; }
    else if(dist < s->d[3]) { s->d[3] = dist; s->v[3] = val; }
  }
}

/* Upstream interpolates COEFFICIENTS across focal, with per-term scaling hooks
 * (__parameter_scales). For every model in the shipping database those scales are
 * identity for distortion and index<2 of TCA uses focal-domain scaling — ported below
 * exactly as the 0.3.4 code has it (lens.cpp:841): distortion: none; TCA terms 0..1:
 * none (the switch body falls through empty for LINEAR/POLY3 v-terms in 0.3.4). */
static int _interp_dist(const ls_lens_t *lens, float focal, ls_calib_dist_t *res)
{
  if(lens->n_dist == 0) return 0;
  spline_t s; _spline_init(&s);
  ls_dist_model_t model = LS_DIST_NONE;

  for(int i = 0; i < lens->n_dist; i++)
  {
    const ls_calib_dist_t *c = &lens->dist[i];
    if(c->model == LS_DIST_NONE) continue;
    if(model == LS_DIST_NONE) model = c->model;
    else if(model != c->model) continue;           /* first model wins, upstream warning case */
    const float df = focal - c->focal;
    if(df == 0.0f) { *res = *c; return 1; }
    _spline_insert(&s, df, c);
  }
  const ls_calib_dist_t *lo = (const ls_calib_dist_t *)s.v[1];
  const ls_calib_dist_t *hi = (const ls_calib_dist_t *)s.v[2];
  if(!lo || !hi)
  {
    if(lo) { *res = *lo; return 1; }
    if(hi) { *res = *hi; return 1; }
    return 0;
  }
  res->model = model;
  res->focal = focal;
  const float t = (focal - lo->focal) / (hi->focal - lo->focal);
  const ls_calib_dist_t *fl = (const ls_calib_dist_t *)s.v[0];
  const ls_calib_dist_t *fh = (const ls_calib_dist_t *)s.v[3];
  /* Upstream __parameter_scales (lens.cpp:841) leaves the FOCALS in for every distortion
   * model: coefficients are interpolated in the term×focal domain, then divided by the
   * target focal. The first harness run against the raw-coefficient version measured up
   * to 4.9 px of divergence at interpolated focals -- this scaling is load-bearing. */
  for(int i = 0; i < 3; i++)
    res->terms[i] = _interpolate(fl ? fl->terms[i] * fl->focal : FLT_MAX,
                                 lo->terms[i] * lo->focal, hi->terms[i] * hi->focal,
                                 fh ? fh->terms[i] * fh->focal : FLT_MAX, t) / focal;
  return 1;
}

static int _interp_tca(const ls_lens_t *lens, float focal, ls_calib_tca_t *res)
{
  if(lens->n_tca == 0) return 0;
  spline_t s; _spline_init(&s);
  ls_tca_model_t model = LS_TCA_NONE;

  for(int i = 0; i < lens->n_tca; i++)
  {
    const ls_calib_tca_t *c = &lens->tca[i];
    if(c->model == LS_TCA_NONE) continue;
    if(model == LS_TCA_NONE) model = c->model;
    else if(model != c->model) continue;
    const float df = focal - c->focal;
    if(df == 0.0f) { *res = *c; return 1; }
    _spline_insert(&s, df, c);
  }
  const ls_calib_tca_t *lo = (const ls_calib_tca_t *)s.v[1];
  const ls_calib_tca_t *hi = (const ls_calib_tca_t *)s.v[2];
  if(!lo || !hi)
  {
    if(lo) { *res = *lo; return 1; }
    if(hi) { *res = *hi; return 1; }
    return 0;
  }
  res->model = model;
  res->focal = focal;
  const float t = (focal - lo->focal) / (hi->focal - lo->focal);
  const ls_calib_tca_t *fl = (const ls_calib_tca_t *)s.v[0];
  const ls_calib_tca_t *fh = (const ls_calib_tca_t *)s.v[3];
  /* Same term×focal domain as distortion, EXCEPT terms 0..1 (vr/vb, kr/kb), which
   * __parameter_scales exempts by setting the scales to 1.0. */
  for(int i = 0; i < 6; i++)
  {
    const float sl0 = (i < 2) ? 1.f : (fl ? fl->focal : 1.f);
    const float sl1 = (i < 2) ? 1.f : lo->focal;
    const float sl2 = (i < 2) ? 1.f : hi->focal;
    const float sl3 = (i < 2) ? 1.f : (fh ? fh->focal : 1.f);
    const float sl4 = (i < 2) ? 1.f : focal;
    res->terms[i] = _interpolate(fl ? fl->terms[i] * sl0 : FLT_MAX,
                                 lo->terms[i] * sl1, hi->terms[i] * sl2,
                                 fh ? fh->terms[i] * sl3 : FLT_MAX, t) / sl4;
  }
  return 1;
}

/* Vignetting — verbatim port of lfLens::InterpolateVignetting() and __vignetting_dist()
 * (lens.cpp): inverse-distance weighting with p = 3.5 over ALL calibration points, in a
 * space where focal is normalized by the lens's whole range, aperture enters as 4/A and
 * distance as 0.1/D. An exact hit (< 1e-4) short-circuits; a nearest point further than
 * 1.0 rejects the whole interpolation. The first version of this file used an ad-hoc
 * metric here; the harness measured multiplier deltas up to 3.65 against upstream and
 * this port is what removed them. */
static float _vig_dist(const ls_lens_t *lens, const ls_calib_vig_t *c,
                       float focal, float aperture, float distance)
{
  float f1 = focal - lens->min_focal;
  float f2 = c->focal - lens->min_focal;
  const float df = lens->max_focal - lens->min_focal;
  if(df != 0.f) { f1 /= df; f2 /= df; }
  const float a1 = 4.f / aperture;
  const float a2 = 4.f / c->aperture;
  const float d1 = 0.1f / distance;
  const float d2 = 0.1f / c->distance;
  return sqrtf((f2 - f1) * (f2 - f1) + (a2 - a1) * (a2 - a1) + (d2 - d1) * (d2 - d1));
}

static int _interp_vig(const ls_lens_t *lens, float focal, float aperture, float distance,
                       ls_calib_vig_t *res)
{
  if(lens->n_vig == 0) return 0;
  ls_vig_model_t model = LS_VIG_NONE;
  float total_weighting = 0.f;
  float smallest = FLT_MAX;
  float terms[3] = { 0.f, 0.f, 0.f };

  for(int i = 0; i < lens->n_vig; i++)
  {
    const ls_calib_vig_t *c = &lens->vig[i];
    if(c->model == LS_VIG_NONE) continue;
    if(model == LS_VIG_NONE) model = c->model;
    else if(model != c->model) continue;

    const float dist = _vig_dist(lens, c, focal, aperture, distance);
    if(dist < 0.0001f) { *res = *c; return 1; }
    if(dist < smallest) smallest = dist;
    /* Upstream computes the weight in DOUBLE (fabs(1.0 / pow(dist, 3.5))) and only then
     * truncates: near an exact hit d^3.5 ~ 1e-14 and the two arithmetics disagree enough
     * to change the mixture -- the harness measured 0.77 of multiplier drift for it. */
    const float w = (float)fabs(1.0 / pow((double)dist, 3.5));
    for(int k = 0; k < 3; k++) terms[k] += w * c->terms[k];
    total_weighting += w;
  }
  if(smallest > 1.f) return 0;
  if(total_weighting <= 0.f || smallest == FLT_MAX) return 0;

  res->model = model;
  res->focal = focal; res->aperture = aperture; res->distance = distance;
  for(int k = 0; k < 3; k++) res->terms[k] = terms[k] / total_weighting;
  return 1;
}

/* ------------------------------------------------------------------------- */
/* Coordinate system — ported from lfModifier::lfModifier() (modifier.cpp:203) */
/* ------------------------------------------------------------------------- */

int ls_modifier_init(ls_modifier_t *mod, const ls_lens_t *lens,
                     float crop, int width, int height,
                     float focal, float aperture, float distance,
                     float scale, int target_type, int flags)
{
  memset(mod, 0, sizeof(*mod));
  if(!lens || crop <= 0.f) return 0;

  /* "The '- 1' is due to the fact that Width and Height are measured at the pixel
   * centres (they are actually transformed) instead of at their outer rims." */
  const float w = (width  >= 2) ? (float)(width  - 1) : 1.f;
  const float h = (height >= 2) ? (float)(height - 1) : 1.f;
  mod->width = w; mod->height = h;

  const float size = (w < h) ? w : h;
  const float image_aspect_ratio = (w < h) ? h / w : w / h;

  const float calibration_cropfactor = lens->crop_factor;
  const float ar = (lens->aspect_ratio > 0.f) ? lens->aspect_ratio : 1.5f;
  const float aspect_ratio_correction = sqrtf(ar * ar + 1.f);

  const float coordinate_correction =
      1.f / sqrtf(image_aspect_ratio * image_aspect_ratio + 1.f)
      * calibration_cropfactor / crop * aspect_ratio_correction;

  mod->norm_scale = 2.f / size * coordinate_correction;
  mod->norm_unscale = size * 0.5f / coordinate_correction;
  mod->aspect_ratio_correction = aspect_ratio_correction;
  mod->center_x = (w / size + lens->center_x) * coordinate_correction;
  mod->center_y = (h / size + lens->center_y) * coordinate_correction;

  /* A projection change is radial for every type except panoramic and equirectangular,
   * which treat the axes differently (see ls_eval_geometry()). Report those rather than
   * approximate them. */
  const int from = (int)lens->type;
  const int to = (target_type == LS_LENS_UNKNOWN) ? from : target_type;
  const int radial = (from != LS_LENS_PANORAMIC && from != LS_LENS_EQUIRECTANGULAR
                      && to != LS_LENS_PANORAMIC && to != LS_LENS_EQUIRECTANGULAR
                      && from != LS_LENS_UNKNOWN && to != LS_LENS_UNKNOWN);
  mod->geom_from = from;
  mod->geom_to = to;
  /* The projection focal, in the same normalized units as a radius.
   *
   * MEASURED, not derived: liblensfun's own geometry callback was fitted over lenses
   * spanning crop factors 1.0 to 7.66 and aspect ratios 1.143 to 1.5, and
   *
   *     f_norm = focal * lens_crop_factor * sqrt(ar^2 + 1) / 21.633307
   *
   * reproduces all of them to five digits. The denominator is half the diagonal of a
   * 36x24 frame, so the numerator is that diagonal reduced to the CALIBRATION sensor and
   * then to its short side -- which is what normalized radius 1.0 means here.
   *
   * The SHOOTING crop deliberately does not appear: it is already inside norm_scale, and
   * a first version that used a bare focal/12 (fitted on 3:2 lenses alone, where the two
   * happen to coincide) was out by 284 px on a 4:3 compact. */
  mod->geom_focal = focal * lens->crop_factor * aspect_ratio_correction
                    / LS_EVAL_FULL_FRAME_HALF_DIAG_MM;
  mod->geometry_unsupported = (from != to) && !radial;

  /* A projection change is offered where the geometry FOCAL is known to be the nominal
   * one, which is the whole database except the handful of lenses carrying
   * <real-focal-length> data.
   *
   * lensfun's geometry callback runs on GetRealFocalLength(focal) divided by
   * get_hugin_focal_correction(focal) (modifier.cpp). Without real-focal data the first
   * multiplies the hugin factor in and the second divides it straight back out, so the
   * result is exactly `focal` -- measured across every fisheye in the database, ratio
   * 1.0000. With real-focal data they do not cancel and the focal becomes
   * real_focal / hugin, which ls_lens_t does not carry the calibration points to compute.
   *
   * Getting that focal wrong is not subtle: on the Sigma 4.5mm circular fisheye, whose
   * real focal is 4.533 against a nominal 4.5 and whose hugin factor is 2.13, using the
   * nominal focal is 28 px out at the centre of the frame and 135 px at the edge. So those
   * lenses fall back instead.
   *
   * The composition ORDER is settled and is the one below: lensfun's correction chain is
   * scale (100), geometry (500), distortion (750, ModifyCoord_Dist_* -- the 250 UnDist_*
   * variants are the reverse direction), then TCA as a subpixel callback after every
   * coordinate callback. An earlier reading of this file's own header suggested the
   * reverse and was wrong. */
  if(from != to && radial && lens->has_real_focal) mod->geometry_unsupported = 1;

  int enabled = 0;
  if(from != to && radial && !lens->has_real_focal && focal > 0.f) enabled |= LS_ENABLE_GEOMETRY;
  if((flags & LS_ENABLE_DISTORTION) && _interp_dist(lens, focal, &mod->dist))
    enabled |= LS_ENABLE_DISTORTION;
  if((flags & LS_ENABLE_TCA) && _interp_tca(lens, focal, &mod->tca))
    enabled |= LS_ENABLE_TCA;
  if((flags & LS_ENABLE_VIGNETTING) && _interp_vig(lens, focal, aperture, distance, &mod->vig))
    enabled |= LS_ENABLE_VIGNETTING;
  if((flags & LS_ENABLE_SCALE) && scale != 1.0f && scale > 0.f)
  {
    mod->scale = 1.f / scale;   /* upstream stores the reciprocal for the correction pass */
    enabled |= LS_ENABLE_SCALE;
  }
  else
    mod->scale = 1.f;

  mod->enabled = enabled;
  return enabled;
}

/* ------------------------------------------------------------------------- */
/* The map. The closed forms themselves live in include/lensserious_eval.h,   */
/* which the OpenCL kernel includes verbatim: the CPU and the GPU evaluate    */
/* one source text, so they cannot drift apart between releases.              */
/* ------------------------------------------------------------------------- */

/* ls_eval_from_modifier() casts the model enums straight to int, and the kernel compares
 * the result against the LS_EVAL_* mirrors in lensserious_eval.h -- which cannot name the
 * enums, since it must also compile as OpenCL C. Renumbering an enum without touching its
 * mirror would make every GPU render evaluate the wrong model, silently and only on the
 * GPU. Cheapest possible place to catch that is here, at compile time. */
_Static_assert((int)LS_DIST_NONE   == LS_EVAL_DIST_NONE,   "distortion model mirror drifted");
_Static_assert((int)LS_DIST_POLY3  == LS_EVAL_DIST_POLY3,  "distortion model mirror drifted");
_Static_assert((int)LS_DIST_POLY5  == LS_EVAL_DIST_POLY5,  "distortion model mirror drifted");
_Static_assert((int)LS_DIST_PTLENS == LS_EVAL_DIST_PTLENS, "distortion model mirror drifted");
_Static_assert((int)LS_TCA_NONE    == LS_EVAL_TCA_NONE,    "TCA model mirror drifted");
_Static_assert((int)LS_TCA_LINEAR  == LS_EVAL_TCA_LINEAR,  "TCA model mirror drifted");
_Static_assert((int)LS_TCA_POLY3   == LS_EVAL_TCA_POLY3,   "TCA model mirror drifted");
_Static_assert((int)LS_VIG_NONE    == LS_EVAL_VIG_NONE,    "vignetting model mirror drifted");
_Static_assert((int)LS_VIG_PA      == LS_EVAL_VIG_PA,      "vignetting model mirror drifted");

_Static_assert(LS_ENABLE_DISTORTION == LS_EVAL_ENABLE_DISTORTION, "enable bit mirror drifted");
_Static_assert(LS_ENABLE_TCA        == LS_EVAL_ENABLE_TCA,        "enable bit mirror drifted");
_Static_assert(LS_ENABLE_VIGNETTING == LS_EVAL_ENABLE_VIGNETTING, "enable bit mirror drifted");
_Static_assert(LS_ENABLE_SCALE      == LS_EVAL_ENABLE_SCALE,      "enable bit mirror drifted");
_Static_assert(LS_ENABLE_GEOMETRY   == LS_EVAL_ENABLE_GEOMETRY,   "enable bit mirror drifted");
_Static_assert((int)LS_LENS_RECTILINEAR == LS_EVAL_LENS_RECTILINEAR, "lens type mirror drifted");
_Static_assert((int)LS_LENS_FISHEYE_THOBY == LS_EVAL_LENS_FISHEYE_THOBY, "lens type mirror drifted");

/* ls_eval_t crosses to the device as a by-value kernel argument, so host and device must
 * lay it out identically. Scalars only is what guarantees that; this pins the consequence
 * so adding a float2 (or a double, or a bool) fails to build instead of corrupting every
 * field after it. */
_Static_assert(sizeof(ls_eval_t) == 8 * sizeof(float)      /* the coordinate system */
                                  + 4 * sizeof(int)        /* model ids and enable bits */
                                  + 2 * sizeof(int)        /* geom_from, geom_to */
                                  + 2 * sizeof(float)      /* geom_focal, _pad */
                                  + 12 * sizeof(float),    /* the terms */
               "ls_eval_t gained padding or a member: check it is still scalar-only");
_Static_assert(_Alignof(ls_eval_t) == _Alignof(float), "ls_eval_t alignment is no longer 4");

int ls_eval_from_modifier(const ls_modifier_t *mod, ls_eval_t *out)
{
  if(!mod || !out) return 0;
  memset(out, 0, sizeof(*out));

  out->norm_scale = mod->norm_scale;
  out->norm_unscale = mod->norm_unscale;
  out->center_x = mod->center_x;
  out->center_y = mod->center_y;
  {
    const float arc = (mod->aspect_ratio_correction > 0.f) ? mod->aspect_ratio_correction : 1.f;
    const float inv_arc = 1.f / arc;
    out->vig_scale = mod->norm_scale * inv_arc;
    out->vig_center_x = mod->center_x * inv_arc;
    out->vig_center_y = mod->center_y * inv_arc;
  }
  out->scale = mod->scale;
  out->enabled = mod->enabled;

  out->dist_model = (int)mod->dist.model;
  out->tca_model = (int)mod->tca.model;
  out->vig_model = (int)mod->vig.model;
  out->geom_from = mod->geom_from;
  out->geom_to = mod->geom_to;
  out->geom_focal = mod->geom_focal;

  for(int i = 0; i < 3; i++) out->dist_terms[i] = mod->dist.terms[i];
  for(int i = 0; i < 6; i++) out->tca_terms[i] = mod->tca.terms[i];
  for(int i = 0; i < 3; i++) out->vig_terms[i] = mod->vig.terms[i];

  return 1;
}

int ls_modifier_apply_subpixel_geometry(const ls_modifier_t *mod,
                                        float xu, float yu, int width, int height,
                                        float *res)
{
  if(!mod || !res || width <= 0 || height <= 0) return 0;

  ls_eval_t p;
  ls_eval_from_modifier(mod, &p);

  /* Exact per-pixel evaluation, identical to the OpenCL kernel because it IS the OpenCL
   * kernel's code, and deliberately NOT bit-matched to upstream's row walker: lensfun's
   * SSE path computes sqrt as _mm_rcp_ps(_mm_rsqrt_ps(r2)) -- two chained 12-bit
   * approximations with no Newton step (mod-coord-sse.cpp) -- and disagrees with its own
   * scalar math by up to 0.19 px at the end of a 6016-wide row. Matching that would mean
   * reimplementing an approximation error. LensSerious matches upstream's SCALAR semantics
   * to < 0.01 px; against upstream's SSE rows the residual is bounded by upstream's own
   * approximation, and the parity harness asserts both bounds separately. */
  for(int row = 0; row < height; row++)
  {
    float *out = res + (size_t)row * width * 6;
    const float y = yu + (float)row;
    for(int col = 0; col < width; col++, out += 6)
      ls_eval_map(&p, xu + (float)col, y, out);
  }
  return 1;
}

int ls_modifier_apply_vignetting(const ls_modifier_t *mod,
                                 float xu, float yu, int width, int height,
                                 float *rgba, int row_stride_bytes)
{
  if(!mod || !rgba || width <= 0 || height <= 0) return 0;
  if(!(mod->enabled & LS_ENABLE_VIGNETTING)) return 0;

  ls_eval_t p;
  ls_eval_from_modifier(mod, &p);

  const size_t stride = row_stride_bytes ? (size_t)row_stride_bytes / sizeof(float)
                                         : (size_t)width * 4;
  const float vs = p.vig_scale;
  const float vcx = p.vig_center_x;

  /* Two passes over a block, and the reason is the divide.
   *
   * Written as one loop -- derive the multiplier, then scale the pixel -- the compiler
   * vectorises only the four-component store, because the store is what looks like a
   * vector. The 1/c stays scalar, one divide per pixel, and a single-precision divide is
   * the most expensive thing in this function by a wide margin.
   *
   * Computing the multipliers for a block FIRST gives a loop with no stores to the image
   * in it, which vectorises across PIXELS: four divides per instruction. The second pass
   * is then pure multiply-and-store. Measured on 24 Mpx, this is what took the function
   * from 1.8x slower than lensfun's hand-written SSE2 to parity with it.
   *
   * The block is sized to stay in L1 alongside the pixels it is about to scale. */
  enum { LS_VIG_BLOCK = 256 };
  float mbuf[LS_VIG_BLOCK];

  for(int row = 0; row < height; row++)
  {
    float *LS_RESTRICT px = rgba + (size_t)row * stride;

    /* Everything constant along the row is computed once. The x coordinate is still
     * derived from the absolute column rather than stepped incrementally, so this stays
     * bit-identical to what ls_eval_vignette_factor() produces in a kernel -- upstream
     * steps r2 by a recurrence, which is cheaper still and not reproducible per work-item. */
    const float y = (yu + (float)row) * vs - p.vig_center_y;
    const float yy = y * y;

    for(int col0 = 0; col0 < width; col0 += LS_VIG_BLOCK)
    {
      const int n = (width - col0 < LS_VIG_BLOCK) ? (width - col0) : LS_VIG_BLOCK;

      for(int i = 0; i < n; i++)
      {
        const float x = (xu + (float)(col0 + i)) * vs - vcx;
        mbuf[i] = ls_eval_vignette_from_r2(&p, x * x + yy);
      }

      for(int i = 0; i < n; i++, px += 4)
      {
        const float m = mbuf[i];
        /* All FOUR components, alpha included: upstream's apply_multiplier under
         * LF_CR_4(RED, GREEN, BLUE, UNKNOWN) multiplies the UNKNOWN channel too, and the
         * harness caught the difference on literally every sampled pixel (k=3,7,11...).
         * Whether that is wise is not this library's question; parity is. */
        px[0] *= m;
        px[1] *= m;
        px[2] *= m;
        px[3] *= m;
      }
    }
  }

  return 1;
}
