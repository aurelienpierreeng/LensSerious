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
                     float scale, int flags)
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

  mod->geometry_unsupported = (lens->type != LS_LENS_RECTILINEAR && lens->type != LS_LENS_UNKNOWN);

  int enabled = 0;
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
/* The map. Composition order is lensfun's callback-priority order for the     */
/* correction direction: scaling (100) → [projection: unsupported, reported]  */
/* → distortion (750), then the TCA subpixel stage per channel.               */
/* ------------------------------------------------------------------------- */

static inline void _dist_eval(const ls_calib_dist_t *d, float *x, float *y)
{
  const float xu = *x, yu = *y;
  const float ru2 = xu * xu + yu * yu;
  float m = 1.f;
  switch(d->model)
  {
    case LS_DIST_POLY3:
    {   /* mod-coord.cpp: Rd = Ru · (1 − k1 + k1·Ru²) */
      const float k1 = d->terms[0];
      m = (1.f - k1) + k1 * ru2;
      break;
    }
    case LS_DIST_POLY5:
    {   /* Rd = Ru · (1 + k1·Ru² + k2·Ru⁴) */
      m = 1.f + d->terms[0] * ru2 + d->terms[1] * ru2 * ru2;
      break;
    }
    case LS_DIST_PTLENS:
    {   /* Rd = Ru · (a·Ru³ + b·Ru² + c·Ru + d), d = 1−a−b−c */
      const float a = d->terms[0], b = d->terms[1], c = d->terms[2];
      const float r = sqrtf(ru2);
      m = a * ru2 * r + b * ru2 + c * r + (1.f - a - b - c);
      break;
    }
    default: break;
  }
  *x = xu * m;
  *y = yu * m;
}

static inline void _tca_eval(const ls_calib_tca_t *t, float *xr, float *yr, float *xb, float *yb)
{
  switch(t->model)
  {
    case LS_TCA_LINEAR:
    {   /* mod-subpix.cpp: per-channel radial scale */
      const float kr = t->terms[0], kb = t->terms[1];
      *xr *= kr; *yr *= kr;
      *xb *= kb; *yb *= kb;
      break;
    }
    case LS_TCA_POLY3:
    {   /* Rd = Ru · (b·Ru² + c·Ru + v), params packed vr vb cr cb br bb */
      const float vr = t->terms[0], vb = t->terms[1];
      const float cr = t->terms[2], cb = t->terms[3];
      const float br = t->terms[4], bb = t->terms[5];
      float x = *xr, y = *yr;
      float ru2 = x * x + y * y;
      float m = br * ru2 + vr + (cr != 0.f ? cr * sqrtf(ru2) : 0.f);
      *xr = x * m; *yr = y * m;
      x = *xb; y = *yb;
      ru2 = x * x + y * y;
      m = bb * ru2 + vb + (cb != 0.f ? cb * sqrtf(ru2) : 0.f);
      *xb = x * m; *yb = y * m;
      break;
    }
    default: break;
  }
}

int ls_modifier_apply_subpixel_geometry(const ls_modifier_t *mod,
                                        float xu, float yu, int width, int height,
                                        float *res)
{
  if(!mod || !res || width <= 0 || height <= 0) return 0;
  const int do_dist = (mod->enabled & LS_ENABLE_DISTORTION);
  const int do_tca = (mod->enabled & LS_ENABLE_TCA);
  const int do_scale = (mod->enabled & LS_ENABLE_SCALE);

  /* Exact per-column evaluation, identical to the OpenCL kernel, and deliberately NOT
   * bit-matched to upstream's row walker: lensfun's SSE path computes sqrt as
   * _mm_rcp_ps(_mm_rsqrt_ps(r2)) -- two chained 12-bit approximations with no Newton
   * step (mod-coord-sse.cpp) -- and disagrees with its own scalar math by up to 0.19 px
   * at the end of a 6016-wide row. Matching that would mean reimplementing an
   * approximation error. LensSerious matches upstream's SCALAR semantics to < 0.01 px;
   * against upstream's SSE rows the residual is bounded by upstream's own approximation,
   * and the parity harness asserts both bounds separately. */
  for(int row = 0; row < height; row++)
  {
    float *out = res + (size_t)row * width * 6;
    const float y0 = (yu + row) * mod->norm_scale - mod->center_y;
    for(int col = 0; col < width; col++, out += 6)
    {
      float x = (xu + col) * mod->norm_scale - mod->center_x;
      float y = y0;

      if(do_scale) { x *= mod->scale; y *= mod->scale; }
      if(do_dist) _dist_eval(&mod->dist, &x, &y);

      float xr = x, yr = y, xb = x, yb = y;
      if(do_tca) _tca_eval(&mod->tca, &xr, &yr, &xb, &yb);

      out[0] = (xr + mod->center_x) * mod->norm_unscale;
      out[1] = (yr + mod->center_y) * mod->norm_unscale;
      out[2] = (x  + mod->center_x) * mod->norm_unscale;
      out[3] = (y  + mod->center_y) * mod->norm_unscale;
      out[4] = (xb + mod->center_x) * mod->norm_unscale;
      out[5] = (yb + mod->center_y) * mod->norm_unscale;
    }
  }
  return 1;
}

int ls_modifier_apply_vignetting(const ls_modifier_t *mod,
                                 float xu, float yu, int width, int height,
                                 float *rgba, int row_stride_bytes)
{
  if(!mod || !rgba || width <= 0 || height <= 0) return 0;
  if(!(mod->enabled & LS_ENABLE_VIGNETTING)) return 0;
  const float k1 = mod->vig.terms[0], k2 = mod->vig.terms[1], k3 = mod->vig.terms[2];
  const size_t stride = row_stride_bytes ? (size_t)row_stride_bytes / sizeof(float)
                                         : (size_t)width * 4;

  /* "Damn! Hugin uses two different 'normalized' coordinate systems: for distortions it
   * uses 1.0 = min(half width, half height) and for vignetting it uses 1.0 = half
   * diagonal length." (mod-color.cpp) -- so this radius is the geometry one divided by
   * the aspect-ratio correction. And CORRECTION divides: ModifyColor_DeVignetting_PA
   * applies 1/c; multiplying by c is what APPLIES vignetting. The first harness run got
   * both of these wrong and measured multiplier deltas up to 37.9 for it. */
  const float arc = (mod->aspect_ratio_correction > 0.f) ? mod->aspect_ratio_correction : 1.f;

  for(int row = 0; row < height; row++)
  {
    float *px = rgba + (size_t)row * stride;
    const float y = ((yu + row) * mod->norm_scale - mod->center_y) / arc;
    for(int col = 0; col < width; col++, px += 4)
    {
      const float x = ((xu + col) * mod->norm_scale - mod->center_x) / arc;
      const float r2 = x * x + y * y;
      const float r4 = r2 * r2;
      const float c = 1.f + k1 * r2 + k2 * r4 + k3 * r4 * r2;
      const float inv = (c != 0.f) ? 1.f / c : 1.f;
      /* All FOUR components, alpha included: upstream's apply_multiplier under
       * LF_CR_4(RED,GREEN,BLUE,UNKNOWN) multiplies the UNKNOWN channel too, and the
       * harness caught the difference on literally every sampled pixel (k=3,7,11...).
       * Whether that is wise is not this library's question; parity is. */
      px[0] *= inv; px[1] *= inv; px[2] *= inv; px[3] *= inv;
    }
  }
  return 1;
}
