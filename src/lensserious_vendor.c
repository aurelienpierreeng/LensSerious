/*
    LensSerious — a maker's own correction tables, decoded into the library's generic knots.

    Copyright (C) 2026 Aurélien PIERRE.
    Ported from darktable's embedded-metadata lens correction (src/iop/lens.cc and
    src/common/dng_opcode.c), Copyright (C) the darktable contributors.

    License: GPL-3.0-or-later — see the note in lensserious_vendor.h.
*/

/** @file lensserious_vendor.c
 *
 * @brief The per-vendor decoders behind ls_vendor_resolve(), and the OpcodeList3 parser.
 *
 * @details Every decoder ends in the same place: radii and per-channel scale factors in
 * #ls_knots_t layout, half-diagonal normalized. What differs per vendor is the road there —
 * Sony's fixed-point tables at implicit radii, Fujifilm's percentage tables at explicit
 * radii needing a resample, a DNG's closed-form polynomials sampled onto knots, Olympus's
 * likewise. The maths is a faithful port of the darktable/Ansel implementation; a pixel
 * rendered through these knots must not move because the code changed repositories.
 */

#include "lensserious_vendor.h"

#include <math.h>
#include <string.h>

/* ------------------------------------------------------------------------- */
/* Shared helpers                                                             */
/* ------------------------------------------------------------------------- */

/** Piecewise-linear interpolation over (xi, yi), clamped at both ends. The resampler the
 * Fujifilm decoder and the autoscale sweep share. */
static float _linear_spline(const float *xi, const float *yi, const int ni, const float x)
{
  if(ni <= 0) return 1.0f;
  if(x < xi[0]) return yi[0];

  for(int i = 1; i < ni; i++)
  {
    if(x >= xi[i - 1] && x <= xi[i])
    {
      const float denom = xi[i] - xi[i - 1];
      if(denom == 0.0f) return yi[i - 1];
      const float dydx = (yi[i] - yi[i - 1]) / denom;
      return yi[i - 1] + (x - xi[i - 1]) * dydx;
    }
  }

  return yi[ni - 1];
}

static const ls_vendor_finetune_t _ft_identity = { 1.f, 1.f, 1.f, 1.f };

/* ------------------------------------------------------------------------- */
/* Sony                                                                       */
/* ------------------------------------------------------------------------- */

static int _sony_has_data(const ls_vendor_data_t *d)
{
  return d->u.sony.nc >= 2 && d->u.sony.nc <= LS_MAX_KNOTS;
}

static int _sony_axes(const ls_vendor_data_t *d)
{
  if(!_sony_has_data(d)) return 0;
  const ls_vendor_sony_t *const sony = &d->u.sony;

  int axes = 0;
  for(int i = 0; i < sony->nc; i++)
  {
    if(sony->distortion[i] != 0) axes |= LS_ENABLE_DISTORTION;
    if(sony->ca_r[i] != 0 || sony->ca_b[i] != 0) axes |= LS_ENABLE_TCA;
    if(sony->vignetting[i] != 0) axes |= LS_ENABLE_VIGNETTING;
  }
  return axes;
}

static int _sony_populate(const ls_vendor_data_t *d, const ls_vendor_finetune_t *ft,
                          ls_knots_t *knots)
{
  static const float SONY_DIST_SCALE = 1.0f / 16384.0f;
  static const float SONY_CA_SCALE = 1.0f / 2097152.0f;
  static const float SONY_VIG_SCALE = 1.0f / 8192.0f;

  const ls_vendor_sony_t *const sony = &d->u.sony;
  const int nc = sony->nc;

  for(int i = 0; i < nc; i++)
  {
    const float frac = (float)(i + 0.5) / (float)(nc - 1);
    knots->radius[i] = frac;
    knots->vig_radius[i] = frac;

    const float dist_cor = ft->distortion * ((float)sony->distortion[i] * SONY_DIST_SCALE) + 1.0f;
    knots->cor_rgb[0][i] = dist_cor;
    knots->cor_rgb[1][i] = dist_cor;
    knots->cor_rgb[2][i] = dist_cor;

    knots->cor_rgb[0][i] *= ft->ca_red * ((float)sony->ca_r[i] * SONY_CA_SCALE) + 1.0f;
    knots->cor_rgb[2][i] *= ft->ca_blue * ((float)sony->ca_b[i] * SONY_CA_SCALE) + 1.0f;

    const float val = ft->vignette * ((float)sony->vignetting[i] * SONY_VIG_SCALE);
    knots->vig[i] = powf(2.0f, 0.5f - powf(2.0f, val - 1.0f));
  }

  return nc;
}

/* ------------------------------------------------------------------------- */
/* Fujifilm                                                                   */
/* ------------------------------------------------------------------------- */

static int _fuji_has_data(const ls_vendor_data_t *d)
{
  return d->u.fuji.nc > 0 && d->u.fuji.nc <= 11;
}

static int _fuji_axes(const ls_vendor_data_t *d)
{
  if(!_fuji_has_data(d)) return 0;
  const ls_vendor_fuji_t *const fuji = &d->u.fuji;

  int axes = 0;
  for(int i = 0; i < fuji->nc; i++)
  {
    if(fuji->distortion[i] != 0.0f) axes |= LS_ENABLE_DISTORTION;
    if(fuji->ca_r[i] != 0.0f || fuji->ca_b[i] != 0.0f) axes |= LS_ENABLE_TCA;
    if(fuji->vignetting[i] != 0.0f) axes |= LS_ENABLE_VIGNETTING;
  }
  return axes;
}

/* Fujifilm publishes the tables at ITS knots, in percent, indexed by SOURCE radius; the
 * knot model wants them at OUR radii, as factors, indexed by DESTINATION radius. So:
 * resample onto a regular grid, convert each sample to a factor, then move the radius from
 * the source system to the destination one by dividing the correction back out. cropf is
 * applied to the maker's radii first — see the contract in lensserious_vendor.h. */
static int _fuji_populate(const ls_vendor_data_t *d, const ls_vendor_finetune_t *ft,
                          ls_knots_t *knots)
{
  const ls_vendor_fuji_t *const fuji = &d->u.fuji;
  const int ncsrc = fuji->nc;

  float knots_in[LS_MAX_KNOTS] = { 0.f };
  float cor_rgb_in[LS_MAX_KNOTS] = { 0.f };
  float cor_ca_r_in[LS_MAX_KNOTS] = { 0.f };
  float cor_ca_b_in[LS_MAX_KNOTS] = { 0.f };

  int j = 0;
  if(fuji->knots[0] > 0.0f)
  {
    knots_in[j] = 0.0f;
    cor_rgb_in[j] = 1.0f;
    cor_ca_r_in[j] = 0.0f;
    cor_ca_b_in[j] = 0.0f;
    knots->vig_radius[j] = 0.0f;
    knots->vig[j] = 1.0f;
    j++;
  }

  for(int i = 0; i < ncsrc; i++, j++)
  {
    knots_in[j] = fuji->cropf * fuji->knots[i];
    cor_rgb_in[j] = ft->distortion * (fuji->distortion[i] / 100.0f) + 1.0f;
    cor_ca_r_in[j] = ft->ca_red * fuji->ca_r[i];
    cor_ca_b_in[j] = ft->ca_blue * fuji->ca_b[i];

    knots->vig_radius[j] = fuji->cropf * fuji->knots[i];
    knots->vig[j] = 1.0f - ft->vignette * (1.0f - fuji->vignetting[i] / 100.0f);
  }
  const int ncin = j;
  if(ncin <= 0) return 0;

  for(int k = ncin; k < LS_MAX_KNOTS; k++)
  {
    knots->vig_radius[k] = knots->vig_radius[ncin - 1] + (float)(k - ncin + 1);
    knots->vig[k] = knots->vig[ncin - 1];
  }

  const int nc = LS_MAX_KNOTS;
  for(int i = 0; i < nc; i++)
  {
    const float rin = (float)i / (float)(nc - 1);
    const float m = _linear_spline(knots_in, cor_rgb_in, ncin, rin);
    const float r = (fabsf(m) > 1e-6f) ? rin / m : rin;
    knots->radius[i] = r;

    knots->cor_rgb[0][i] = m;
    knots->cor_rgb[1][i] = m;
    knots->cor_rgb[2][i] = m;

    const float mcar = _linear_spline(knots_in, cor_ca_r_in, ncin, rin);
    const float mcab = _linear_spline(knots_in, cor_ca_b_in, ncin, rin);
    knots->cor_rgb[0][i] *= mcar + 1.0f;
    knots->cor_rgb[2][i] *= mcab + 1.0f;
  }

  return nc;
}

/* ------------------------------------------------------------------------- */
/* DNG                                                                        */
/* ------------------------------------------------------------------------- */

static int _dng_axes(const ls_vendor_data_t *d)
{
  const ls_vendor_dng_t *const dng = &d->u.dng;
  int axes = 0;
  if(dng->has_warp) axes |= LS_ENABLE_DISTORTION;
  if(dng->has_warp && dng->warp_planes > 1) axes |= LS_ENABLE_TCA;
  if(dng->has_vignette) axes |= LS_ENABLE_VIGNETTING;
  return axes;
}

/** WarpRectilinear's radial part: kr0..kr3 evaluated at r^2, by Horner. */
static double _dng_warp_radial(const double coeffs[6], const double r2)
{
  return coeffs[0] + r2 * (coeffs[1] + r2 * (coeffs[2] + r2 * coeffs[3]));
}

static int _dng_populate(const ls_vendor_data_t *d, const ls_vendor_finetune_t *ft,
                         ls_knots_t *knots)
{
  const ls_vendor_dng_t *const dng = &d->u.dng;
  const int nc = LS_MAX_KNOTS;
  const uint32_t three = 3;
  const int nplanes = (int)(dng->warp_planes < three ? dng->warp_planes : three);
  const int canonical_plane = (dng->warp_planes > 1) ? 1 : 0;
  const int apply_tca = dng->warp_planes > 1;

  for(int i = 0; i < nc; i++)
  {
    const float r = (float)i / (float)(nc - 1);
    knots->radius[i] = r;
    knots->vig_radius[i] = r;
    const double r2 = (double)r * (double)r;

    if(dng->has_warp)
    {
      for(int c = 0; c < 3; c++)
      {
        const int plane = apply_tca ? (c < nplanes - 1 ? c : nplanes - 1) : canonical_plane;
        const double r_cor = _dng_warp_radial(dng->warp_coeffs[plane], r2);
        knots->cor_rgb[c][i] = (float)(ft->distortion * (r_cor - 1.0) + 1.0);
      }
    }
    else
    {
      knots->cor_rgb[0][i] = 1.0f;
      knots->cor_rgb[1][i] = 1.0f;
      knots->cor_rgb[2][i] = 1.0f;
    }

    if(dng->has_vignette)
    {
      const double dvig = r2
          * (dng->vig_coeffs[0]
             + r2 * (dng->vig_coeffs[1]
                     + r2 * (dng->vig_coeffs[2] + r2 * (dng->vig_coeffs[3] + r2 * dng->vig_coeffs[4]))));
      knots->vig[i] = (float)(1.0 / (1.0 + ft->vignette * dvig));
    }
    else
    {
      knots->vig[i] = 1.0f;
    }
  }

  return nc;
}

/* ------------------------------------------------------------------------- */
/* Olympus                                                                    */
/* ------------------------------------------------------------------------- */

static int _olympus_axes(const ls_vendor_data_t *d)
{
  const ls_vendor_olympus_t *const oly = &d->u.olympus;
  int axes = 0;
  if(oly->has_dist) axes |= LS_ENABLE_DISTORTION;
  if(oly->has_ca) axes |= LS_ENABLE_TCA;
  /* No vignetting, ever: Olympus does not publish one. */
  return axes;
}

static int _olympus_populate(const ls_vendor_data_t *d, const ls_vendor_finetune_t *ft,
                             ls_knots_t *knots)
{
  const ls_vendor_olympus_t *const oly = &d->u.olympus;
  const int nc = LS_MAX_KNOTS;

  float drs = 1.0f;
  float dk2 = 0.0f;
  float dk4 = 0.0f;
  float dk6 = 0.0f;
  if(oly->has_dist)
  {
    dk2 = oly->dist[0];
    dk4 = oly->dist[1];
    dk6 = oly->dist[2];
    drs = oly->dist[3];
  }
  float car0 = 0.0f;
  float car2 = 0.0f;
  float car4 = 0.0f;
  float cab0 = 0.0f;
  float cab2 = 0.0f;
  float cab4 = 0.0f;
  if(oly->has_ca)
  {
    car0 = oly->ca[0];
    car2 = oly->ca[1];
    car4 = oly->ca[2];
    cab0 = oly->ca[3];
    cab2 = oly->ca[4];
    cab4 = oly->ca[5];
  }

  for(int i = 0; i < nc; i++)
  {
    const float r = (float)i / (float)(nc - 1);
    knots->radius[i] = r;
    knots->vig_radius[i] = r;
    knots->vig[i] = 1.0f;

    float base = 1.0f;
    if(oly->has_dist)
    {
      const float rs2 = (r * drs) * (r * drs);
      const float r_cor = drs * (1.0f + rs2 * (dk2 + rs2 * (dk4 + rs2 * dk6)));
      base = ft->distortion * (r_cor - 1.0f) + 1.0f;
    }
    knots->cor_rgb[0][i] = base;
    knots->cor_rgb[1][i] = base;
    knots->cor_rgb[2][i] = base;

    if(oly->has_ca && r > 0.0f)
    {
      const float rd = base * r;
      const float rd2 = rd * rd;
      knots->cor_rgb[0][i] += ft->ca_red * (rd * (car0 + rd2 * (car2 + rd2 * car4))) / r;
      knots->cor_rgb[2][i] += ft->ca_blue * (rd * (cab0 + rd2 * (cab2 + rd2 * cab4))) / r;
    }
  }

  return nc;
}

/* ------------------------------------------------------------------------- */
/* The public entry points                                                    */
/* ------------------------------------------------------------------------- */

int ls_vendor_axes(const ls_vendor_data_t *d)
{
  if(!d) return 0;

  switch(d->type)
  {
    case LS_VENDOR_SONY:    return _sony_axes(d);
    case LS_VENDOR_FUJI:    return _fuji_axes(d);
    case LS_VENDOR_DNG:     return _dng_axes(d);
    case LS_VENDOR_OLYMPUS: return _olympus_axes(d);
    case LS_VENDOR_NONE:
    default:                return 0;
  }
}

int ls_vendor_resolve(const ls_vendor_data_t *d, const ls_vendor_finetune_t *ft,
                      int w, int h, ls_knots_t *knots, float *autoscale)
{
  if(autoscale) *autoscale = 1.0f;
  if(!knots) return LS_VENDOR_EBADDATA;
  memset(knots, 0, sizeof(*knots));
  if(!d) return 0;

  /* NULL means "as the maker measured", not "reject". The opposite convention already
   * cost a consumer a silent no-op that took a pixel-diff to find. */
  if(!ft) ft = &_ft_identity;

  /* An all-zero table with a valid type is "carries nothing", not "malformed": several
   * bodies write zero-filled tables for lenses they choose not to correct, and telling
   * their users the file could not be decoded would be noise about a non-event. Malformed
   * is reserved for tables that CLAIM an axis and then fail to convert, below. */
  const int axes = ls_vendor_axes(d);
  if(axes == 0) return 0;

  int nc = 0;
  switch(d->type)
  {
    case LS_VENDOR_SONY:    nc = _sony_populate(d, ft, knots); break;
    case LS_VENDOR_FUJI:    nc = _fuji_populate(d, ft, knots); break;
    case LS_VENDOR_DNG:     nc = _dng_populate(d, ft, knots); break;
    case LS_VENDOR_OLYMPUS: nc = _olympus_populate(d, ft, knots); break;
    default: break;
  }
  if(nc <= 0)
  {
    memset(knots, 0, sizeof(*knots));
    return LS_VENDOR_EBADDATA;
  }

  /* The autoscale: the largest correction factor over the region a border can appear in —
   * between the inscribed circle and the corner. Dividing it out of the factors and into
   * the radii re-normalises the table so the returned scale is exactly the zoom that
   * clears the borders, and the knots are neutral without it. */
  const float iwd2 = 0.5f * (float)w;
  const float iht2 = 0.5f * (float)h;
  const float diag = hypotf(iwd2, iht2);
  const float sr = fminf(iwd2, iht2);
  const float srr = (diag > 1e-6f) ? sr / diag : 0.0f;

  const int tested = 200;
  float scale = 0.0f;
  for(int i = 0; i < tested; i++)
  {
    const float x = srr + (1.0f - srr) * (float)i / (float)(tested - 1);
    for(int c = 0; c < 3; c++)
      scale = fmaxf(scale, _linear_spline(knots->radius, knots->cor_rgb[c], nc, x));
  }
  if(scale <= 1e-6f) scale = 1.0f;

  for(int i = 0; i < nc; i++)
  {
    knots->radius[i] *= scale;
    for(int c = 0; c < 3; c++) knots->cor_rgb[c][i] /= scale;
  }

  knots->n = nc;
  knots->vn = (axes & LS_ENABLE_VIGNETTING) ? nc : 0;

  if(autoscale) *autoscale = scale;
  return axes;
}

/* ------------------------------------------------------------------------- */
/* DNG OpcodeList3                                                            */
/* ------------------------------------------------------------------------- */

/* The DNG 1.3 specification, chapter 7: an opcode list is a big-endian count followed by
 * entries of [id, DNG version, flags, param_size, params]. */

#define DNG_OPCODE_ID_WARP_RECTILINEAR 1u
#define DNG_OPCODE_ID_VIGNETTE_RADIAL 3u

#define DNG_WARP_PLANES_MIN 1u
#define DNG_WARP_PLANES_MAX 3u
#define DNG_WARP_HEADER_SIZE 4u          /* the plane count */
#define DNG_WARP_PLANE_SIZE (6u * 8u)    /* six doubles per plane */
#define DNG_WARP_CENTER_SIZE (2u * 8u)   /* cx, cy */
#define DNG_VIGNETTE_COEFFS_SIZE (5u * 8u)
#define DNG_VIGNETTE_CENTER_SIZE (2u * 8u)

static double _get_be_double(const uint8_t *p)
{
  uint64_t v = 0;
  for(int i = 0; i < 8; i++) v = (v << 8) | p[i];
  double out;
  memcpy(&out, &v, sizeof(out));
  return out;
}

static uint32_t _get_be_long(const uint8_t *p)
{
  return ((uint32_t)p[0] << 24) | ((uint32_t)p[1] << 16) | ((uint32_t)p[2] << 8) | (uint32_t)p[3];
}

/* Parses a WarpRectilinear param body already known to be param_size bytes long and fully
 * inside the caller's buffer. Writes out->warp_* only on success. */
static void _parse_warp_rectilinear(const uint8_t *param, uint32_t param_size, ls_vendor_dng_t *out)
{
  if(param_size < DNG_WARP_HEADER_SIZE) return;

  const uint32_t planes = _get_be_long(&param[0]);
  if(planes < DNG_WARP_PLANES_MIN || planes > DNG_WARP_PLANES_MAX) return;
  const uint64_t min_size
      = (uint64_t)DNG_WARP_HEADER_SIZE + (uint64_t)planes * DNG_WARP_PLANE_SIZE + DNG_WARP_CENTER_SIZE;
  if(param_size < min_size) return;

  uint32_t off = DNG_WARP_HEADER_SIZE;
  for(uint32_t p = 0; p < planes; p++)
  {
    for(int c = 0; c < 6; c++)
    {
      out->warp_coeffs[p][c] = _get_be_double(&param[off]);
      off += 8;
    }
  }
  out->warp_cx = _get_be_double(&param[off]);
  out->warp_cy = _get_be_double(&param[off + 8]);
  out->warp_planes = planes;
  out->has_warp = 1;
}

/* Parses a VignetteRadial param body under the same contract. */
static void _parse_vignette_radial(const uint8_t *param, uint32_t param_size, ls_vendor_dng_t *out)
{
  if(param_size < DNG_VIGNETTE_COEFFS_SIZE + DNG_VIGNETTE_CENTER_SIZE) return;

  for(int c = 0; c < 5; c++) out->vig_coeffs[c] = _get_be_double(&param[c * 8]);
  out->vig_cx = _get_be_double(&param[DNG_VIGNETTE_COEFFS_SIZE]);
  out->vig_cy = _get_be_double(&param[DNG_VIGNETTE_COEFFS_SIZE + 8]);
  out->has_vignette = 1;
}

int ls_vendor_parse_dng_opcodelist3(const uint8_t *blob, size_t len, ls_vendor_dng_t *out)
{
  if(!out) return LS_VENDOR_EBADDATA;
  memset(out, 0, sizeof(*out));
  if(!blob || len < 4) return LS_VENDOR_EBADDATA;

  /* Bounds guard, applied before every read: the 4-byte count first; then each opcode's
   * 16-byte header against len before any header field; then the header's declared
   * param_size against len before the body is touched; then each opcode class enforces
   * its own minimum size before indexing the bytes it uses. Unknown ids are skipped —
   * no correction is invented for a class the data does not actually hold. */
  uint32_t count = _get_be_long(&blob[0]);
  uint64_t offset = 4;

  int broke = 0;
  while(count > 0)
  {
    if(offset + 16 > len) { broke = 1; break; }

    const uint32_t opcode_id = _get_be_long(&blob[offset]);
    const uint32_t param_size = _get_be_long(&blob[offset + 12]);
    const uint8_t *param = &blob[offset + 16];

    if(offset + 16 + (uint64_t)param_size > len) { broke = 1; break; }

    if(opcode_id == DNG_OPCODE_ID_WARP_RECTILINEAR)
      _parse_warp_rectilinear(param, param_size, out);
    else if(opcode_id == DNG_OPCODE_ID_VIGNETTE_RADIAL)
      _parse_vignette_radial(param, param_size, out);

    offset += 16 + (uint64_t)param_size;
    count--;
  }

  const ls_vendor_data_t probe = { .type = LS_VENDOR_DNG, .u.dng = *out };
  const int axes = ls_vendor_axes(&probe);

  /* A list that breaks after yielding a usable opcode returns what it yielded: partial
   * data is data. Only a structure broken before anything was recognised is an error. */
  if(axes == 0 && broke) return LS_VENDOR_EBADDATA;
  return axes;
}
