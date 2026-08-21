/*
    LensSerious — the evaluators, as one source text for the CPU and the GPU.

    Copyright (C) 2026 Aurélien PIERRE.
    Model evaluators, coordinate conventions and interpolation semantics ported from
    Lensfun 0.3.4 (libs/lensfun/{modifier,mod-coord,mod-subpix,mod-color}.cpp),
    Copyright (C) 2007 Andrew Zabolotny and the Lensfun contributors.

    This library is free software: you can redistribute it and/or modify it under the
    terms of the GNU Lesser General Public License as published by the Free Software
    Foundation, either version 3 of the License, or (at your option) any later version.

    It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
    without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
    PURPOSE. See the GNU Lesser General Public License for more details.
*/

/** @file lensserious_eval.h
 *
 * @brief The closed forms, written once, compiled as C99 and as OpenCL C.
 *
 * @details This file is included by src/lensserious.c and by opencl/lensserious.cl, and
 * it is the reason the two agree. The previous arrangement -- the same six models typed
 * out in both -- had already drifted: the kernel computed its square roots with
 * `native_sqrt` (implementation-defined precision, no accuracy guarantee) against the
 * library's `sqrtf`, and never grew a vignetting evaluator at all. A drift like that is
 * invisible to a parity harness that only ever tests the CPU side, which is exactly what
 * tests/parity_lensfun.c does.
 *
 * The rule for anything added here: it must compile unmodified under both toolchains.
 * That means no `float2`/`float4` (their host equivalents are not portable C), no
 * address-space qualifiers on the arguments (@ref ls_eval_t travels as a by-value kernel
 * argument, so it lands in private memory on the device and is an ordinary local on the
 * host), and arithmetic written so both sides fold it identically -- see the note on
 * argument form in @ref ls_eval_map.
 */

#ifndef LENSSERIOUS_EVAL_H
#define LENSSERIOUS_EVAL_H

/* Which toolchain is compiling this. All three spellings, because no single one is
 * dependable: OpenCL mandates __OPENCL_VERSION__ for RUNTIME compilation (clBuildProgram,
 * which is how Ansel loads this), while clang's offline OpenCL mode -- the only way to
 * syntax-check the kernel without a device, and what tests/ uses -- defines
 * __OPENCL_C_VERSION__ and not the other. Guarding on one alone means the file silently
 * takes the host branch under the other and tries to #include <math.h> into a kernel. */
#if defined(__OPENCL_VERSION__) || defined(__OPENCL_C_VERSION__) || defined(__OPENCL__)
  #define LS_EVAL_IS_OPENCL 1
#endif

#ifdef LS_EVAL_IS_OPENCL
  /* OpenCL C: sqrt() is overloaded on float and is correctly rounded to <= 3 ulp.
   * Deliberately NOT native_sqrt(): these coordinates address a resampler, and a
   * 12-bit approximation -- which is all native_sqrt() promises -- moves pixels. The
   * kernel this file replaces used native_sqrt() against the library's sqrtf(), and
   * nothing in the harness could see the difference. */
  #define LS_SQRT(x) sqrt(x)
#else
  #include <math.h>
  #define LS_SQRT(x) sqrtf(x)
#endif

/* Mirrors of the model enumerations in lensserious.h. Spelled as plain integers so this
 * file stays free of any host-only declaration; ls_eval_from_modifier() is what converts,
 * and tests/parity_lensfun.c asserts the two agree. */
#define LS_EVAL_DIST_NONE   0
#define LS_EVAL_DIST_POLY3  1
#define LS_EVAL_DIST_POLY5  2
#define LS_EVAL_DIST_PTLENS 3

#define LS_EVAL_TCA_NONE   0
#define LS_EVAL_TCA_LINEAR 1
#define LS_EVAL_TCA_POLY3  2

#define LS_EVAL_VIG_NONE 0
#define LS_EVAL_VIG_PA   1

/** Half the diagonal of a 36x24 frame, in mm: sqrt(36^2+24^2)/2. The projection focal is
 * expressed against it -- see ls_modifier_init(), which is where the conversion lives. */
#define LS_EVAL_FULL_FRAME_HALF_DIAG_MM 21.633307f

/* Mirrors of ls_lens_type_t, for the same reason as the model mirrors above. */
#define LS_EVAL_LENS_UNKNOWN                0
#define LS_EVAL_LENS_RECTILINEAR            1
#define LS_EVAL_LENS_FISHEYE                2
#define LS_EVAL_LENS_PANORAMIC              3
#define LS_EVAL_LENS_EQUIRECTANGULAR        4
#define LS_EVAL_LENS_FISHEYE_ORTHOGRAPHIC   5
#define LS_EVAL_LENS_FISHEYE_STEREOGRAPHIC  6
#define LS_EVAL_LENS_FISHEYE_EQUISOLID      7
#define LS_EVAL_LENS_FISHEYE_THOBY          8

#define LS_EVAL_ENABLE_DISTORTION (1 << 0)
#define LS_EVAL_ENABLE_TCA        (1 << 1)
#define LS_EVAL_ENABLE_VIGNETTING (1 << 2)
#define LS_EVAL_ENABLE_SCALE      (1 << 3)
#define LS_EVAL_ENABLE_GEOMETRY   (1 << 4)

/**
 * @brief One lens resolved at one shooting configuration, as a flat block of scalars.
 *
 * @details Scalars only, and no vector types: a struct of `float` and `int` has the same
 * layout under every host C ABI and under OpenCL C, so this can be memcpy'd into a kernel
 * argument. Anything richer (a `float2` centre, say) would lay out differently and corrupt
 * every field after it, silently. tests/parity_lensfun.c pins the size.
 *
 * Produced by ls_eval_from_modifier(); consumed by ls_eval_map() and
 * ls_eval_vignette_factor(), on either side.
 */
typedef struct ls_eval_t
{
  float norm_scale;              /**< centred pixel coordinates -> calibration-normalized */
  float norm_unscale;            /**< and back */
  float center_x, center_y;      /**< optical centre, in normalized coordinates */
  float aspect_ratio_correction; /**< vignetting works in the half-diagonal system */
  float scale;                   /**< linear scaling, already reciprocal; 1.0 when off */
  int   dist_model;              /**< LS_EVAL_DIST_* */
  int   tca_model;               /**< LS_EVAL_TCA_* */
  int   vig_model;               /**< LS_EVAL_VIG_* */
  int   enabled;                 /**< LS_EVAL_ENABLE_* actually resolved */
  int   geom_from;               /**< the lens's own projection, LS_EVAL_LENS_* */
  int   geom_to;                 /**< the projection asked for */
  float geom_focal;              /**< focal / LS_EVAL_GEOM_HALF_HEIGHT_MM */
  float _pad;                    /**< keeps the struct 4-aligned and its size stated */
  float dist_terms[3];
  float tca_terms[6];
  float vig_terms[3];
} ls_eval_t;

/** @brief Distortion, in place, in normalized coordinates. */
static inline void ls_eval_dist(const ls_eval_t *p, float *x, float *y)
{
  const float xu = *x, yu = *y;
  const float ru2 = xu * xu + yu * yu;
  float m = 1.f;

  if(p->dist_model == LS_EVAL_DIST_POLY3)
  {   /* mod-coord.cpp: Rd = Ru · (1 − k1 + k1·Ru²) */
    const float k1 = p->dist_terms[0];
    m = (1.f - k1) + k1 * ru2;
  }
  else if(p->dist_model == LS_EVAL_DIST_POLY5)
  {   /* Rd = Ru · (1 + k1·Ru² + k2·Ru⁴) */
    m = 1.f + p->dist_terms[0] * ru2 + p->dist_terms[1] * ru2 * ru2;
  }
  else if(p->dist_model == LS_EVAL_DIST_PTLENS)
  {   /* Rd = Ru · (a·Ru³ + b·Ru² + c·Ru + d), d = 1−a−b−c */
    const float a = p->dist_terms[0], b = p->dist_terms[1], c = p->dist_terms[2];
    const float r = LS_SQRT(ru2);
    m = a * ru2 * r + b * ru2 + c * r + (1.f - a - b - c);
  }

  *x = xu * m;
  *y = yu * m;
}

/** @brief Transverse chromatic aberration: the red and blue coordinates diverge from green. */
static inline void ls_eval_tca(const ls_eval_t *p, float *xr, float *yr, float *xb, float *yb)
{
  if(p->tca_model == LS_EVAL_TCA_LINEAR)
  {   /* mod-subpix.cpp: per-channel radial scale */
    const float kr = p->tca_terms[0], kb = p->tca_terms[1];
    *xr *= kr; *yr *= kr;
    *xb *= kb; *yb *= kb;
  }
  else if(p->tca_model == LS_EVAL_TCA_POLY3)
  {   /* Rd = Ru · (b·Ru² + c·Ru + v), terms packed vr vb cr cb br bb */
    const float vr = p->tca_terms[0], vb = p->tca_terms[1];
    const float cr = p->tca_terms[2], cb = p->tca_terms[3];
    const float br = p->tca_terms[4], bb = p->tca_terms[5];

    float x = *xr, y = *yr;
    float ru2 = x * x + y * y;
    float m = br * ru2 + vr + ((cr != 0.f) ? cr * LS_SQRT(ru2) : 0.f);
    *xr = x * m; *yr = y * m;

    x = *xb; y = *yb;
    ru2 = x * x + y * y;
    m = bb * ru2 + vb + ((cb != 0.f) ? cb * LS_SQRT(ru2) : 0.f);
    *xb = x * m; *yb = y * m;
  }
}

/**
 * @brief The map for ONE output pixel: six floats, source coordinates for R, G, B.
 *
 * @param xu absolute output column, in pixels. @param yu absolute output row.
 * @param out six floats: xr yr xg yg xb yb, the layout of
 * lfModifier::ApplySubpixelGeometryDistortion()'s buffer.
 *
 * @details @p xu and @p yu are ABSOLUTE, already summed by the caller -- a row walker
 * passes `xu + col`, a work-item passes `xu + get_global_id(0)`. Both then evaluate the
 * identical expression on the identical float, which is what makes the CPU row loop and
 * the kernel agree bit for bit rather than merely closely. Do not "optimise" this into
 * taking an origin plus an index: the two would then round differently.
 */
/**
 * @brief Field angle for a radius, under one projection. Negative if @p model has none.
 *
 * @details The destination half of a projection change: the output image is in @p model's
 * geometry, and this recovers the angle the ray came in at. Radii outside a projection's
 * domain (a fisheye only maps so far) return -1 so the caller can leave the pixel black
 * rather than fold the image back on itself.
 */
static inline float ls_eval_geom_angle(const int model, const float f, const float r)
{
  if(f <= 0.f) return -1.f;
  switch(model)
  {
    case LS_EVAL_LENS_RECTILINEAR:           return atan(r / f);
    case LS_EVAL_LENS_FISHEYE:               return r / f;
    case LS_EVAL_LENS_FISHEYE_ORTHOGRAPHIC:  return (r <= f) ? asin(r / f) : -1.f;
    case LS_EVAL_LENS_FISHEYE_STEREOGRAPHIC: return 2.f * atan(r / (2.f * f));
    case LS_EVAL_LENS_FISHEYE_EQUISOLID:     return (r <= 2.f * f) ? 2.f * asin(r / (2.f * f)) : -1.f;
    case LS_EVAL_LENS_FISHEYE_THOBY:
      return (r <= 1.47f * f) ? asin(r / (1.47f * f)) / 0.713f : -1.f;
    default: return -1.f;   /* panoramic and equirectangular are not radial; see the header */
  }
}

/** @brief Radius for a field angle, under one projection. The inverse of the above. */
static inline float ls_eval_geom_radius(const int model, const float f, const float theta)
{
  if(f <= 0.f || theta < 0.f) return -1.f;
  switch(model)
  {
    case LS_EVAL_LENS_RECTILINEAR:
      /* Beyond a right angle a rectilinear lens images nothing: tan() would silently wrap
       * a point behind the camera round to the front. */
      return (theta < 1.5707963f) ? f * tan(theta) : -1.f;
    case LS_EVAL_LENS_FISHEYE:               return f * theta;
    case LS_EVAL_LENS_FISHEYE_ORTHOGRAPHIC:  return f * sin(theta);
    case LS_EVAL_LENS_FISHEYE_STEREOGRAPHIC: return 2.f * f * tan(theta * 0.5f);
    case LS_EVAL_LENS_FISHEYE_EQUISOLID:     return 2.f * f * sin(theta * 0.5f);
    case LS_EVAL_LENS_FISHEYE_THOBY:         return 1.47f * f * sin(0.713f * theta);
    default: return -1.f;
  }
}

/**
 * @brief Reproject one point from the target geometry into the lens's own.
 *
 * @details Radially symmetric, so it scales x and y by one factor rather than mapping them
 * separately -- which is exactly why panoramic and equirectangular are excluded: those two
 * treat the axes differently and cannot be expressed this way.
 *
 * @return 0 if the point has no source (outside the projection's domain), 1 otherwise.
 */
static inline int ls_eval_geometry(const ls_eval_t *p, float *x, float *y)
{
  const float r = LS_SQRT((*x) * (*x) + (*y) * (*y));
  if(r <= 0.f) return 1;                       /* the centre maps to itself */

  const float theta = ls_eval_geom_angle(p->geom_to, p->geom_focal, r);
  if(theta < 0.f) return 0;
  const float r_src = ls_eval_geom_radius(p->geom_from, p->geom_focal, theta);
  if(r_src < 0.f) return 0;

  const float k = r_src / r;
  *x *= k;
  *y *= k;
  return 1;
}

static inline void ls_eval_map(const ls_eval_t *p, float xu, float yu, float *out)
{
  float x = xu * p->norm_scale - p->center_x;
  float y = yu * p->norm_scale - p->center_y;

  /* lensfun's callback-priority order for the correction direction: scaling (100) ->
   * projection (500) -> distortion (750), then the TCA subpixel stage per channel. */
  if(p->enabled & LS_EVAL_ENABLE_SCALE)
  {
    x *= p->scale;
    y *= p->scale;
  }
  if(p->enabled & LS_EVAL_ENABLE_GEOMETRY)
  {
    if(!ls_eval_geometry(p, &x, &y))
    {
      /* No source pixel. NaN is what lensfun writes here too, and every consumer in this
       * project already checks for it (do_nan_checks in the kernels, isfinite() on the
       * CPU) before sampling. */
      for(int k = 0; k < 6; k++) out[k] = (float)(0.0f / 0.0f);
      return;
    }
  }
  if(p->enabled & LS_EVAL_ENABLE_DISTORTION) ls_eval_dist(p, &x, &y);

  float xr = x, yr = y, xb = x, yb = y;
  if(p->enabled & LS_EVAL_ENABLE_TCA) ls_eval_tca(p, &xr, &yr, &xb, &yb);

  out[0] = (xr + p->center_x) * p->norm_unscale;
  out[1] = (yr + p->center_y) * p->norm_unscale;
  out[2] = (x  + p->center_x) * p->norm_unscale;
  out[3] = (y  + p->center_y) * p->norm_unscale;
  out[4] = (xb + p->center_x) * p->norm_unscale;
  out[5] = (yb + p->center_y) * p->norm_unscale;
}

/**
 * @brief The vignetting multiplier for ONE output pixel. Multiply the pixel by it.
 *
 * @details Returns 1.0 when no vignetting is resolved, so a caller may apply it
 * unconditionally. Absolute pixel coordinates, as in ls_eval_map().
 *
 * Two conventions that a reimplementation gets wrong and a harness catches late:
 * "Damn! Hugin uses two different 'normalized' coordinate systems: for distortions it uses
 * 1.0 = min(half width, half height) and for vignetting it uses 1.0 = half diagonal
 * length" (mod-color.cpp) -- hence the division by the aspect-ratio correction; and
 * CORRECTING vignetting DIVIDES by the polynomial (ModifyColor_DeVignetting_PA applies
 * 1/c), where multiplying by c is what applies it.
 */
static inline float ls_eval_vignette_factor(const ls_eval_t *p, float xu, float yu)
{
  if(!(p->enabled & LS_EVAL_ENABLE_VIGNETTING)) return 1.f;

  const float arc = (p->aspect_ratio_correction > 0.f) ? p->aspect_ratio_correction : 1.f;
  const float x = (xu * p->norm_scale - p->center_x) / arc;
  const float y = (yu * p->norm_scale - p->center_y) / arc;

  const float r2 = x * x + y * y;
  const float r4 = r2 * r2;
  const float c = 1.f + p->vig_terms[0] * r2 + p->vig_terms[1] * r4
                      + p->vig_terms[2] * r4 * r2;
  return (c != 0.f) ? 1.f / c : 1.f;
}

#endif /* LENSSERIOUS_EVAL_H */
