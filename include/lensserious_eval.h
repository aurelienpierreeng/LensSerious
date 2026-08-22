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

/* Every transcendental goes through one of these, and the reason is single precision.
 *
 * In OpenCL C these names ARE the float overloads. In C they are the DOUBLE functions:
 * `atan(some_float)` promotes its argument to double, evaluates in double, and truncates
 * the result back on assignment. That is slower than the float form, and it stops the
 * loops below from vectorising, since a `double` temporary halves the lane count. Writing
 * the `f`-suffixed name in the host branch is what keeps this file honestly single
 * precision on both sides.
 *
 * Deliberately NOT the native_* family on the OpenCL side: these results address a
 * resampler, and native_sqrt() promises only ~12 bits. The kernel this file replaced used
 * native_sqrt() against the library's sqrtf(), and nothing in the harness could see it. */
#ifdef LS_EVAL_IS_OPENCL
  #define LS_SQRT(x)  sqrt(x)
  #define LS_ATAN(x)  atan(x)
  #define LS_ASIN(x)  asin(x)
  #define LS_SIN(x)   sin(x)
  #define LS_TAN(x)   tan(x)
  #define LS_FABS(x)  fabs(x)
#else
  #include <math.h>
  #define LS_SQRT(x)  sqrtf(x)
  #define LS_ATAN(x)  atanf(x)
  #define LS_ASIN(x)  asinf(x)
  #define LS_SIN(x)   sinf(x)
  #define LS_TAN(x)   tanf(x)
  #define LS_FABS(x)  fabsf(x)
#endif

/* Mirrors of the model enumerations in lensserious.h. Spelled as plain integers so this
 * file stays free of any host-only declaration; ls_eval_from_modifier() is what converts,
 * and tests/parity_lensfun.c asserts the two agree. */
/* A radial correction given as a TABLE rather than a polynomial: vendors embed their own
 * profiles in the raw's metadata that way, as a handful of knots the manufacturer measured.
 * Everything else about the model is unchanged -- it is still "scale the centred coordinate
 * by a factor of its radius", which is what every model here does -- so this rides through
 * the same evaluator, the same kernels and the same by-value transport. */
#define LS_EVAL_DIST_KNOTS  4
#define LS_EVAL_VIG_KNOTS   2

/* Sixteen matches what the vendor formats actually carry: Fuji ships nine distortion knots,
 * Sony and Olympus fewer. It brings ls_eval_t to 632 bytes, which the by-value kernel
 * argument still absorbs -- the smallest CL_DEVICE_MAX_PARAMETER_SIZE that OpenCL 1.2
 * guarantees is 1024, for a kernel's whole argument list, and the measured floor on this
 * project's devices is the same. src/lensserious.c asserts both the size and that headroom,
 * so raising this fails the build rather than the device. */
#define LS_MAX_KNOTS 16

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
  /* Vignetting works in the half-diagonal system, i.e. the geometry coordinates divided by
   * the aspect-ratio correction. Folded into a scale and an offset here rather than stored
   * as the correction itself, so the per-pixel path needs one multiply and one subtract
   * instead of a divide -- per axis, per pixel. That was three divides per pixel against
   * lensfun's one. */
  float vig_scale;               /**< norm_scale / aspect_ratio_correction */
  float vig_center_x, vig_center_y; /**< centre / aspect_ratio_correction */
  float scale;                   /**< linear scaling, already reciprocal; 1.0 when off */
  int   dist_model;              /**< LS_EVAL_DIST_* */
  int   tca_model;               /**< LS_EVAL_TCA_* */
  int   vig_model;               /**< LS_EVAL_VIG_* */
  int   enabled;                 /**< LS_EVAL_ENABLE_* actually resolved */
  int   geom_from;               /**< the lens's own projection, LS_EVAL_LENS_* */
  int   geom_to;                 /**< the projection asked for */
  float geom_focal;              /**< focal / LS_EVAL_GEOM_HALF_HEIGHT_MM */
  /** Non-zero when this block describes the REVERSE direction (lensfun's `reverse` flag).
   *
   * It is not "the same chain, backwards". Upstream registers different callbacks at
   * different priorities, so the composition order itself changes: forward runs
   * scale(100) -> geometry(500) -> distortion(750), reverse runs undistortion(250) ->
   * geometry(500) -> scale(900). The projection endpoints are swapped and the scale factor
   * is stored un-reciprocated by the resolver, so only the ORDER and the choice of
   * dist/undist live here. */
  int   reverse;
  float dist_terms[3];
  float tca_terms[6];
  float vig_terms[3];

  /* Knot tables, used when dist_model / vig_model say so. Per CHANNEL for the geometry,
   * because that is what a vendor profile is: one radial curve per channel, which is
   * distortion and TCA expressed together rather than as two stages.
   *
   * The radii are per-channel too. Going forwards that is redundant -- a vendor gives ONE
   * set of knots shared by all three curves -- and it is what makes the REVERSE direction
   * exact rather than merely close. Inverting r_src = r * cor(r) has no closed form in
   * general, but a piecewise-linear curve inverts exactly at its own knots: (r*cor(r),
   * 1/cor(r)) is the same curve read the other way. Each channel's inverse lands on its own
   * radii, and forcing all three back onto a shared axis would mean resampling two of them.
   * The resolver builds whichever direction was asked for, so the evaluator below never
   * inverts anything -- no Newton, no per-pixel search, unlike the polynomial models. */
  int   knot_n;
  float knot_r[3][LS_MAX_KNOTS];
  float knot_c[3][LS_MAX_KNOTS];
  int   knot_vn;
  float knot_vr[LS_MAX_KNOTS];
  float knot_v[LS_MAX_KNOTS];
} ls_eval_t;

/**
 * @brief Piecewise-linear lookup over a knot table.
 *
 * @param xs the knot positions, ascending, @p n of them.
 * @param ys the value at each of them.
 * @param n how many knots.
 * @param x where to sample. Clamped to the ends rather than extrapolated.
 *
 * @details The vendors' own convention, and deliberately not something smoother: their
 * profiles are DEFINED as the piecewise-linear interpolant of these points, so a spline
 * through them would be a different correction, not a better one.
 */
static inline float ls_eval_knot_lookup(const float *xs, const float *ys, const int n,
                                        const float x)
{
  if(n <= 0) return 1.f;
  if(x <= xs[0]) return ys[0];

  for(int i = 1; i < n; i++)
  {
    if(x <= xs[i])
    {
      const float d = xs[i] - xs[i - 1];
      if(d <= 0.f) return ys[i - 1];
      return ys[i - 1] + (x - xs[i - 1]) * (ys[i] - ys[i - 1]) / d;
    }
  }
  return ys[n - 1];
}

/**
 * @brief The per-channel radial scale a knot table gives at this point.
 *
 * @param p the lens resolved at one shooting configuration.
 * @param c the channel, 0 red / 1 green / 2 blue.
 * @param x, y a point in normalized coordinates.
 *
 * @details No conversion between radius conventions, because there is none to make:
 * ls_modifier_init_knots() defines this modifier's normalized system to BE the one the
 * vendors index their tables by -- the distance from the image centre over half the image
 * diagonal, so 1.0 at the corner. A knot table never shares a modifier with lens-database
 * calibration data, so nothing has to reconcile the two.
 */
static inline float ls_eval_knot_factor(const ls_eval_t *p, const int c, const float x,
                                        const float y)
{
  const float r = LS_SQRT(x * x + y * y);
  return ls_eval_knot_lookup(p->knot_r[c], p->knot_c[c], p->knot_n, r);
}

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
    case LS_EVAL_LENS_RECTILINEAR:           return LS_ATAN(r / f);
    case LS_EVAL_LENS_FISHEYE:               return r / f;
    case LS_EVAL_LENS_FISHEYE_ORTHOGRAPHIC:  return (r <= f) ? LS_ASIN(r / f) : -1.f;
    case LS_EVAL_LENS_FISHEYE_STEREOGRAPHIC: return 2.f * LS_ATAN(r / (2.f * f));
    case LS_EVAL_LENS_FISHEYE_EQUISOLID:     return (r <= 2.f * f) ? 2.f * LS_ASIN(r / (2.f * f)) : -1.f;
    case LS_EVAL_LENS_FISHEYE_THOBY:
      return (r <= 1.47f * f) ? LS_ASIN(r / (1.47f * f)) / 0.713f : -1.f;
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
      /* Beyond a right angle a rectilinear lens images nothing: LS_TAN() would silently wrap
       * a point behind the camera round to the front. */
      return (theta < 1.5707963f) ? f * LS_TAN(theta) : -1.f;
    case LS_EVAL_LENS_FISHEYE:               return f * theta;
    case LS_EVAL_LENS_FISHEYE_ORTHOGRAPHIC:  return f * LS_SIN(theta);
    case LS_EVAL_LENS_FISHEYE_STEREOGRAPHIC: return 2.f * f * LS_TAN(theta * 0.5f);
    case LS_EVAL_LENS_FISHEYE_EQUISOLID:     return 2.f * f * LS_SIN(theta * 0.5f);
    case LS_EVAL_LENS_FISHEYE_THOBY:         return 1.47f * f * LS_SIN(0.713f * theta);
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

/* Newton, ported from mod-coord.cpp with ONE deliberate difference, which cost a bench
 * round to find and is the reason this comment is long.
 *
 * Upstream iterates in double and tests an ABSOLUTE residual, |f(ru)| < 1e-5. For poly3 it
 * also divides the whole equation by k1 to make it monic: ru^3 + ru*(1-k1)/k1 - rd/k1. That
 * is free in double and fatal in float. Real lenses have small k1 -- the Beroflex 500mm has
 * k1 = 0.00108 -- so the scaled equation carries terms of magnitude 1/k1 = 926, whose float
 * rounding noise alone is ~5.5e-5. The residual can then NEVER fall below 1e-5, every pixel
 * exhausts the six-step budget, and the coordinate is returned uncorrected: measured as
 * 21177 samples out of tolerance, up to 6 px, on 0.6% of the reverse database.
 *
 * So the equations here are the UNSCALED ones -- same roots, magnitudes near 1 -- and
 * convergence is judged on the relative STEP SIZE rather than an absolute residual. That is
 * scale-free, which is what float needs, and it tests the thing actually wanted: that the
 * iteration has stopped moving. The budget stays at upstream's six steps and a
 * non-converging pixel is still left untouched, so genuinely divergent regions (ultrawide
 * fisheye corners, where upstream gives up too) behave the same on both sides.
 *
 * Iterating in float rather than double is not a shortcut either: one source text has to
 * compile as OpenCL C, where double is an optional extension, and CPU/GPU bit-exactness is
 * a property this project asserts. The cost of that choice is measured over the whole
 * database in both directions by tests/parity_lensfun.c. */
#define LS_NEWTON_STEPS 6
#define LS_NEWTON_RTOL  1e-6f

/**
 * @brief Solve Rd = f(Ru) for Ru, given a radius Rd. The inverse of ls_eval_dist().
 *
 * @param p the lens resolved at one shooting configuration, reverse direction.
 * @param rd the distorted radius, in normalized units.
 * @return the factor Ru/Rd to scale the coordinate by, or 1.0 when the iteration does not
 * converge or lands on a negative radius -- upstream leaves such a coordinate untouched,
 * and returning 1 here is that same decision expressed as a multiplier.
 */
static inline float ls_eval_undist_factor(const ls_eval_t *p, const float rd)
{
  if(rd == 0.f) return 1.f;

  float ru = rd;
  int converged = 0;

  if(p->dist_model == LS_EVAL_DIST_POLY3)
  {
    /* Rd = k1*Ru^3 + (1-k1)*Ru, solved as written -- see the note above on why this is not
     * divided through by k1 the way upstream does it. */
    const float k1 = p->dist_terms[0], one_minus_k1 = 1.f - k1;
    for(int step = 0; step < LS_NEWTON_STEPS && !converged; step++)
    {
      const float f = k1 * ru * ru * ru + one_minus_k1 * ru - rd;
      const float fp = 3.f * k1 * ru * ru + one_minus_k1;
      if(fp == 0.f) return 1.f;
      const float d = f / fp;
      ru -= d;
      if(LS_FABS(d) <= LS_NEWTON_RTOL * LS_FABS(ru)) converged = 1;
    }
  }
  else if(p->dist_model == LS_EVAL_DIST_POLY5)
  {
    const float k1 = p->dist_terms[0], k2 = p->dist_terms[1];
    for(int step = 0; step < LS_NEWTON_STEPS && !converged; step++)
    {
      const float ru2 = ru * ru;
      const float f = ru * (1.f + k1 * ru2 + k2 * ru2 * ru2) - rd;
      const float fp = 1.f + 3.f * k1 * ru2 + 5.f * k2 * ru2 * ru2;
      if(fp == 0.f) return 1.f;
      const float d = f / fp;
      ru -= d;
      if(LS_FABS(d) <= LS_NEWTON_RTOL * LS_FABS(ru)) converged = 1;
    }
  }
  else if(p->dist_model == LS_EVAL_DIST_PTLENS)
  {
    const float a = p->dist_terms[0], b = p->dist_terms[1], c = p->dist_terms[2];
    const float d0 = 1.f - a - b - c;
    for(int step = 0; step < LS_NEWTON_STEPS && !converged; step++)
    {
      const float f = ru * (a * ru * ru * ru + b * ru * ru + c * ru + d0) - rd;
      const float fp = 4.f * a * ru * ru * ru + 3.f * b * ru * ru + 2.f * c * ru + d0;
      if(fp == 0.f) return 1.f;
      const float d = f / fp;
      ru -= d;
      if(LS_FABS(d) <= LS_NEWTON_RTOL * LS_FABS(ru)) converged = 1;
    }
  }
  else
    return 1.f;

  if(!converged || ru < 0.f) return 1.f;  /* a negative radius is not a position */
  return ru / rd;
}

/** @brief Undistort in place. The reverse of ls_eval_dist(). */
static inline void ls_eval_undist(const ls_eval_t *p, float *x, float *y)
{
  const float rd = LS_SQRT(*x * *x + *y * *y);
  const float m = ls_eval_undist_factor(p, rd);
  *x *= m;
  *y *= m;
}

/**
 * @brief The reverse of ls_eval_tca(): recover each channel's undistorted radius.
 *
 * @details Linear is a reciprocal the resolver already took. Poly3 is the same Newton as
 * above on Rd = b·Ru³ + c·Ru² + v·Ru, run once per channel, and -- exactly as upstream --
 * a channel that fails to converge keeps the coordinate it came in with while the other
 * channel may still be corrected.
 */
static inline void ls_eval_untca(const ls_eval_t *p, float *xr, float *yr, float *xb, float *yb)
{
  if(p->tca_model == LS_EVAL_TCA_LINEAR)
  {
    const float kr = p->tca_terms[0], kb = p->tca_terms[1];  /* already reciprocals */
    *xr *= kr; *yr *= kr;
    *xb *= kb; *yb *= kb;
  }
  else if(p->tca_model == LS_EVAL_TCA_POLY3)
  {
    const float vr = p->tca_terms[0], vb = p->tca_terms[1];
    const float cr = p->tca_terms[2], cb = p->tca_terms[3];
    const float br = p->tca_terms[4], bb = p->tca_terms[5];

    for(int ch = 0; ch < 2; ch++)
    {
      const float v = ch ? vb : vr, c = ch ? cb : cr, b = ch ? bb : br;
      float *px = ch ? xb : xr, *py = ch ? yb : yr;
      const float rd = LS_SQRT(*px * *px + *py * *py);
      if(rd == 0.f) continue;

      float ru = rd;
      int ok = 0;
      for(int step = 0; step < LS_NEWTON_STEPS && !ok; step++)
      {
        const float ru2 = ru * ru;
        const float f = b * ru2 * ru + c * ru2 + v * ru - rd;
        const float fp = 3.f * b * ru2 + 2.f * c * ru + v;
        if(fp == 0.f) break;
        const float dd = f / fp;
        ru -= dd;
        if(LS_FABS(dd) <= LS_NEWTON_RTOL * LS_FABS(ru)) ok = 1;
      }
      /* Upstream requires ru STRICTLY positive here, where the coordinate path accepts
       * zero. Kept as it is rather than unified: they are different functions. */
      if(ok && ru > 0.f)
      {
        const float m = ru / rd;
        *px *= m; *py *= m;
      }
    }
  }
}

/**
 * @brief The coordinate chain: scale, projection and distortion, in direction order.
 *
 * @param p the lens resolved at one shooting configuration.
 * @param c which channel's geometry to follow, 0 red / 1 green / 2 blue. Only a knot table
 * distinguishes them -- a vendor profile IS one radial curve per channel, distortion and TCA
 * measured together rather than modelled as two stages. Every polynomial model puts all
 * three on the green curve and expresses the difference as a separate TCA stage afterwards,
 * so for those this argument is ignored. Pass 1 wherever a single geometry is meant.
 * @param x, y a point in normalized coordinates, transformed in place.
 * @return 0 when the point has no source pixel at all, in which case @p x and @p y are
 * untouched and the caller decides what to write.
 *
 * @details Everything lensfun registers as a COORDINATE callback and nothing else -- TCA is
 * a subpixel callback and runs after, vignetting is a colour callback. It exists as its own
 * function because autoscaling has to evaluate exactly this, and a second copy written for
 * that purpose would be a second copy to drift.
 *
 * The order is not symmetric between the directions:
 *   forward   scale (100)        -> projection (500) -> distortion (750)
 *   reverse   undistortion (250) -> projection (500) -> scale (900)
 * The projection sits in the middle either way; scale moves from first to last, and the
 * resolver has already swapped the projection endpoints and un-reciprocated the scale.
 *
 * A knot table occupies the distortion slot, and only that slot: the vendor measured the
 * lens as shipped, so there is no projection change and no scale baked into its numbers, and
 * a consumer that asks for either still gets it composed around the table in the order above.
 * It sits in the FORWARD slot in both directions because the resolver builds the table for
 * the direction asked for -- inverting a piecewise-linear curve is exact at its own knots,
 * which is why nothing here iterates the way ls_eval_undist() has to.
 */
static inline int ls_eval_coord_chain(const ls_eval_t *p, const int c, float *x, float *y)
{
  if(!p->reverse && (p->enabled & LS_EVAL_ENABLE_SCALE))
  {
    *x *= p->scale;
    *y *= p->scale;
  }
  if(p->reverse && p->dist_model != LS_EVAL_DIST_KNOTS
     && (p->enabled & LS_EVAL_ENABLE_DISTORTION))
    ls_eval_undist(p, x, y);

  if((p->enabled & LS_EVAL_ENABLE_GEOMETRY) && !ls_eval_geometry(p, x, y)) return 0;

  if(p->dist_model == LS_EVAL_DIST_KNOTS)
  {
    if(p->enabled & LS_EVAL_ENABLE_DISTORTION)
    {
      const float f = ls_eval_knot_factor(p, c, *x, *y);
      *x *= f;
      *y *= f;
    }
  }
  else if(!p->reverse && (p->enabled & LS_EVAL_ENABLE_DISTORTION))
    ls_eval_dist(p, x, y);
  if(p->reverse && (p->enabled & LS_EVAL_ENABLE_SCALE))
  {
    *x *= p->scale;
    *y *= p->scale;
  }
  return 1;
}

/**
 * @brief The map for ONE output pixel: six floats, source coordinates for R, G, B.
 *
 * @param p the lens resolved at one shooting configuration.
 * @param xu absolute output column, in pixels.
 * @param yu absolute output row, in pixels.
 * @param out six floats: xr yr xg yg xb yb, the layout of
 * lfModifier::ApplySubpixelGeometryDistortion()'s buffer.
 *
 * @details @p xu and @p yu are ABSOLUTE, already summed by the caller -- a row walker
 * passes `xu + col`, a work-item passes `xu + get_global_id(0)`. Both then evaluate the
 * identical expression on the identical float, which is what makes the CPU row loop and
 * the kernel agree bit for bit rather than merely closely. Do not "optimise" this into
 * taking an origin plus an index: the two would then round differently.
 *
 * @see @ref fused-evaluation for why a consumer should usually call this in the loop that
 * needs the coordinates, rather than filling a buffer with it.
 */
static inline void ls_eval_map(const ls_eval_t *p, float xu, float yu, float *out)
{
  float x = xu * p->norm_scale - p->center_x;
  float y = yu * p->norm_scale - p->center_y;

  /* lensfun's callback-priority order, which is NOT symmetric between the two directions:
   *   forward   scale (100)        -> projection (500) -> distortion (750)
   *   reverse   undistortion (250) -> projection (500) -> scale (900)
   * Projection sits in the middle either way; scale moves from first to last, and the
   * resolver has already swapped the projection endpoints and un-reciprocated the scale.
   * The TCA subpixel stage runs after the coordinate chain in both directions. */
  float xr = x, yr = y, xb = x, yb = y;

  if(p->dist_model == LS_EVAL_DIST_KNOTS)
  {
    /* Three chains, one per channel, because a vendor profile has no single geometry to run
     * a TCA correction on top of -- the red, green and blue curves ARE the correction, and
     * the difference between them is the lateral chromatic aberration. Running the chain
     * three times rather than once-plus-a-delta is what keeps that exact.
     *
     * It is also cheap where it is used: a vendor table is measured on the lens as shipped,
     * so the projection stage is an identity that ls_eval_coord_chain() skips outright, and
     * what repeats is a multiply and a table lookup.
     *
     * Only the green chain decides whether there is a source pixel at all. The three differ
     * by the TCA amount, tenths of a percent, so a per-channel verdict here could hand the
     * caller two valid channels and one NaN for the same pixel. */
    ls_eval_coord_chain(p, 0, &xr, &yr);
    ls_eval_coord_chain(p, 2, &xb, &yb);
    if(!ls_eval_coord_chain(p, 1, &x, &y))
    {
      for(int k = 0; k < 6; k++) out[k] = (float)(0.0f / 0.0f);
      return;
    }
  }
  else
  {
    if(!ls_eval_coord_chain(p, 1, &x, &y))
    {
      /* No source pixel. NaN is what lensfun writes here too, and every consumer in this
       * project already checks for it (do_nan_checks in the kernels, isfinite() on the CPU)
       * before sampling. */
      for(int k = 0; k < 6; k++) out[k] = (float)(0.0f / 0.0f);
      return;
    }
    xr = x; yr = y; xb = x; yb = y;
    if(p->enabled & LS_EVAL_ENABLE_TCA)
    {
      if(p->reverse) ls_eval_untca(p, &xr, &yr, &xb, &yb);
      else           ls_eval_tca(p, &xr, &yr, &xb, &yb);
    }
  }

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
/**
 * @brief The vignetting multiplier for a squared radius already in the half-diagonal system.
 *
 * @details Split out so a row walker can hoist everything that does not vary along a row --
 * the y term and its square -- while the polynomial itself, the pole handling and the clamp
 * stay in ONE place for both the CPU and the kernel. Callers that have a pixel coordinate
 * rather than a radius want ls_eval_vignette_factor().
 */
static inline float ls_eval_vignette_from_r2(const ls_eval_t *p, const float r2)
{
  if(p->vig_model == LS_EVAL_VIG_KNOTS)
  {
    /* The table already states the multiplier that CORRECTS the falloff, so the direction
     * is handled the same way the polynomial's is -- reversing puts it back. */
    const float r = LS_SQRT(r2);
    const float v = ls_eval_knot_lookup(p->knot_vr, p->knot_v, p->knot_vn, r);
    /* The table states the falloff -- what the lens DID -- so correcting divides by it and
     * applying multiplies, which is the same way round as the polynomial below and the same
     * way round as lensfun's two callbacks. */
    const float m = p->reverse ? v : ((v != 0.f) ? (1.f / v) : 0.f);
    return (m > 0.f) ? m : 0.f;
  }

  const float r4 = r2 * r2;
  const float c = 1.f + p->vig_terms[0] * r2 + p->vig_terms[1] * r4
                      + p->vig_terms[2] * r4 * r2;

  /* Clamped at zero, because the pa polynomial is not constrained to stay positive and for
   * some lenses it crosses zero INSIDE the frame -- the Canon EF 8-15mm Fisheye at 8mm has
   * k = (-0.625, 5.648, -19.330), whose root sits near r = 0.65. Past that root 1/c is
   * negative, which is not a brightness.
   *
   * Upstream clamps the same thing one step later: apply_multiplier() writes
   * clampd(pixel * c, 0, type_max) (mod-color.cpp), so for the non-negative pixels this
   * ever sees, clamping the factor here is equivalent. Without it this returned -0.203
   * where lensfun returns 0, on every fisheye whose model has a root in frame. */
  /* One divide and one select, and nothing else -- which is exactly what upstream does:
   * ModifyColor_DeVignetting_PA divides by c with no special case for zero, and
   * apply_multiplier() clamps the product at zero afterwards. An earlier version guarded
   * c == 0 here; that guard was not upstream behaviour AND it was a second conditional,
   * which is what stopped a caller's row loop from vectorising. Dividing by zero yields an
   * infinity that the clamp below leaves alone for +0 and flattens for -0, matching
   * upstream on a case that no real calibration reaches anyway.
   *
   * Keeping this to a single select matters more than it looks: it is the difference
   * between one divide per pixel and one per four. */
  /* The DIRECTION decides whether the lens's falloff is removed or re-applied, and
   * upstream implements the two as separate callbacks: ModifyColor_DeVignetting_PA
   * multiplies by 1/c to correct, ModifyColor_Vignetting_PA multiplies by c to put the
   * falloff back (mod-color.cpp, priorities 750 and 250). Ignoring the flag here meant the
   * reverse direction BRIGHTENED the corners it was supposed to darken -- the correction
   * applied twice instead of undone.
   *
   * Still one select on top of the divide: the reciprocal is computed either way, and the
   * branch is on a field that is constant for the whole frame. */
  const float m = p->reverse ? c : (1.f / c);
  return (m > 0.f) ? m : 0.f;
}

/** @brief The vignetting multiplier for ONE output pixel. Multiply the pixel by it.
 *
 * @details Returns 1.0 when no vignetting is resolved, so a caller may apply it
 * unconditionally. Absolute pixel coordinates, as in ls_eval_map(). */
static inline float ls_eval_vignette_factor(const ls_eval_t *p, float xu, float yu)
{
  if(!(p->enabled & LS_EVAL_ENABLE_VIGNETTING)) return 1.f;

  const float x = xu * p->vig_scale - p->vig_center_x;
  const float y = yu * p->vig_scale - p->vig_center_y;
  return ls_eval_vignette_from_r2(p, x * x + y * y);
}

#endif /* LENSSERIOUS_EVAL_H */
