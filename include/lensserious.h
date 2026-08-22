/*
    LensSerious — lens-correction mathematics as data, not as a library of callbacks.

    Copyright (C) 2026 Aurélien PIERRE.
    Model evaluators, coordinate conventions and interpolation semantics ported from
    Lensfun 0.3.4 (libs/lensfun/{modifier,mod-coord,mod-subpix,mod-color,lens,auxfun}.cpp),
    Copyright (C) 2007 Andrew Zabolotny and the Lensfun contributors.

    This library is free software: you can redistribute it and/or modify it under the
    terms of the GNU Lesser General Public License as published by the Free Software
    Foundation, either version 3 of the License, or (at your option) any later version.

    It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
    without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
    PURPOSE. See the GNU Lesser General Public License for more details.
*/

/** @file lensserious.h
 *
 * @brief What this is, and what it deliberately is not.
 *
 * Lensfun's per-frame cost in a raw pipeline is not its database: it is that the
 * displacement map can only be produced by its C++ callback machinery, single-threaded,
 * on the CPU — measured at 278 ms per 24 Mpx frame — after which a GPU pipeline must
 * upload six floats per pixel just to resample. Yet the mathematics is six small closed
 * forms: a census of the complete lensfun database (version_1, 2026) finds distortion
 * ptlens/poly3/poly5 (4810/875/5 entries), TCA poly3/linear (3355/6), vignetting pa
 * (25269, the only model). LensSerious expresses those forms as plain C over a plain
 * struct of coefficients, so the same evaluation can run vectorised on the CPU, inside
 * an OpenCL kernel (opencl/lensserious.cl is the same math, textually), or anywhere a
 * float goes.
 *
 * The lensfun PROJECT is not being forked: its XML database remains the interchange
 * format and its community remains the source of calibrations. What is replaced is the
 * runtime. Every convention below is ported from lensfun 0.3.4 source with its quirks
 * kept on purpose — the (width-1) sizing, the calibration-sensor aspect/crop correction,
 * the Hermite interpolation with one-sided tangents — because a "cleaner" convention
 * would silently shift every render. Parity with liblensfun is not assumed; it is
 * asserted, by tests/parity_lensfun.c against the whole installed database, and by the
 * side-by-side latch in Ansel's iop/lens.cc which logs the live deviation per frame.
 *
 * Scope of this seed: distortion + TCA + linear scaling composed into the geometry map,
 * and pa vignetting. Projection conversion between lens types (fisheye → rectilinear
 * and friends) is not implemented yet: the modifier reports it via
 * ls_modifier_t::geometry_unsupported so a caller can fall back, and the parity harness
 * skips those lenses explicitly rather than pretending.
 */

#ifndef LENSSERIOUS_H
#define LENSSERIOUS_H

#include <stddef.h>

/* The closed forms themselves, shared verbatim with opencl/lensserious.cl. */
#include "lensserious_eval.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum ls_dist_model_t
{
  LS_DIST_NONE = 0,
  LS_DIST_POLY3,   /**< Rd = Ru · (1 − k1 + k1·Ru²)            terms: k1 */
  LS_DIST_POLY5,   /**< Rd = Ru · (1 + k1·Ru² + k2·Ru⁴)        terms: k1 k2 */
  LS_DIST_PTLENS,  /**< Rd = Ru · (a·Ru³ + b·Ru² + c·Ru + d),  d = 1−a−b−c  terms: a b c */
  LS_DIST_KNOTS,   /**< Rd = Ru · cor(Ru), cor read from a table. No terms; see ls_knots_t. */
} ls_dist_model_t;

typedef enum ls_tca_model_t
{
  LS_TCA_NONE = 0,
  LS_TCA_LINEAR,   /**< Rd = Ru · k per channel                 terms: kr kb */
  LS_TCA_POLY3,    /**< Rd = Ru · (b·Ru² + c·Ru + v) per channel terms: vr vb cr cb br bb */
} ls_tca_model_t;

typedef enum ls_vig_model_t
{
  LS_VIG_NONE = 0,
  LS_VIG_PA,       /**< Cd = Cs · (1 + k1·r² + k2·r⁴ + k3·r⁶)  terms: k1 k2 k3 */
  LS_VIG_KNOTS,    /**< Cd = Cs · v(r), v read from a table. No terms; see ls_knots_t. */
} ls_vig_model_t;

/** One calibration entry, at one focal length (vignetting: one focal/aperture/distance). */
typedef struct ls_calib_dist_t { ls_dist_model_t model; float focal; float terms[3]; } ls_calib_dist_t;
typedef struct ls_calib_tca_t  { ls_tca_model_t  model; float focal; float terms[6]; } ls_calib_tca_t;
typedef struct ls_calib_vig_t  { ls_vig_model_t  model; float focal, aperture, distance; float terms[3]; } ls_calib_vig_t;
/** One <real-focal-length> point: the focal engraved on the barrel, and the one the lens
 * actually has there. They differ by up to a factor of two on fisheyes. */
typedef struct ls_calib_real_focal_t { float focal; float real_focal; } ls_calib_real_focal_t;

enum { LS_MAX_CALIB = 512 }; /* the densest lens in the 2026 database has ~300 vignetting points */

typedef enum ls_lens_type_t
{
  LS_LENS_UNKNOWN = 0,
  LS_LENS_RECTILINEAR,
  LS_LENS_FISHEYE,
  LS_LENS_PANORAMIC,
  LS_LENS_EQUIRECTANGULAR,
  LS_LENS_FISHEYE_ORTHOGRAPHIC,
  LS_LENS_FISHEYE_STEREOGRAPHIC,
  LS_LENS_FISHEYE_EQUISOLID,
  LS_LENS_FISHEYE_THOBY,
} ls_lens_type_t;

/** A lens as data: coefficients plus the calibration sensor's identity. */
typedef struct ls_lens_t
{
  ls_lens_type_t type;
  float crop_factor;   /**< of the CALIBRATION sensor: coefficient rescaling depends on it */
  float aspect_ratio;  /**< of the calibration sensor, e.g. 1.5 */
  float min_focal, max_focal; /**< the lens's whole range: the vignetting metric needs it */
  float center_x, center_y; /**< optical centre shift, lensfun convention (fraction of size) */
  int n_dist; ls_calib_dist_t dist[LS_MAX_CALIB];
  int n_tca;  ls_calib_tca_t  tca[LS_MAX_CALIB];
  int n_vig;  ls_calib_vig_t  vig[LS_MAX_CALIB];
  /** <real-focal-length>, which FEEDS the projection stage rather than gating it.
   *
   * lensfun's geometry callback runs on GetRealFocalLength(focal) divided by
   * get_hugin_focal_correction(focal), and those two cancel to exactly the nominal focal
   * when a lens carries no real-focal data -- the correction is multiplied in and divided
   * straight back out (modifier.cpp, and measured as a ratio of 1.0000 across every fisheye
   * in the database). When a lens DOES carry it, they do not cancel: the geometry focal is
   * real_focal / hugin, which on the Sigma 4.5mm circular fisheye is 0.4727x the nominal
   * one. Correcting such a lens with its nominal focal is 28 px out at the centre of the
   * frame and 135 px at the edge, which is why these points have to be carried. */
  int n_real_focal; ls_calib_real_focal_t real_focal[LS_MAX_CALIB];
} ls_lens_t;

#define LS_ENABLE_DISTORTION (1 << 0)
#define LS_ENABLE_TCA        (1 << 1)
#define LS_ENABLE_VIGNETTING (1 << 2)
#define LS_ENABLE_SCALE      (1 << 3)
#define LS_ENABLE_GEOMETRY   (1 << 4)

/**
 * @brief A lens correction the camera maker measured and wrote into the file, as knots.
 *
 * @details The second source of correction data this library resolves, beside the lens
 * database. Most mirrorless makers -- Sony, Fujifilm, Olympus, and anyone writing DNG
 * 1.3 opcodes -- embed a profile of the lens that took the picture in its metadata, and
 * they all express it the same way: a short list of radii and, at each, the factor the
 * coordinate should be scaled by. Not a model with coefficients; the samples themselves,
 * meant to be read back with straight-line interpolation between them.
 *
 * Two things about the shape of it matter to everything downstream:
 *
 * - It is per CHANNEL. `cor_rgb` holds red, green and blue separately, so distortion
 *   and lateral chromatic aberration arrive measured together rather than as a geometry
 *   plus a correction on top of it. There is no TCA stage to run afterwards.
 * - It is measured on the lens as shipped. There is no projection to change and no
 *   calibration sensor to convert from, which is why this needs none of the arguments
 *   ls_modifier_init() takes: no crop factor, no aperture, no subject distance, and no
 *   focal, because the maker already resolved the profile at the focal that was used.
 *
 * `radius`, and `vig_radius` with it, are the distance from the IMAGE CENTRE over
 * half the image diagonal -- 1.0 at the far corner. That is the makers' own convention and
 * this library adopts it wholesale rather than converting, so the numbers a decoder lifts
 * out of the file go in unmodified.
 */
typedef struct ls_knots_t
{
  /** How many distortion/TCA knots, up to LS_MAX_KNOTS. Zero for none. */
  int   n;
  /** Ascending radii, #n of them, shared by all three channels. */
  float radius[LS_MAX_KNOTS];
  /** At each radius, source_radius / this_radius, per channel: 0 red, 1 green, 2 blue.
   * The direction is the CORRECTING one -- given a point in the corrected image, where in
   * the source image it came from -- which is what the makers publish and what a resampler
   * consumes. Ask ls_modifier_init_knots() for the other one and it inverts the table. */
  float cor_rgb[3][LS_MAX_KNOTS];

  /** How many vignetting knots. Zero for none; independent of #n. */
  int   vn;
  /** Ascending radii for #vig, in the same units as #radius. */
  float vig_radius[LS_MAX_KNOTS];
  /** The falloff itself, as the maker measured it: below 1 where the lens darkens. A pixel
   * is CORRECTED by dividing by it, which is upstream lensfun's convention too. */
  float vig[LS_MAX_KNOTS];
} ls_knots_t;

/**
 * @brief A modifier: the lens resolved at one shooting configuration.
 *
 * @details Everything here is a plain value: build it once per commit, hand copies to
 * threads, embed the resolved terms in a kernel's argument block. Nothing is owned,
 * nothing is freed, nothing is locked.
 */
typedef struct ls_modifier_t
{
  /* Coordinate system — ported from lfModifier::lfModifier() (modifier.cpp).
   * norm_scale converts centred pixel coordinates to the calibration sensor's
   * normalized system; centre includes the lens's optical-centre shift. */
  float norm_scale, norm_unscale;
  float aspect_ratio_correction; /* vignetting works in the half-diagonal system */
  float center_x, center_y;   /* in normalized coordinates */
  float width, height;        /* the (dimension − 1) values, lensfun's own convention */

  int   enabled;              /* LS_ENABLE_* actually resolved (calibration present) */
  int   reverse;              /* non-zero for the reverse direction (lensfun's `reverse`) */
  /* Set when the requested projection change cannot be expressed radially -- panoramic or
   * equirectangular on either side. Those two map x and y differently; everything else
   * (rectilinear, the four fisheyes, thoby) goes through ls_eval_geometry(). */
  int   geometry_unsupported;

  float scale;                /* linear scaling factor, applied first (priority 100) */
  int   geom_from, geom_to;   /* projection change, ls_lens_type_t (priority 500) */
  float geom_focal;           /* focal / 12mm, the normalized-radius reference */
  ls_calib_dist_t dist;       /* resolved (interpolated) at the shooting focal */
  ls_calib_tca_t  tca;
  ls_calib_vig_t  vig;

  /* The alternative to the three above: a maker's own measured table, already turned to
   * face the requested direction. Filled by ls_modifier_init_knots(), left zeroed by
   * ls_modifier_init(), and never both.
   *
   * Shaped the way ls_eval_t holds it rather than the way ls_knots_t states it, because
   * turning it round is what this struct is for and doing it twice is how the two would
   * drift. The radii are per channel HERE and shared THERE: inverting sends each channel to
   * its own set of radii, so only the direction that needs them can carry them. */
  int   knot_n, knot_vn;
  float knot_r[3][LS_MAX_KNOTS];
  float knot_c[3][LS_MAX_KNOTS];
  float knot_vr[LS_MAX_KNOTS];
  float knot_v[LS_MAX_KNOTS];
} ls_modifier_t;

/**
 * @brief Resolve a lens at one (crop, geometry, focal, aperture, distance, scale).
 *
 * @param mod filled in by this call; nothing in it is owned or must be freed.
 * @param lens the lens, as data.
 * @param crop the crop factor of the sensor the picture was TAKEN with.
 * @param width, height the image dimensions, in pixels.
 * @param focal, aperture, distance the shooting configuration.
 * @param scale a linear scaling factor; 1.0 for none.
 * @param target_type the projection the output should be in, as ls_lens_type_t. Pass the
 * lens's own type (or LS_LENS_UNKNOWN) for no projection change; LS_ENABLE_GEOMETRY is
 * only raised when it actually differs and the pair is radially expressible.
 * @param flags which LS_ENABLE_* axes to attempt.
 * @param reverse non-zero to resolve the REVERSE direction -- the transform that takes a
 * corrected coordinate back to where it came from in the source image, which is what a
 * consumer needs to place masks and drawn shapes on an image it is correcting. It is not
 * the same chain read backwards: upstream registers different callbacks at different
 * priorities, so the composition order changes too (see ls_eval_t::reverse). An axis whose
 * model cannot be inverted at these coefficients -- poly3 distortion with k1 = 0, linear
 * TCA with a zero term -- is dropped from the returned flags exactly as upstream drops it.
 * @return the LS_ENABLE_* flags actually in effect, mirroring lensfun's oflags.
 */
int ls_modifier_init(ls_modifier_t *mod, const ls_lens_t *lens,
                     float crop, int width, int height,
                     float focal, float aperture, float distance,
                     float scale, int target_type, int flags, int reverse);

/**
 * @brief Resolve a maker's embedded profile, in place of a database lens.
 *
 * @param mod filled in by this call; nothing in it is owned or must be freed.
 * @param knots the table, as lifted from the file's metadata. Copied, not retained.
 * @param width, height the image dimensions, in pixels.
 * @param scale a linear scaling factor; 1.0 for none.
 * @param flags which LS_ENABLE_* axes to attempt. LS_ENABLE_TCA has no meaning here and is
 * ignored: a maker's table is already per channel, so asking for distortion asks for the
 * chromatic part of it too. LS_ENABLE_GEOMETRY likewise -- there is no projection to change.
 * @param reverse non-zero for the reverse direction, as in ls_modifier_init().
 * @return the LS_ENABLE_* flags actually in effect. An axis with no knots is dropped.
 *
 * @details The counterpart to ls_modifier_init(), and deliberately the only difference
 * between the two paths: what comes out is an ls_modifier_t like any other, which
 * ls_eval_from_modifier() flattens like any other, which the same evaluator and the same
 * kernels consume like any other. A consumer that already corrects lenses from the database
 * gains embedded profiles by choosing a different resolver -- not by growing a second pixel
 * path, and not by touching the one it has.
 *
 * That is possible because the models agree on what a correction IS. Every one of them,
 * polynomial or tabulated, answers the same question: at this radius, by what factor is the
 * coordinate scaled. So the table takes the distortion slot in the chain and nothing else
 * moves -- @p scale still composes around it in priority order, and the caller can still
 * ask ls_modifier_autoscale() what scale removes the borders, because that measures the
 * chain rather than the model.
 *
 * @note @p reverse costs nothing per pixel, unlike the polynomial models. Inverting
 * r -> r*cor(r) has no closed form for a polynomial, so those models pay a Newton iteration
 * at every pixel; a piecewise-linear curve inverts by reading the same points the other way
 * round -- (r*cor(r), 1/cor(r)) -- once, at this call, into a table the same size.
 *
 * That inverse is exact AT the knots and second-order between them, since a segment that is
 * straight going forwards is not straight coming back. Measured on a nine-knot profile with
 * 2% of distortion at the corner, over a 6000x4000 frame: correcting a point and then
 * un-correcting it returns it to within **0.13 px** at worst, and the worst case sits at
 * r = 0.94, mid-segment, exactly where the reasoning says it should. For scale, that is a
 * sixth of the 0.8 px by which liblensfun's own vectorised row walk differs from
 * liblensfun's own scalar answer. tests/knots.c is where the number comes from and will
 * fail if it grows past a quarter of a pixel.
 *
 * @note The three channels' inverses land on three different sets of radii, which is why
 * ls_eval_t carries a radius axis per channel where this input carries one shared. Forcing
 * them back onto a common axis would mean resampling two of the three, and there is nothing
 * to gain by it -- the tables are already resident in the block either way.
 */
int ls_modifier_init_knots(ls_modifier_t *mod, const ls_knots_t *knots,
                           int width, int height, float scale, int flags, int reverse);

/**
 * @brief The scale that just removes the black borders a correction leaves behind.
 *
 * @param mod a resolved modifier. Its own scale factor, if it has one, is part of the
 * transform being measured -- exactly as upstream measures whatever callbacks are
 * registered at the time.
 * @return the linear factor to feed back as ls_modifier_init()'s @p scale, or 1.0 when the
 * modifier transforms no coordinates at all.
 *
 * @details Ported from lfModifier::GetAutoScale(). Eight points around the frame -- the
 * four edge midpoints and the four corners -- are pushed through the coordinate chain, and
 * the largest ratio of where a point should be to where it landed is the scale. It is a
 * measurement, not a formula: there is no closed form for "where does the corner of a
 * distorted, reprojected frame end up", so each point is found by Newton iteration on the
 * chain itself with a numeric derivative.
 *
 * Two upstream constants are kept rather than tidied, because they are what its renders
 * were produced with: a flat 1.001 margin ("1 permille is our limit of accuracy"), and a
 * second 1.001 when TCA is active, since the per-channel radii extend slightly past the
 * green one that was measured.
 *
 * @note Where a point cannot be found -- ultrawide fisheye corners extending to infinity --
 * upstream lets the iteration fail and the point simply loses the max(). Same here.
 */
float ls_modifier_autoscale(const ls_modifier_t *mod);

/**
 * @brief Flatten a resolved modifier into the scalar block a kernel can take by value.
 *
 * @details This is the whole point of the exercise: a correction crosses to the GPU as one
 * fixed 632-byte block that every work-item evaluates for itself, instead of a
 * six-float-per-pixel map that the CPU builds single-threaded and then uploads (measured
 * at 278 ms plus 576 MB of transfer for a 24 Mpx frame). Fixed regardless of which kind of
 * correction it carries -- a maker's table is a hundred-odd more floats, and the whole
 * block is still smaller than one row of that map.
 *
 * @return 0 if either pointer is NULL, 1 otherwise.
 */
int ls_eval_from_modifier(const ls_modifier_t *mod, ls_eval_t *out);

/**
 * @brief Move @p src's vignetting into @p dst, leaving @p dst's geometry untouched.
 *
 * @details For a consumer whose axes do not all come from the same place -- an Olympus body
 * embeds distortion and lateral CA but no vignetting, so the falloff has to come from the
 * database while the geometry comes from the file. Resolve one modifier per source, flatten
 * both, and graft the vignetting across.
 *
 * It is a graft rather than a merge because the two halves of ::ls_eval_t are genuinely
 * independent: vignetting reads vig_scale, vig_center_x/y, vig_model and either vig_terms or
 * the vignetting knot table, and NOTHING a coordinate transform touches. In particular each
 * half carries its own normalization, which is what makes this safe across resolvers that do
 * not share one -- the table path measures radius against the half diagonal, the database
 * path against lensfun's short side, and neither has to know about the other.
 *
 * Cheaper than the alternative, too: a second block would double what a kernel receives by
 * value, and ::ls_eval_t is already 632 bytes against the 1024 that OpenCL 1.2 guarantees
 * for a kernel's whole argument list.
 *
 * @return 0 if either pointer is NULL, 1 otherwise.
 */
int ls_eval_adopt_vignetting(ls_eval_t *dst, const ls_eval_t *src);

/**
 * @brief The geometry map: for @p count output pixels starting at (@p xu, @p yu),
 * write 6 floats per pixel — source coordinates for R, G, B — in pixel space.
 * Bit-for-bit the contract of lfModifier::ApplySubpixelGeometryDistortion(), so a
 * caller can compare the two buffers element-wise.
 */
int ls_modifier_apply_subpixel_geometry(const ls_modifier_t *mod,
                                        float xu, float yu, int width, int height,
                                        float *res);

/**
 * @brief pa vignetting, multiplied in place over RGBA float rows.
 * Contract of lfModifier::ApplyColorModification(LF_CR_4(RED,GREEN,BLUE,UNKNOWN), F32).
 */
int ls_modifier_apply_vignetting(const ls_modifier_t *mod,
                                 float xu, float yu, int width, int height,
                                 float *rgba, int row_stride_bytes);

#ifdef __cplusplus
}
#endif

#endif /* LENSSERIOUS_H */
