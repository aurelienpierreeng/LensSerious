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
} ls_vig_model_t;

/** One calibration entry, at one focal length (vignetting: one focal/aperture/distance). */
typedef struct ls_calib_dist_t { ls_dist_model_t model; float focal; float terms[3]; } ls_calib_dist_t;
typedef struct ls_calib_tca_t  { ls_tca_model_t  model; float focal; float terms[6]; } ls_calib_tca_t;
typedef struct ls_calib_vig_t  { ls_vig_model_t  model; float focal, aperture, distance; float terms[3]; } ls_calib_vig_t;

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
} ls_lens_t;

#define LS_ENABLE_DISTORTION (1 << 0)
#define LS_ENABLE_TCA        (1 << 1)
#define LS_ENABLE_VIGNETTING (1 << 2)
#define LS_ENABLE_SCALE      (1 << 3)
#define LS_ENABLE_GEOMETRY   (1 << 4)

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
} ls_modifier_t;

/**
 * @brief Resolve a lens at one (crop, geometry, focal, aperture, distance, scale).
 *
 * @param target_type the projection the output should be in, as ls_lens_type_t. Pass the
 * lens's own type (or LS_LENS_UNKNOWN) for no projection change; LS_ENABLE_GEOMETRY is
 * only raised when it actually differs and the pair is radially expressible.
 * @return the LS_ENABLE_* flags actually in effect, mirroring lensfun's oflags.
 */
int ls_modifier_init(ls_modifier_t *mod, const ls_lens_t *lens,
                     float crop, int width, int height,
                     float focal, float aperture, float distance,
                     float scale, int target_type, int flags);

/**
 * @brief Flatten a resolved modifier into the scalar block a kernel can take by value.
 *
 * @details This is the whole point of the exercise: a correction crosses to the GPU as
 * ~80 bytes of coefficients that every work-item evaluates for itself, instead of a
 * six-float-per-pixel map that the CPU builds single-threaded and then uploads (measured
 * at 278 ms plus 576 MB of transfer for a 24 Mpx frame).
 *
 * @return 0 if either pointer is NULL, 1 otherwise.
 */
int ls_eval_from_modifier(const ls_modifier_t *mod, ls_eval_t *out);

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
