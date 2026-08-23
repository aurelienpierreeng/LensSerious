/*
    LensSerious — a maker's own correction tables, decoded into the library's generic knots.

    Copyright (C) 2026 Aurélien PIERRE.
    Ported from darktable's embedded-metadata lens correction (src/iop/lens.cc and
    src/common/dng_opcode.c), Copyright (C) the darktable contributors.

    License: GPL-3.0-or-later — NOT the library's usual LGPL. The decoders derive from
    GPL-licensed darktable code whose copyright this project does not own, so this module
    is a separate build target (`lensserious_vendor`) a consumer links only if the GPL
    suits it. Nothing in the core library depends on anything here.
*/

/** @file lensserious_vendor.h
 *
 * @brief Turn a camera maker's embedded lens-correction tables into #ls_knots_t.
 *
 * @details Some cameras write their own lens profile into every raw file — measured on the
 * actual body and the actual lens, which makes it the best-informed correction available.
 * All the formats this module reads express it the same way underneath: a short list of
 * radii and, at each, the factor a coordinate is scaled by — per channel, so distortion
 * and lateral chromatic aberration arrive measured together. This module normalises
 * whichever format the file holds into one #ls_knots_t, so past ls_vendor_resolve()
 * nothing knows or cares which maker wrote the profile.
 *
 * @section vendor_boundary Where this API sits — after metadata decoding
 *
 * This library links no metadata reader, deliberately. The boundary is drawn where the
 * *numbers* appear: the consumer's EXIF library (exiv2, libraw, anything) reads the
 * maker-note tags and fills one of the per-vendor structs below with the decoded values;
 * this module owns everything after that — the vendor-specific coefficient semantics, the
 * per-channel untangling, the resampling onto a common radius axis. Re-parsing maker-note
 * *bytes* here would mean re-doing tag layout and endianness that the consumer's metadata
 * library already resolved, so the input structs carry values, not bytes.
 *
 * The one exception proves the rule: a DNG's OpcodeList3 payload is an opaque,
 * self-contained big-endian blob that no metadata library decodes for you. Its format is
 * correction-domain knowledge, so ls_vendor_parse_dng_opcodelist3() lives here — the
 * consumer hands the raw tag payload and gets the filled struct back.
 *
 * @section vendor_extraction What the consumer extracts, per vendor
 *
 * Worked examples from Ansel's `src/metadata/exif.cc` (exiv2), for future implementors.
 * Tag spellings are exiv2's; other libraries name the same tags differently.
 *
 * **Sony** — three SubIFD arrays of int16, count first:
 * `Exif.SubImage1.DistortionCorrParams`, `.ChromaticAberrationCorrParams` (2×nc values:
 * red then blue), `.VignettingCorrParams` (fallback spellings `Exif.Image.0x7037` /
 * `0x7035` / `0x7032` on older bodies). Validate `nc` in [2, 16], that the CA array holds
 * 2×nc, then copy values 1..nc of each into #ls_vendor_sony_t verbatim. The radii are
 * implicit (evenly spaced); the fixed-point scale factors are this module's business, not
 * the extractor's.
 *
 * **Fujifilm** — three maker-note float arrays sharing one knot list:
 * `Exif.Fujifilm.GeometricDistortionParams` (19 or 23 values: 9 or 11 knots then as many
 * coefficients), `.ChromaticAberrationParams` (29 or 31: knots, red, blue),
 * `.VignettingParams`. The extractor checks the three knot lists agree and fills
 * #ls_vendor_fuji_t with knots and coefficients separately.
 *
 * @anchor vendor_fuji_cropf
 * `cropf` is the extractor's responsibility, by contract. It is CAMERA knowledge, not
 * format knowledge: the tables are calibrated against the full sensor, so when the body
 * shot in a cropped mode the knot radii must be scaled by the crop ratio — and whether
 * the file was shot cropped is written nowhere near the correction tables. Ansel decides
 * it from `Exif.Fujifilm.CropMode`: values 2 ("sports finder") and 4 (electronic
 * teleconverter) mean the 1.25× crop, so
 * @code
 *   fuji.cropf = (crop_mode == 2 || crop_mode == 4) ? 1.25f : 1.0f;
 * @endcode
 * An extractor that cannot determine the mode must pass 1.0f — never 0, which would
 * collapse every radius to the image centre. This module applies `cropf` blindly and has
 * no opinion about where it came from.
 *
 * **DNG** — one opaque byte blob: `Exif.SubImage1.OpcodeList3` (or `Exif.Image.OpcodeList3`
 * on files without an embedded preview). Hand it, undecoded, to
 * ls_vendor_parse_dng_opcodelist3(). A WarpRectilinear opcode with one plane is
 * distortion alone; with three planes it is distortion and TCA measured together.
 * VignetteRadial is the falloff. Note OpcodeList*2* is a different list (sensor-domain
 * corrections, gain maps) and is none of this module's business.
 *
 * **Olympus** — two maker-note arrays of already-decoded rationals/floats:
 * `Exif.OlympusIp.0x150c` (distortion: k2, k4, k6, scale) and `Exif.OlympusIp.0x150e`
 * (CA: three coefficients each for red and blue). Olympus publishes no vignetting;
 * that axis is simply absent from what ls_vendor_axes() reports, and a consumer offering
 * per-axis sources falls back to its other ones (Ansel uses its lens database there).
 *
 * @section vendor_conventions Output conventions
 *
 * The emitted #ls_knots_t follows the library's knot conventions (see its own
 * documentation): radii normalized to half the image diagonal, per-channel correcting
 * scale factors where green IS the distortion and red/blue carry the CA relative to it,
 * vignetting as the measured falloff (below 1 where the lens darkens, corrected by
 * dividing). The autoscale returned alongside is the zoom that just clears the borders
 * the geometric correction leaves — apply it as (or composed into) the `scale` argument
 * of ls_modifier_init_knots(); it is not baked into the table beyond the normalisation
 * the resolver already performed.
 *
 * @section vendor_limits Known limits
 *
 * The DNG opcodes carry an optical-centre offset (`warp_cx/cy`, `vig_cx/cy`); the
 * conversion currently assumes a centred system and ignores them — wrong by the offset
 * amount for the rare decentred profile, faithful to what every shipping consumer of this
 * code has done so far. The struct carries the values so fixing it needs no ABI change.
 */

#ifndef LENSSERIOUS_VENDOR_H
#define LENSSERIOUS_VENDOR_H

#include <stddef.h>
#include <stdint.h>

#include "lensserious.h"

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Which maker's format a #ls_vendor_data_t holds. NONE == 0, so a zero-initialised
 *  struct safely means "no data". Values are not a wire format — do not persist them. */
typedef enum ls_vendor_type_t
{
  LS_VENDOR_NONE = 0,
  LS_VENDOR_SONY,
  LS_VENDOR_FUJI,
  LS_VENDOR_DNG,
  LS_VENDOR_OLYMPUS
} ls_vendor_type_t;

/** @brief Sony maker-note tables, verbatim: fixed-point int16 at #nc implicit, evenly
 *  spaced radii. The fixed-point scales (1/16384 distortion, 1/2^21 CA, 1/8192
 *  vignetting) are applied by ls_vendor_resolve(), never by the extractor. */
typedef struct ls_vendor_sony_t
{
  int nc;                 /**< knots; valid data has 2..16 */
  int16_t distortion[16];
  int16_t ca_r[16];
  int16_t ca_b[16];
  int16_t vignetting[16];
} ls_vendor_sony_t;

/** @brief Fujifilm maker-note tables, verbatim: explicit knot radii shared by the three
 *  corrections, percentage-domain coefficients. */
typedef struct ls_vendor_fuji_t
{
  int nc;            /**< knots; valid data has 1..11 (9 or 11 in the wild) */
  float cropf;       /**< crop-mode radius factor, decided by the EXTRACTOR — see
                          @ref vendor_fuji_cropf "the cropf contract". 1.0f when unknown,
                          never 0. */
  float knots[11];   /**< radii, half-diagonal normalized, ascending */
  float distortion[11]; /**< percent deviation at each knot */
  float ca_r[11];    /**< red scale offset at each knot */
  float ca_b[11];    /**< blue scale offset at each knot */
  float vignetting[11]; /**< percent transmission at each knot */
} ls_vendor_fuji_t;

/** @brief A DNG's OpcodeList3 corrections, parsed. Fill it with
 *  ls_vendor_parse_dng_opcodelist3() — or by hand, if the consumer's DNG reader already
 *  decoded the opcodes itself. */
typedef struct ls_vendor_dng_t
{
  int has_warp;             /**< a WarpRectilinear (opcode 1) was present */
  int has_vignette;         /**< a VignetteRadial (opcode 3) was present */
  uint32_t warp_planes;     /**< 1 = distortion only; >1 (typically 3) = per-channel, i.e. TCA too */
  double warp_coeffs[3][6]; /**< per plane: kr0..kr3 radial, kt0..kt1 tangential */
  double warp_cx;           /**< normalized optical centre [0..1] — carried, not yet consumed */
  double warp_cy;
  double vig_coeffs[5];     /**< k0..k4 radial vignette polynomial */
  double vig_cx;            /**< normalized optical centre [0..1] — carried, not yet consumed */
  double vig_cy;
} ls_vendor_dng_t;

/** @brief Olympus maker-note polynomials, verbatim. No vignetting: Olympus does not
 *  publish one, and no flag here can invent it. */
typedef struct ls_vendor_olympus_t
{
  int has_dist;
  float dist[4]; /**< k2, k4, k6, radial scale */
  int has_ca;
  float ca[6];   /**< red: c0, c2, c4; blue: c0, c2, c4 */
} ls_vendor_olympus_t;

/** @brief One image's embedded correction data: the type tag and the matching member.
 *  A plain value — embed it in an image struct, memcpy it, zero it for "none". */
typedef struct ls_vendor_data_t
{
  ls_vendor_type_t type;
  union
  {
    ls_vendor_sony_t sony;
    ls_vendor_fuji_t fuji;
    ls_vendor_dng_t dng;
    ls_vendor_olympus_t olympus;
  } u;
} ls_vendor_data_t;

/** @brief Per-axis strength blend, 1.0 = exactly as the maker measured. The vendors' own
 *  GUIs expose these as "fine-tune" sliders; a consumer that does not can simply pass
 *  NULL to ls_vendor_resolve(), which means all-1.0. */
typedef struct ls_vendor_finetune_t
{
  float distortion;
  float vignette;
  float ca_red;
  float ca_blue;
} ls_vendor_finetune_t;

/** @brief Resolve/parse result: the data is structurally malformed. Distinct from 0
 *  ("carries nothing") so the consumer can SAY SO — a user whose file claims a profile
 *  that cannot be read is entitled to hear why the fallback engaged. */
#define LS_VENDOR_EBADDATA (-1)

/**
 * @brief Which correction axes this data carries.
 *
 * @return A mask of LS_ENABLE_DISTORTION, LS_ENABLE_TCA and LS_ENABLE_VIGNETTING; 0 when
 * the data carries nothing usable (including @p d NULL or type NONE).
 *
 * Cheap — table scans, no knots work — so a GUI may ask per repaint. The same mask
 * vocabulary as ls_modifier_init() returns, one axis language across the library.
 * An all-zero table is reported as ABSENT, not as a present-but-identity correction:
 * several bodies write zero-filled tables for lenses they choose not to correct.
 */
int ls_vendor_axes(const ls_vendor_data_t *d);

/**
 * @brief The conversion: vendor tables in, generic knots out.
 *
 * @param d the extracted data. See @ref vendor_extraction "what the consumer extracts".
 * @param ft per-axis strength, or NULL for "as the maker measured" (all 1.0).
 * @param w image width in pixels — the PROCESSED dimensions the radii are measured
 * against (after the sensor's fixed border trim, before any user crop). Ansel passes
 * `dt_image_t.p_width`/`p_height`.
 * @param h image height in pixels.
 * @param knots out: the normalized table, ready for ls_modifier_init_knots(). Fully
 * written on success (n and vn included); zeroed on failure.
 * @param autoscale out, may be NULL: the zoom that just clears the borders the correction
 * leaves, measured over the region between the inscribed circle and the corner. Pass it
 * (times any user scale) as ls_modifier_init_knots()'s scale argument. 1.0 on failure.
 *
 * @return The resolved axes mask (same bits as ls_vendor_axes()); 0 when the data carries
 * nothing — a zero-filled table with a valid type included, see ls_vendor_axes();
 * #LS_VENDOR_EBADDATA when a claimed correction fails to convert.
 *
 * Pure and thread-safe: no allocation, no globals, callable from a pipeline thread per
 * commit. Not per pixel — resolve once, keep the knots.
 */
int ls_vendor_resolve(const ls_vendor_data_t *d, const ls_vendor_finetune_t *ft,
                      int w, int h, ls_knots_t *knots, float *autoscale);

/**
 * @brief Parse a DNG OpcodeList3 payload: WarpRectilinear (1) and VignetteRadial (3).
 *
 * @param blob the tag's raw bytes, exactly as the file carries them (big-endian,
 * per the DNG 1.3 specification). Not written to.
 * @param len the payload size in bytes.
 * @param out filled with whatever was recognised; zeroed first, so a partial parse
 * leaves the untouched corrections cleanly absent.
 *
 * @return The axes mask the parsed opcodes provide; 0 when the list holds none of the
 * recognised opcodes (other opcode ids are skipped, as the specification intends);
 * #LS_VENDOR_EBADDATA when the list structure itself is broken before anything was
 * recognised. A list that breaks AFTER yielding a usable opcode returns what it yielded:
 * partial data is data.
 *
 * Every read is bounds-checked against @p len; a truncated or hostile blob cannot read
 * past the buffer.
 */
int ls_vendor_parse_dng_opcodelist3(const uint8_t *blob, size_t len, ls_vendor_dng_t *out);

#ifdef __cplusplus
}
#endif

#endif /* LENSSERIOUS_VENDOR_H */
