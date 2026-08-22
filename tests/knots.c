/*
    LensSerious — the maker-profile path, held to properties rather than to an oracle.

    Copyright (C) 2026 Aurélien PIERRE.  License: LGPL-3.0-or-later.
*/

/** @file knots.c
 *
 * @brief What can be asserted about a tabulated correction when nothing else implements it.
 *
 * Every other test here has liblensfun to compare against. This one does not: lensfun has
 * no knot model, so there is no second implementation to disagree with. That does not leave
 * the path unverifiable, it changes what is verified -- from "the same numbers as upstream"
 * to properties the correction has to satisfy whatever the numbers are:
 *
 *  1. **A table of ones is an identity.** Every pixel must come back exactly where it
 *     started, in both directions, with the table active rather than switched off. This is
 *     what catches an off-by-one in the radius normalisation, which no amount of staring at
 *     a distorted render will show.
 *  2. **The inverse is an inverse.** Correct a point, then un-correct it, and measure how
 *     far it landed from where it began -- IN PIXELS, because that is the unit the error
 *     matters in. The doc comment on ls_modifier_init_knots() claims the inversion is exact
 *     at the knots and second-order between them; this is the measurement behind the claim.
 *  3. **The channels separate by the amount the table says.** A profile whose red curve
 *     differs from its green one must move red by exactly that difference, and by nothing
 *     else -- the failure mode being a TCA stage that quietly also runs.
 *  4. **Refusal beats a wrong answer.** A table that cannot be inverted must come back with
 *     its axis dropped from the returned flags, not inverted approximately.
 *  5. **The rest of the chain still composes.** Scale, and the autoscale that measures it,
 *     work on a knot modifier because they operate on the chain rather than on the model.
 *  6. **Vignetting divides, and reversing multiplies back.**
 *
 * The profile used throughout is shaped like the ones cameras actually write: nine knots,
 * about 2% of barrel at the corner, red and blue a couple of tenths of a percent off green.
 */

#include "lensserious.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

static int failures = 0;

#define CHECK(cond, ...)                                                                   \
  do {                                                                                     \
    if(!(cond)) { printf("FAIL %s:%d: ", __FILE__, __LINE__); printf(__VA_ARGS__);          \
                  printf("\n"); failures++; }                                              \
  } while(0)

#define IMG_W 6000
#define IMG_H 4000

/* A plausible mirrorless wide-angle: barrel, so a point at the corner of the corrected
 * image came from further IN, and cor < 1 out there. Red bends slightly more than green and
 * blue slightly less, which is the ordinary sign of lateral chromatic aberration. */
static void _make_profile(ls_knots_t *k)
{
  memset(k, 0, sizeof(*k));

  static const float r[9]   = { 0.00f, 0.125f, 0.25f, 0.375f, 0.50f, 0.625f, 0.75f, 0.875f, 1.00f };
  static const float cg[9]  = { 1.0000f, 1.0004f, 1.0012f, 1.0018f, 1.0010f,
                                0.9985f, 0.9940f, 0.9875f, 0.9790f };
  /* The chromatic offset grows with radius, as it does optically: zero on axis. */
  k->n = 9;
  for(int i = 0; i < 9; i++)
  {
    const float chroma = 0.0025f * r[i];
    k->radius[i] = r[i];
    k->cor_rgb[0][i] = cg[i] * (1.f + chroma);
    k->cor_rgb[1][i] = cg[i];
    k->cor_rgb[2][i] = cg[i] * (1.f - chroma);
  }

  static const float vr[7] = { 0.00f, 0.20f, 0.40f, 0.55f, 0.70f, 0.85f, 1.00f };
  static const float vv[7] = { 1.000f, 0.981f, 0.930f, 0.878f, 0.812f, 0.734f, 0.648f };
  k->vn = 7;
  for(int i = 0; i < 7; i++) { k->vig_radius[i] = vr[i]; k->vig[i] = vv[i]; }
}

/* An identity table: the same nine radii, every factor exactly 1. */
static void _make_identity(ls_knots_t *k)
{
  memset(k, 0, sizeof(*k));
  k->n = 9;
  for(int i = 0; i < 9; i++)
  {
    k->radius[i] = (float)i / 8.f;
    for(int c = 0; c < 3; c++) k->cor_rgb[c][i] = 1.f;
  }
  k->vn = 5;
  for(int i = 0; i < 5; i++) { k->vig_radius[i] = (float)i / 4.f; k->vig[i] = 1.f; }
}

/* A 5x5 lattice over the frame, corners included, as the sample set for every property. */
#define GRID 5
static void _sample(const int i, const int j, float *x, float *y)
{
  *x = (float)i * (float)(IMG_W - 1) / (float)(GRID - 1);
  *y = (float)j * (float)(IMG_H - 1) / (float)(GRID - 1);
}

static int _resolve(ls_eval_t *p, const ls_knots_t *k, const float scale, const int flags,
                    const int reverse)
{
  ls_modifier_t mod;
  const int got = ls_modifier_init_knots(&mod, k, IMG_W, IMG_H, scale, flags, reverse);
  ls_eval_from_modifier(&mod, p);
  return got;
}

/* ------------------------------------------------------------------ */

static void test_identity(void)
{
  ls_knots_t k;
  _make_identity(&k);

  for(int rev = 0; rev < 2; rev++)
  {
    ls_eval_t p;
    const int got = _resolve(&p, &k, 1.f, LS_ENABLE_DISTORTION | LS_ENABLE_VIGNETTING, rev);
    CHECK(got & LS_ENABLE_DISTORTION, "identity: distortion not resolved (rev %d)", rev);
    CHECK(!(got & LS_ENABLE_TCA), "identity: TCA reported, and it never should be");

    float worst = 0.f;
    for(int i = 0; i < GRID; i++)
      for(int j = 0; j < GRID; j++)
      {
        float x, y, out[6];
        _sample(i, j, &x, &y);
        ls_eval_map(&p, x, y, out);
        for(int c = 0; c < 3; c++)
        {
          const float dx = out[2 * c] - x, dy = out[2 * c + 1] - y;
          const float d = sqrtf(dx * dx + dy * dy);
          if(d > worst) worst = d;
        }
      }
    /* Not a tolerance chosen to pass: the arithmetic really is exact here except for the
     * round trip through the normalized system, which costs a few ulps of 6000. */
    CHECK(worst < 1e-2f, "identity: moved a pixel by %.4f px (rev %d)", worst, rev);
    printf("knots: identity rev=%d worst displacement %.6f px\n", rev, worst);

    const float v = ls_eval_vignette_factor(&p, IMG_W * 0.5f, IMG_H * 0.5f);
    CHECK(fabsf(v - 1.f) < 1e-5f, "identity: vignetting factor %.6f, expected 1", v);
  }
}

static void test_roundtrip(void)
{
  ls_knots_t k;
  _make_profile(&k);

  ls_eval_t fwd, rev;
  const int gf = _resolve(&fwd, &k, 1.f, LS_ENABLE_DISTORTION, 0);
  const int gr = _resolve(&rev, &k, 1.f, LS_ENABLE_DISTORTION, 1);
  CHECK(gf & LS_ENABLE_DISTORTION, "roundtrip: forward not resolved");
  CHECK(gr & LS_ENABLE_DISTORTION, "roundtrip: reverse not resolved");

  /* A denser lattice than the property tests use: the inversion error is smallest at the
   * knots and largest between them, so a coarse grid would sample mostly the good case. */
  float worst = 0.f, worst_at = 0.f;
  for(int i = 0; i <= 40; i++)
    for(int j = 0; j <= 40; j++)
    {
      const float x = (float)i * (float)(IMG_W - 1) / 40.f;
      const float y = (float)j * (float)(IMG_H - 1) / 40.f;

      float a[6], b[6];
      ls_eval_map(&fwd, x, y, a);      /* corrected point -> where it came from */
      ls_eval_map(&rev, a[2], a[3], b); /* and back */
      const float dx = b[2] - x, dy = b[3] - y;
      const float d = sqrtf(dx * dx + dy * dy);
      if(d > worst)
      {
        worst = d;
        const float cx = x - IMG_W * 0.5f, cy = y - IMG_H * 0.5f;
        worst_at = sqrtf(cx * cx + cy * cy) / sqrtf(IMG_W * IMG_W + IMG_H * IMG_H) * 2.f;
      }
    }

  printf("knots: forward/reverse round trip, worst %.4f px (at r = %.3f)\n", worst, worst_at);
  /* The bound the documentation is allowed to claim. A quarter of a pixel is below what a
   * Mitchell resampler resolves, and below the 0.8 px liblensfun's own SSE row walk differs
   * from its own scalar answer by -- so this is not the weakest link in any render. */
  CHECK(worst < 0.25f, "roundtrip: %.4f px is too far to call the inverse an inverse", worst);
}

static void test_channels(void)
{
  ls_knots_t k;
  _make_profile(&k);

  ls_eval_t p;
  _resolve(&p, &k, 1.f, LS_ENABLE_DISTORTION, 0);

  /* At the corner the table says red is 0.25% further out than green. The map must
   * reproduce that ratio on the radius, and reproduce it on BOTH axes -- a per-channel
   * scale that leaked into only one of them would still look plausible on a test that
   * measured radius alone. */
  const float x = (float)(IMG_W - 1), y = (float)(IMG_H - 1);
  float out[6];
  ls_eval_map(&p, x, y, out);

  const float cx = IMG_W * 0.5f, cy = IMG_H * 0.5f;
  const float rr = hypotf(out[0] - cx, out[1] - cy);
  const float rg = hypotf(out[2] - cx, out[3] - cy);
  const float rb = hypotf(out[4] - cx, out[5] - cy);

  const float expect = 1.f + 0.0025f;   /* the chroma term at r = 1 */
  CHECK(fabsf(rr / rg - expect) < 2e-4f, "channels: R/G radius ratio %.6f, expected %.6f",
        rr / rg, expect);
  CHECK(fabsf(rb / rg - (2.f - expect)) < 2e-4f, "channels: B/G ratio %.6f, expected %.6f",
        rb / rg, 2.f - expect);
  CHECK(rr > rg && rg > rb, "channels: R %.2f G %.2f B %.2f are not in order", rr, rg, rb);

  /* Same scale on x and y: the displacement is radial, not axial. */
  const float sx = (out[0] - cx) / (out[2] - cx);
  const float sy = (out[1] - cy) / (out[3] - cy);
  CHECK(fabsf(sx - sy) < 1e-5f, "channels: R scaled x by %.6f and y by %.6f", sx, sy);
  printf("knots: corner radii R %.2f G %.2f B %.2f px\n", rr, rg, rb);
}

static void test_refusal(void)
{
  ls_knots_t k;
  _make_identity(&k);

  /* Out of order going in: the lookup walks the axis and would silently read the wrong
   * segment, so this has to be declined at resolve time. */
  ls_knots_t bad = k;
  bad.radius[4] = bad.radius[2];
  ls_eval_t p;
  int got = _resolve(&p, &bad, 1.f, LS_ENABLE_DISTORTION, 0);
  CHECK(!(got & LS_ENABLE_DISTORTION), "refusal: a non-ascending axis was accepted");

  /* A correction so strong the image folds back on itself: the forward axis is fine, and
   * only the CONSTRUCTED reverse one is not monotone. Forwards it must still be accepted. */
  ls_knots_t fold = k;
  for(int c = 0; c < 3; c++) { fold.cor_rgb[c][7] = 0.30f; fold.cor_rgb[c][8] = 0.10f; }
  got = _resolve(&p, &fold, 1.f, LS_ENABLE_DISTORTION, 0);
  CHECK(got & LS_ENABLE_DISTORTION, "refusal: a foldable profile was refused going forwards");
  got = _resolve(&p, &fold, 1.f, LS_ENABLE_DISTORTION, 1);
  CHECK(!(got & LS_ENABLE_DISTORTION), "refusal: an uninvertible profile was inverted anyway");

  /* A zero or negative factor is not a radius scale. */
  ls_knots_t zero = k;
  zero.cor_rgb[1][5] = 0.f;
  got = _resolve(&p, &zero, 1.f, LS_ENABLE_DISTORTION, 0);
  CHECK(!(got & LS_ENABLE_DISTORTION), "refusal: a zero factor was accepted");

  /* No knots at all is not an error, it is an absent axis. */
  ls_knots_t empty;
  memset(&empty, 0, sizeof(empty));
  got = _resolve(&p, &empty, 1.f, LS_ENABLE_DISTORTION | LS_ENABLE_VIGNETTING, 0);
  CHECK(got == 0, "refusal: an empty table resolved to flags %d", got);
}

static void test_chain(void)
{
  ls_knots_t k;
  _make_profile(&k);

  /* PINCUSHION, deliberately, and not the profile the other tests use: a barrel one leaves
   * no borders to remove when corrected -- its corners sample from further in -- so
   * autoscale would come back at the flat 1.001 margin and the test would pass without
   * having measured anything. Reversing the sign makes the corrected corner reach outside
   * the source frame, which is the case autoscale exists for. */
  for(int i = 0; i < k.n; i++)
    for(int c = 0; c < 3; c++)
      k.cor_rgb[c][i] = 2.f - k.cor_rgb[c][i];

  /* Autoscale has to work here for the same reason it works anywhere: it measures the
   * coordinate chain, and a knot table is a stage in that chain like any other. */
  ls_modifier_t mod;
  ls_modifier_init_knots(&mod, &k, IMG_W, IMG_H, 1.f, LS_ENABLE_DISTORTION, 0);
  const float s = ls_modifier_autoscale(&mod);
  CHECK(s > 1.01f && s < 1.2f,
        "chain: autoscale returned %.5f, which is not a plausible crop", s);
  printf("knots: autoscale %.5f\n", s);

  /* And feeding it back has to actually pull the corner in to the frame edge. */
  ls_eval_t p;
  const int got = _resolve(&p, &k, s, LS_ENABLE_DISTORTION | LS_ENABLE_SCALE, 0);
  CHECK(got & LS_ENABLE_SCALE, "chain: scale not resolved");

  float out[6];
  ls_eval_map(&p, (float)(IMG_W - 1), (float)(IMG_H - 1), out);
  CHECK(out[2] <= IMG_W && out[3] <= IMG_H,
        "chain: after autoscale the corner still samples (%.1f, %.1f)", out[2], out[3]);
  CHECK(out[2] > IMG_W * 0.9f, "chain: autoscale overshot, corner at %.1f", out[2]);
}

static void test_vignetting(void)
{
  ls_knots_t k;
  _make_profile(&k);

  ls_eval_t corr, appl;
  CHECK(_resolve(&corr, &k, 1.f, LS_ENABLE_VIGNETTING, 0) & LS_ENABLE_VIGNETTING,
        "vignetting: not resolved");
  _resolve(&appl, &k, 1.f, LS_ENABLE_VIGNETTING, 1);

  /* Correcting brightens where the table says the lens darkens, and the two directions
   * are reciprocal at every point. */
  for(int i = 0; i < GRID; i++)
    for(int j = 0; j < GRID; j++)
    {
      float x, y;
      _sample(i, j, &x, &y);
      const float a = ls_eval_vignette_factor(&corr, x, y);
      const float b = ls_eval_vignette_factor(&appl, x, y);
      CHECK(a >= 1.f - 1e-5f, "vignetting: correcting darkens to %.5f at (%.0f, %.0f)", a, x, y);
      CHECK(fabsf(a * b - 1.f) < 1e-4f, "vignetting: %.6f and %.6f are not reciprocal", a, b);
    }

  /* The corner is the 0.648 knot: correcting it is a lift of 1/0.648. */
  const float corner = ls_eval_vignette_factor(&corr, (float)(IMG_W - 1), (float)(IMG_H - 1));
  CHECK(fabsf(corner - 1.f / 0.648f) < 2e-3f,
        "vignetting: corner lift %.5f, expected %.5f", corner, 1.f / 0.648f);
  printf("knots: corner vignetting lift %.4f (%.2f EV)\n", corner, log2f(corner));
}

/* The projection grafted onto a table-resolved modifier has to land in that modifier's
 * OWN units, and the two resolvers do not share units: ls_modifier_init() normalizes radius
 * against the half SHORT SIDE, ls_modifier_init_knots() against the half DIAGONAL.
 *
 * The invariant that survives both is physical: the angle subtended at the frame corner is a
 * property of the lens and the sensor, not of anyone's normalization. Compute it each way
 * and the two must agree -- which is the whole claim ls_modifier_set_projection() rests on,
 * and the one that would silently mis-scale every reprojected pixel if the
 * aspect_ratio_correction factor were dropped or applied twice. */
static void test_projection_units(void)
{
  const float focal = 24.f, crop = 1.53f;
  const float ar = (float)IMG_W / (float)IMG_H;
  const float ar_corr = sqrtf(ar * ar + 1.f);

  /* A table-resolved modifier: aspect_ratio_correction is 1, radius 1 is the half diagonal. */
  ls_knots_t k;
  _make_identity(&k);
  ls_modifier_t mk;
  ls_modifier_init_knots(&mk, &k, IMG_W, IMG_H, 1.f, LS_ENABLE_DISTORTION, 0);
  CHECK(fabsf(mk.aspect_ratio_correction - 1.f) < 1e-6f,
        "projection: a table modifier should normalize against the diagonal, got arc=%.4f",
        mk.aspect_ratio_correction);

  ls_modifier_set_projection(&mk, LS_LENS_FISHEYE, LS_LENS_RECTILINEAR, focal, crop);
  CHECK((mk.enabled & LS_ENABLE_GEOMETRY) != 0, "projection: geometry did not enable");

  /* The corner sits at radius 1 in these coordinates, and at ar_corr in the other's. */
  const float theta_knots = atanf(1.f / mk.geom_focal);
  const float geom_focal_db = focal * crop * ar_corr / 21.633307f;
  const float theta_db = atanf(ar_corr / geom_focal_db);

  CHECK(fabsf(theta_knots - theta_db) < 1e-5f,
        "projection: corner half-angle differs between normalizations: %.6f vs %.6f rad",
        theta_knots, theta_db);
  printf("knots: corner half-angle %.5f rad, same either way\n", theta_knots);

  /* A pair that is not a function of radius alone is declined, not approximated. */
  ls_modifier_t mp;
  ls_modifier_init_knots(&mp, &k, IMG_W, IMG_H, 1.f, LS_ENABLE_DISTORTION, 0);
  CHECK(ls_modifier_set_projection(&mp, LS_LENS_RECTILINEAR, LS_LENS_PANORAMIC, focal, crop) == 0,
        "projection: panoramic should be declined");
  CHECK(mp.geometry_unsupported, "projection: panoramic should raise geometry_unsupported");
  CHECK((mp.enabled & LS_ENABLE_GEOMETRY) == 0,
        "projection: a declined change must leave the stage off");

  /* Same projection on both sides is not a change, and must not switch the stage on. */
  ls_modifier_t mi;
  ls_modifier_init_knots(&mi, &k, IMG_W, IMG_H, 1.f, LS_ENABLE_DISTORTION, 0);
  CHECK(ls_modifier_set_projection(&mi, LS_LENS_RECTILINEAR, LS_LENS_RECTILINEAR, focal, crop) == 0,
        "projection: an identity change should not enable geometry");
  CHECK(!mi.geometry_unsupported, "projection: an identity change is supported, not refused");
}

int main(void)
{
  test_identity();
  test_roundtrip();
  test_channels();
  test_refusal();
  test_chain();
  test_vignetting();
  test_projection_units();

  printf("knots: %d failure(s)\n", failures);
  return failures ? 1 : 0;
}
