/*
    LensSerious — the vendor decoders, held to ranges and to their own contract.

    Copyright (C) 2026 Aurélien PIERRE.  License: GPL-3.0-or-later (this test exercises
    the GPL vendor module — see lensserious_vendor.h).
*/

/** @file vendor.c
 *
 * @brief What can be asserted about the vendor decoders without a camera in hand.
 *
 * The decoders' ground truth is each maker's firmware, which nothing here can consult, so
 * exact numbers are not what is verified. What is:
 *
 *  1. **Valid data resolves, and stays finite.** Each vendor's plausible table yields the
 *     axes it carries, a positive autoscale, and knots that are finite, non-negative and
 *     bounded — the corruption class that would poison a lookup table downstream.
 *  2. **Empty data is absent, not an error.** A zero table means the body chose not to
 *     correct this lens; the resolve must say 0, never EBADDATA, and must leave the
 *     autoscale at exactly 1.
 *  3. **NULL finetune means "as measured".** The opposite convention — NULL rejected —
 *     shipped once and cost a consumer a silent no-op that took a pixel-diff to find, so
 *     it is pinned here: NULL and {1,1,1,1} must produce identical tables.
 *  4. **The finetune actually tunes.** Half strength must move the knots and shrink the
 *     autoscale; a slider that changes nothing is the bug the finetune tests exist for.
 *  5. **The OpcodeList3 parser survives hostile input.** Truncated headers, oversized
 *     param claims and unknown opcodes must skip or fail cleanly, never read past the
 *     buffer; and partial data yielded before a break is kept, because partial data is
 *     data.
 */

#include "lensserious_vendor.h"

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

/* Fujifilm pads vig_radius past the last knot by +1 per step, so with knots[8] = 0.9 the
 * upper entries reach ~6.9. Every other vendor keeps radii and factors near 1. This bound
 * is loose enough for the padding while still catching NaN/Inf/overflow. */
#define KNOT_HI 10.0f

static void check_knots_sane(const float *k, int n, const char *what)
{
  for(int i = 0; i < n; i++)
  {
    CHECK(!isnan(k[i]) && !isinf(k[i]), "%s[%d] is not finite", what, i);
    CHECK(k[i] >= 0.0f && k[i] < KNOT_HI, "%s[%d] = %f out of range", what, i, k[i]);
  }
}

static void check_table_sane(const ls_knots_t *k)
{
  check_knots_sane(k->radius, LS_MAX_KNOTS, "radius");
  check_knots_sane(k->vig_radius, LS_MAX_KNOTS, "vig_radius");
  check_knots_sane(k->cor_rgb[0], LS_MAX_KNOTS, "cor_rgb[0]");
  check_knots_sane(k->cor_rgb[1], LS_MAX_KNOTS, "cor_rgb[1]");
  check_knots_sane(k->cor_rgb[2], LS_MAX_KNOTS, "cor_rgb[2]");
  check_knots_sane(k->vig, LS_MAX_KNOTS, "vig");
}

/* ------------------------------------------------------------------------- */

static ls_vendor_data_t _sony_data(void)
{
  ls_vendor_data_t d;
  memset(&d, 0, sizeof(d));
  d.type = LS_VENDOR_SONY;
  d.u.sony.nc = 16;
  for(int i = 0; i < 16; i++)
  {
    d.u.sony.distortion[i] = (int16_t)(i + 1);
    d.u.sony.ca_r[i] = (int16_t)(i + 1);
    d.u.sony.ca_b[i] = (int16_t)(i + 1);
    d.u.sony.vignetting[i] = (int16_t)(i + 1);
  }
  return d;
}

static void test_sony(void)
{
  const ls_vendor_data_t d = _sony_data();

  const int axes = ls_vendor_axes(&d);
  CHECK(axes == (LS_ENABLE_DISTORTION | LS_ENABLE_TCA | LS_ENABLE_VIGNETTING),
        "sony: axes 0x%x, expected all three", axes);

  ls_knots_t k;
  float scale = 0.0f;
  const int got = ls_vendor_resolve(&d, NULL, IMG_W, IMG_H, &k, &scale);
  CHECK(got == axes, "sony: resolve 0x%x != axes 0x%x", got, axes);
  CHECK(k.n == 16, "sony: n = %d, expected 16", k.n);
  CHECK(k.vn == 16, "sony: vn = %d, expected 16", k.vn);
  CHECK(scale > 1.0f, "sony: autoscale %f, expected > 1 for a positive table", scale);
  check_table_sane(&k);
}

static void test_sony_empty(void)
{
  ls_vendor_data_t d;
  memset(&d, 0, sizeof(d));
  d.type = LS_VENDOR_SONY;
  d.u.sony.nc = 0;

  CHECK(ls_vendor_axes(&d) == 0, "sony empty: axes should be 0");

  ls_knots_t k;
  float scale = -1.0f;
  const int got = ls_vendor_resolve(&d, NULL, IMG_W, IMG_H, &k, &scale);
  CHECK(got == 0, "sony empty: resolve %d, expected 0 (absent, not an error)", got);
  CHECK(fabsf(scale - 1.0f) < 1e-6f, "sony empty: autoscale %f, expected exactly 1", scale);
  CHECK(k.n == 0 && k.vn == 0, "sony empty: table not zeroed");
}

static void test_sony_zero_tables(void)
{
  /* Valid nc, all coefficients zero: a body declining to correct this lens. Absent. */
  ls_vendor_data_t d;
  memset(&d, 0, sizeof(d));
  d.type = LS_VENDOR_SONY;
  d.u.sony.nc = 16;

  CHECK(ls_vendor_axes(&d) == 0, "sony zeros: axes should be 0");

  ls_knots_t k;
  float scale = -1.0f;
  const int got = ls_vendor_resolve(&d, NULL, IMG_W, IMG_H, &k, &scale);
  CHECK(got == 0, "sony zeros: resolve %d, expected 0 — a zero table is not malformed", got);
}

static void test_fuji(void)
{
  ls_vendor_data_t d;
  memset(&d, 0, sizeof(d));
  d.type = LS_VENDOR_FUJI;
  d.u.fuji.nc = 9;
  d.u.fuji.cropf = 1.0f;
  for(int i = 0; i < 9; i++)
  {
    d.u.fuji.knots[i] = 0.1f * (float)(i + 1);
    d.u.fuji.distortion[i] = 1.0f;
    d.u.fuji.ca_r[i] = 0.001f;
    d.u.fuji.ca_b[i] = 0.001f;
    d.u.fuji.vignetting[i] = 95.0f;
  }

  const int axes = ls_vendor_axes(&d);
  CHECK(axes == (LS_ENABLE_DISTORTION | LS_ENABLE_TCA | LS_ENABLE_VIGNETTING),
        "fuji: axes 0x%x, expected all three", axes);

  ls_knots_t k;
  float scale = 0.0f;
  const int got = ls_vendor_resolve(&d, NULL, IMG_W, IMG_H, &k, &scale);
  CHECK(got == axes, "fuji: resolve 0x%x != axes 0x%x", got, axes);
  CHECK(k.n == LS_MAX_KNOTS, "fuji: n = %d, expected %d (resampled)", k.n, LS_MAX_KNOTS);
  CHECK(scale > 0.0f, "fuji: autoscale %f, expected > 0", scale);
  check_table_sane(&k);
}

static void test_fuji_empty(void)
{
  ls_vendor_data_t d;
  memset(&d, 0, sizeof(d));
  d.type = LS_VENDOR_FUJI;

  ls_knots_t k;
  float scale = -1.0f;
  CHECK(ls_vendor_resolve(&d, NULL, IMG_W, IMG_H, &k, &scale) == 0, "fuji empty: expected 0");
  CHECK(fabsf(scale - 1.0f) < 1e-6f, "fuji empty: autoscale %f, expected exactly 1", scale);
}

static ls_vendor_data_t _dng_data(void)
{
  ls_vendor_data_t d;
  memset(&d, 0, sizeof(d));
  d.type = LS_VENDOR_DNG;
  d.u.dng.has_warp = 1;
  d.u.dng.warp_planes = 1;
  d.u.dng.warp_coeffs[0][0] = 1.002;
  d.u.dng.warp_coeffs[0][1] = -0.01;
  d.u.dng.warp_coeffs[0][2] = 0.02;
  d.u.dng.warp_coeffs[0][3] = -0.003;
  d.u.dng.has_vignette = 1;
  d.u.dng.vig_coeffs[0] = -0.3;
  d.u.dng.vig_coeffs[1] = 0.1;
  d.u.dng.vig_coeffs[2] = -0.05;
  d.u.dng.vig_coeffs[3] = 0.01;
  d.u.dng.vig_coeffs[4] = -0.001;
  return d;
}

static void test_dng(void)
{
  const ls_vendor_data_t d = _dng_data();

  /* One warp plane: distortion, no TCA. */
  const int axes = ls_vendor_axes(&d);
  CHECK(axes == (LS_ENABLE_DISTORTION | LS_ENABLE_VIGNETTING),
        "dng: axes 0x%x, expected distortion+vignetting, no TCA at 1 plane", axes);

  ls_knots_t k;
  float scale = 0.0f;
  const int got = ls_vendor_resolve(&d, NULL, IMG_W, IMG_H, &k, &scale);
  CHECK(got == axes, "dng: resolve 0x%x != axes 0x%x", got, axes);
  CHECK(k.n == LS_MAX_KNOTS, "dng: n = %d, expected %d (sampled)", k.n, LS_MAX_KNOTS);
  CHECK(scale > 0.0f, "dng: autoscale %f, expected > 0", scale);
  check_table_sane(&k);

  /* Three planes: the same warp now carries TCA too. */
  ls_vendor_data_t d3 = d;
  d3.u.dng.warp_planes = 3;
  d3.u.dng.warp_coeffs[1][0] = 1.001;
  d3.u.dng.warp_coeffs[2][0] = 1.003;
  CHECK(ls_vendor_axes(&d3) & LS_ENABLE_TCA, "dng: 3 planes must report TCA");
}

static void test_dng_empty(void)
{
  ls_vendor_data_t d;
  memset(&d, 0, sizeof(d));
  d.type = LS_VENDOR_DNG;

  ls_knots_t k;
  float scale = -1.0f;
  CHECK(ls_vendor_resolve(&d, NULL, IMG_W, IMG_H, &k, &scale) == 0, "dng empty: expected 0");
  CHECK(fabsf(scale - 1.0f) < 1e-6f, "dng empty: autoscale %f, expected exactly 1", scale);
}

static void test_olympus(void)
{
  ls_vendor_data_t d;
  memset(&d, 0, sizeof(d));
  d.type = LS_VENDOR_OLYMPUS;
  d.u.olympus.has_dist = 1;
  d.u.olympus.dist[0] = 0.01f;
  d.u.olympus.dist[1] = -0.02f;
  d.u.olympus.dist[2] = 0.005f;
  d.u.olympus.dist[3] = 1.0f;
  d.u.olympus.has_ca = 1;
  d.u.olympus.ca[0] = 0.002f;
  d.u.olympus.ca[1] = -0.001f;
  d.u.olympus.ca[2] = 0.0005f;
  d.u.olympus.ca[3] = -0.001f;
  d.u.olympus.ca[4] = 0.0005f;
  d.u.olympus.ca[5] = -0.0001f;

  /* Olympus publishes no vignetting; the axis must be absent, not identity-present. */
  const int axes = ls_vendor_axes(&d);
  CHECK(axes == (LS_ENABLE_DISTORTION | LS_ENABLE_TCA),
        "olympus: axes 0x%x, expected distortion+TCA and NO vignetting", axes);

  ls_knots_t k;
  float scale = 0.0f;
  const int got = ls_vendor_resolve(&d, NULL, IMG_W, IMG_H, &k, &scale);
  CHECK(got == axes, "olympus: resolve 0x%x != axes 0x%x", got, axes);
  CHECK(k.vn == 0, "olympus: vn = %d, expected 0 — no vignetting to serve", k.vn);
  CHECK(scale > 0.0f, "olympus: autoscale %f, expected > 0", scale);
  check_table_sane(&k);
}

static void test_olympus_empty(void)
{
  ls_vendor_data_t d;
  memset(&d, 0, sizeof(d));
  d.type = LS_VENDOR_OLYMPUS;

  ls_knots_t k;
  float scale = -1.0f;
  CHECK(ls_vendor_resolve(&d, NULL, IMG_W, IMG_H, &k, &scale) == 0, "olympus empty: expected 0");
  CHECK(fabsf(scale - 1.0f) < 1e-6f, "olympus empty: autoscale %f, expected 1", scale);
}

static void test_none(void)
{
  ls_vendor_data_t d;
  memset(&d, 0, sizeof(d));

  ls_knots_t k;
  float scale = -1.0f;
  CHECK(ls_vendor_axes(&d) == 0, "none: axes should be 0");
  CHECK(ls_vendor_resolve(&d, NULL, IMG_W, IMG_H, &k, &scale) == 0, "none: expected 0");
  CHECK(fabsf(scale - 1.0f) < 1e-6f, "none: autoscale %f, expected exactly 1", scale);
  CHECK(ls_vendor_resolve(NULL, NULL, IMG_W, IMG_H, &k, &scale) == 0, "NULL data: expected 0");
}

static void test_finetune(void)
{
  const ls_vendor_data_t d = _sony_data();

  ls_knots_t k_null, k_one, k_half;
  float s_null = 0.f, s_one = 0.f, s_half = 0.f;
  const ls_vendor_finetune_t one = { 1.f, 1.f, 1.f, 1.f };
  const ls_vendor_finetune_t half = { 0.5f, 0.5f, 0.5f, 0.5f };

  ls_vendor_resolve(&d, NULL, IMG_W, IMG_H, &k_null, &s_null);
  ls_vendor_resolve(&d, &one, IMG_W, IMG_H, &k_one, &s_one);
  ls_vendor_resolve(&d, &half, IMG_W, IMG_H, &k_half, &s_half);

  /* NULL is {1,1,1,1}, bit for bit. */
  CHECK(memcmp(&k_null, &k_one, sizeof(k_null)) == 0, "NULL finetune differs from all-1.0");
  CHECK(s_null == s_one, "NULL finetune autoscale %f != all-1.0 autoscale %f", s_null, s_one);

  /* Half strength shrinks the correction, and therefore the zoom that clears it. */
  CHECK(s_half > 0.0f && s_half < s_one, "half finetune: autoscale %f, expected in (0, %f)",
        s_half, s_one);
  CHECK(memcmp(&k_half, &k_one, sizeof(k_half)) != 0, "half finetune changed nothing");
}

/* ------------------------------------------------------------------------- */
/* OpcodeList3                                                                */
/* ------------------------------------------------------------------------- */

static void _put_be_long(uint8_t *p, uint32_t v)
{
  p[0] = (uint8_t)(v >> 24); p[1] = (uint8_t)(v >> 16); p[2] = (uint8_t)(v >> 8); p[3] = (uint8_t)v;
}

static void _put_be_double(uint8_t *p, double d)
{
  uint64_t v;
  memcpy(&v, &d, sizeof(v));
  for(int i = 0; i < 8; i++) p[i] = (uint8_t)(v >> (56 - 8 * i));
}

/* One opcode entry: id, DNG version, flags, param_size, then the params. Returns bytes
 * written. */
static size_t _put_opcode(uint8_t *p, uint32_t id, const uint8_t *param, uint32_t param_size)
{
  _put_be_long(p, id);
  _put_be_long(p + 4, 0x01030000); /* DNG 1.3, as writers stamp it */
  _put_be_long(p + 8, 1);          /* optional */
  _put_be_long(p + 12, param_size);
  memcpy(p + 16, param, param_size);
  return 16 + param_size;
}

static size_t _make_warp_param(uint8_t *p, uint32_t planes)
{
  _put_be_long(p, planes);
  size_t off = 4;
  for(uint32_t pl = 0; pl < planes; pl++)
    for(int c = 0; c < 6; c++)
    {
      _put_be_double(p + off, c == 0 ? 1.002 : (c == 1 ? -0.01 : 0.0));
      off += 8;
    }
  _put_be_double(p + off, 0.5); off += 8;
  _put_be_double(p + off, 0.5); off += 8;
  return off;
}

static void test_opcodelist(void)
{
  uint8_t param[4 + 3 * 48 + 16];
  uint8_t blob[512];

  /* A list with one 1-plane warp and one vignette. */
  const size_t warp_len = _make_warp_param(param, 1);
  _put_be_long(blob, 2);
  size_t off = 4 + _put_opcode(blob + 4, 1, param, (uint32_t)warp_len);

  uint8_t vparam[7 * 8];
  for(int c = 0; c < 5; c++) _put_be_double(vparam + c * 8, -0.1 * (c + 1));
  _put_be_double(vparam + 40, 0.5);
  _put_be_double(vparam + 48, 0.5);
  off += _put_opcode(blob + off, 3, vparam, sizeof(vparam));

  ls_vendor_dng_t out;
  int got = ls_vendor_parse_dng_opcodelist3(blob, off, &out);
  CHECK(got == (LS_ENABLE_DISTORTION | LS_ENABLE_VIGNETTING),
        "opcodes: 0x%x, expected warp+vignette", got);
  CHECK(out.warp_planes == 1, "opcodes: %u planes, expected 1", out.warp_planes);
  CHECK(fabs(out.warp_coeffs[0][0] - 1.002) < 1e-12, "opcodes: kr0 %f", out.warp_coeffs[0][0]);
  CHECK(fabs(out.vig_coeffs[0] - (-0.1)) < 1e-12, "opcodes: k0 %f", out.vig_coeffs[0]);

  /* An unknown opcode id between the two is skipped, not fatal. */
  uint8_t blob2[600];
  _put_be_long(blob2, 3);
  size_t off2 = 4 + _put_opcode(blob2 + 4, 1, param, (uint32_t)warp_len);
  const uint8_t junk[8] = { 0 };
  off2 += _put_opcode(blob2 + off2, 9 /* GainMap */, junk, sizeof(junk));
  off2 += _put_opcode(blob2 + off2, 3, vparam, sizeof(vparam));
  got = ls_vendor_parse_dng_opcodelist3(blob2, off2, &out);
  CHECK(got == (LS_ENABLE_DISTORTION | LS_ENABLE_VIGNETTING),
        "opcodes+junk: 0x%x, expected warp+vignette", got);

  /* Truncated before anything usable: an error, not silence. */
  got = ls_vendor_parse_dng_opcodelist3(blob, 3, &out);
  CHECK(got == LS_VENDOR_EBADDATA, "3-byte blob: %d, expected EBADDATA", got);
  _put_be_long(blob2, 2);
  got = ls_vendor_parse_dng_opcodelist3(blob2, 10, &out);
  CHECK(got == LS_VENDOR_EBADDATA, "truncated header: %d, expected EBADDATA", got);

  /* Truncated AFTER a full opcode: the yielded data is kept. */
  _put_be_long(blob, 2);
  got = ls_vendor_parse_dng_opcodelist3(blob, 4 + 16 + warp_len + 5, &out);
  CHECK(got == LS_ENABLE_DISTORTION, "partial list: 0x%x, expected the warp kept", got);

  /* A param_size claiming past the buffer must not be read. */
  uint8_t blob3[32];
  _put_be_long(blob3, 1);
  _put_be_long(blob3 + 4, 1);
  _put_be_long(blob3 + 8, 0x01030000);
  _put_be_long(blob3 + 12, 1);
  _put_be_long(blob3 + 16, 0xFFFFFFFFu); /* param_size: 4 GB */
  got = ls_vendor_parse_dng_opcodelist3(blob3, sizeof(blob3), &out);
  CHECK(got == LS_VENDOR_EBADDATA, "hostile param_size: %d, expected EBADDATA", got);

  /* An empty list carries nothing, and that is not an error. */
  _put_be_long(blob, 0);
  got = ls_vendor_parse_dng_opcodelist3(blob, 4, &out);
  CHECK(got == 0, "empty list: %d, expected 0", got);
}

int main(void)
{
  test_sony();
  test_sony_empty();
  test_sony_zero_tables();
  test_fuji();
  test_fuji_empty();
  test_dng();
  test_dng_empty();
  test_olympus();
  test_olympus_empty();
  test_none();
  test_finetune();
  test_opcodelist();

  if(failures) { printf("%d FAILURE(S)\n", failures); return 1; }
  printf("vendor: all checks passed\n");
  return 0;
}
