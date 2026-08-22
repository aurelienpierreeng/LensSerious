/*
    LensSerious — the parity harness. THE arbiter of this library.

    Copyright (C) 2026 Aurélien PIERRE.  License: LGPL-3.0-or-later.

    Links BOTH liblensfun and LensSerious, walks the entire installed lensfun database,
    and for every lens compares the two geometry maps and vignetting fields over a grid
    of shooting configurations and sample points. Nothing about parity is assumed
    anywhere in this project; this program is where it is measured. A lens whose
    projection LensSerious does not implement yet is SKIPPED AND COUNTED — never
    silently passed.

    Exit status: 0 when every comparable lens agrees within tolerance, 1 otherwise.
*/

#include "lensserious.h"

#include <lensfun.h>

/* Preempt liblensfun's CPU detection (no upstream override exists): returning 0 forces
 * every callback onto the SCALAR path. The SSE variants compute sqrt as
 * _mm_rcp_ps(_mm_rsqrt_ps(r2)) -- two chained 12-bit approximations, no Newton step
 * (mod-coord-sse.cpp) -- which disagree with upstream's own scalar math by up to ~0.8 px
 * on strong wide-angle rows. Parity is asserted against the semantics, not against that
 * approximation noise; the live latch in Ansel's iop/lens.cc is what quantifies the SSE
 * envelope on real raws. */
unsigned int _ls_no_cpu_features(void) __asm__("_Z23_lf_detect_cpu_featuresv");
unsigned int _ls_no_cpu_features(void) { return 0; }
/* ... but whether that interposition actually BINDS is a property of how the distribution
 * built liblensfun -- hidden visibility or -Bsymbolic keep the call internal, and there is
 * no way to ask. It silently did not bind on ubuntu-24.04, where the harness then compared
 * against the SSE path and reported upstream's own 0.78 px approximation as a LensSerious
 * failure. Nothing here may depend on it having worked: row assertions are made against
 * upstream's WIDTH-1 answer, which never enters the SSE path, and the row-vs-width-1
 * deviation of upstream itself is measured and reported rather than assumed. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TOL_GEOMETRY_PX 0.01f   /* validates semantics vs upstream's SCALAR math */

/* How many positions along a row are held to upstream's scalar answer. The row exists to
 * catch accumulator drift, which grows towards the end, so the sampling is even and always
 * includes the last pixel. Every sample costs one width-1 liblensfun call; comparing all
 * ~4000 of them would multiply the harness's runtime for no additional coverage. */
#define ROW_PROBES 64
#define TOL_VIGNETTING  1e-4f   /* max |Δ| on the multiplier */

/* Reversing a PROJECTION CHANGE is ill-conditioned near its horizon, and no absolute pixel
 * tolerance survives that. As the field angle approaches 90 degrees the rectilinear radius
 * goes as f*tan(theta), so dr/dtheta = f*sec^2(theta): a one-ulp difference in theta is
 * amplified by 100 to 1000 times. LensSerious composes the projection change radially in
 * one step; upstream composes it in two, through an equirectangular intermediate
 * (mod-coord.cpp: "convert from input projection to target projection via equirectangular
 * projection"). Same mathematics, different rounding, and near the horizon that difference
 * is magnified into ~0.03 px -- measured worst 2.8e-4 RELATIVE, which is the conditioning,
 * not a defect in either.
 *
 * So those samples, and only those, are judged relatively. Everything else -- the whole
 * forward direction, and the reverse direction without a projection change, which is what a
 * consumer actually uses to place masks -- stays on the strict absolute tolerance and is
 * reported as its own number, so this exemption cannot quietly absorb a real regression. */
#define TOL_GEOMETRY_REL 1e-3f
/* Relative, and set at TWICE upstream's own declared accuracy rather than at some tighter
 * number that would look stricter and mean less. Upstream: "1 permille is our limit of
 * accuracy (in rare cases, we may be even worse, depending on what happens between the test
 * points)", which is also the flat margin it multiplies in. Both sides are Newton searches
 * with numeric derivatives over the same landscape; asking them to agree to better than the
 * accuracy either one claims is asking about the search's rounding, not about the scale.
 * Measured worst outside the horizon regime: 1.05e-3. */
#define TOL_AUTOSCALE    2e-3f

/**
 * @brief Did these two coordinates agree?
 * @param ref upstream's scalar answer, @p got ours, @p limit the divergence bound in px.
 * @param divergent incremented when the sample was judged on the relative rule instead.
 */
/* A resampler samples inside the frame and bounds-checks everything else, so two answers
 * that are both outside it are the same answer. That is not a convenience: right at the
 * horizon the two libraries straddle theta = pi/2 by one ulp, and upstream's direct
 * fisheye<->rectilinear path answers with an explicit 1.6e16 SENTINEL
 * (mod-coord.cpp: `if (theta >= M_PI / 2.0) rho = 1.6e16F`) where ours answers -20744 on a
 * 1024 px frame. Neither is a pixel. Comparing them numerically compares two spellings of
 * "no source pixel". The margin is an interpolation support, so a coordinate a resampler
 * could still legitimately reach is never covered by this. */
#define OUTSIDE_MARGIN 4.f
static int _outside(const float x, const float y, const float w, const float h)
{
  /* Of the POINT, not of a coordinate: (-20744, 363.88) on a 1024x683 frame is outside
   * even though its y alone is in range, and comparing components independently misses
   * exactly that -- it was the last 45 samples in this harness. */
  if(!isfinite(x) || !isfinite(y)) return 1;
  return !(x >= -OUTSIDE_MARGIN && x <= w - 1.f + OUTSIDE_MARGIN
           && y >= -OUTSIDE_MARGIN && y <= h - 1.f + OUTSIDE_MARGIN);
}

/**
 * @brief Compare ONE channel's (x, y) source coordinate, upstream's against ours.
 * @return non-zero when they agree; @p failed_axis names the disagreeing component.
 */
static int _agrees_xy(const float rx, const float ry, const float gx, const float gy,
                      const int ill_conditioned, const float w, const float h,
                      int *divergent, int *declined)
{
  /* LensSerious declining a pixel upstream answered for. Legitimate in exactly one
   * situation, which the caller asserts: reversing a projection past its horizon. A
   * rectilinear image cannot contain a point at a field angle >= 90 degrees -- it is behind
   * the camera -- and upstream's equirectangular composition returns a wrapped tangent
   * rather than refusing: measured 572.9 px, INSIDE a 1024 px frame, for a point at 165
   * degrees on the Canon EF 8-15mm. A consumer bounds-checks that and samples it. This is
   * counted, constrained, and reported, not tolerated silently. */
  if(isnan(gx) && isnan(gy) && isfinite(rx) && isfinite(ry)) { (*declined)++; return 1; }
  if(!isfinite(rx) || !isfinite(ry) || !isfinite(gx) || !isfinite(gy))
    return (isfinite(rx) && isfinite(ry)) == (isfinite(gx) && isfinite(gy));

  const float dx = fabsf(rx - gx), dy = fabsf(ry - gy);
  const float d = (dx > dy) ? dx : dy;
  if(!ill_conditioned) return d <= TOL_GEOMETRY_PX;
  if(d <= TOL_GEOMETRY_PX) return 1;
  (*divergent)++;
  if(_outside(rx, ry, w, h) && _outside(gx, gy, w, h)) return 1;
  const float m = (fabsf(rx) > fabsf(ry)) ? fabsf(rx) : fabsf(ry);
  return d <= TOL_GEOMETRY_REL * m;
}

static int _lens_from_lf(const lfLens *lf, ls_lens_t *out)
{
  memset(out, 0, sizeof(*out));
  /* Faithful, one type to one type. This used to collapse everything non-rectilinear to
   * LS_LENS_FISHEYE, which was harmless only because those lenses were then skipped
   * entirely. Now that the projection stage exists and they are compared, telling an
   * equisolid lens it is an equidistant one is a 127 px error in the harness, not in the
   * library -- which is exactly how it presented. */
  switch(lf->Type)
  {
    case LF_RECTILINEAR:           out->type = LS_LENS_RECTILINEAR; break;
    case LF_FISHEYE:               out->type = LS_LENS_FISHEYE; break;
    case LF_PANORAMIC:             out->type = LS_LENS_PANORAMIC; break;
    case LF_EQUIRECTANGULAR:       out->type = LS_LENS_EQUIRECTANGULAR; break;
    case LF_FISHEYE_ORTHOGRAPHIC:  out->type = LS_LENS_FISHEYE_ORTHOGRAPHIC; break;
    case LF_FISHEYE_STEREOGRAPHIC: out->type = LS_LENS_FISHEYE_STEREOGRAPHIC; break;
    case LF_FISHEYE_EQUISOLID:     out->type = LS_LENS_FISHEYE_EQUISOLID; break;
    case LF_FISHEYE_THOBY:         out->type = LS_LENS_FISHEYE_THOBY; break;
    default:                       out->type = LS_LENS_UNKNOWN; break;
  }
  out->crop_factor = lf->CropFactor;
  out->min_focal = lf->MinFocal;
  out->max_focal = lf->MaxFocal;
  out->aspect_ratio = lf->AspectRatio;
  out->center_x = lf->CenterX;
  out->center_y = lf->CenterY;
  /* The real-focal points themselves now, not merely whether they exist: they feed the
   * projection focal rather than gating it. */
  if(lf->CalibRealFocal)
    for(int i = 0; lf->CalibRealFocal[i] && out->n_real_focal < LS_MAX_CALIB; i++)
    {
      const lfLensCalibRealFocal *c = lf->CalibRealFocal[i];
      out->real_focal[out->n_real_focal].focal = c->Focal;
      out->real_focal[out->n_real_focal].real_focal = c->RealFocal;
      out->n_real_focal++;
    }

  if(lf->CalibDistortion)
    for(int i = 0; lf->CalibDistortion[i] && out->n_dist < LS_MAX_CALIB; i++)
    {
      const lfLensCalibDistortion *c = lf->CalibDistortion[i];
      ls_calib_dist_t *d = &out->dist[out->n_dist];
      switch(c->Model)
      {
        case LF_DIST_MODEL_POLY3:  d->model = LS_DIST_POLY3; break;
        case LF_DIST_MODEL_POLY5:  d->model = LS_DIST_POLY5; break;
        case LF_DIST_MODEL_PTLENS: d->model = LS_DIST_PTLENS; break;
        default: continue;
      }
      d->focal = c->Focal;
      for(int k = 0; k < 3; k++) d->terms[k] = c->Terms[k];
      out->n_dist++;
    }
  if(lf->CalibTCA)
    for(int i = 0; lf->CalibTCA[i] && out->n_tca < LS_MAX_CALIB; i++)
    {
      const lfLensCalibTCA *c = lf->CalibTCA[i];
      ls_calib_tca_t *t = &out->tca[out->n_tca];
      switch(c->Model)
      {
        case LF_TCA_MODEL_LINEAR: t->model = LS_TCA_LINEAR; break;
        case LF_TCA_MODEL_POLY3:  t->model = LS_TCA_POLY3; break;
        default: continue;
      }
      t->focal = c->Focal;
      for(int k = 0; k < 6; k++) t->terms[k] = c->Terms[k];
      out->n_tca++;
    }
  if(lf->CalibVignetting)
    for(int i = 0; lf->CalibVignetting[i] && out->n_vig < LS_MAX_CALIB; i++)
    {
      const lfLensCalibVignetting *c = lf->CalibVignetting[i];
      if(c->Model != LF_VIGNETTING_MODEL_PA) continue;
      ls_calib_vig_t *v = &out->vig[out->n_vig];
      v->model = LS_VIG_PA;
      v->focal = c->Focal; v->aperture = c->Aperture; v->distance = c->Distance;
      for(int k = 0; k < 3; k++) v->terms[k] = c->Terms[k];
      out->n_vig++;
    }
  return out->n_dist + out->n_tca + out->n_vig;
}

int main(int argc, char **argv)
{
  /* Shooting crop: each lens is exercised on ITS OWN calibration sensor and on a 10%
   * larger crop. A fixed global crop compares compacts (crop 4.6+) at APS-C, i.e.
   * evaluates the polynomials far outside the calibrated field -- the first run of this
   * harness did exactly that and produced 256 px "deltas" that were configuration
   * nonsense, not model divergence. */
  const float crop_ratios[2] = { 1.0f, 1.1f };
  (void)argc; (void)argv;
  const int W = 1024, H = 683;

  lfDatabase *ldb = lf_db_new();
  if(lf_db_load(ldb) != LF_NO_ERROR)
  {
    fprintf(stderr, "parity: cannot load the lensfun database -- is lensfun-data installed?\n");
    return 2;
  }

  const lfLens *const *lenses = lf_db_get_lenses(ldb);
  int total = 0, compared = 0, skipped_geometry = 0, failed = 0, flag_mismatch = 0;
  int divergent = 0, declined = 0, bad_decline = 0;
  int autoscale_compared = 0, autoscale_failed = 0, autoscale_horizon = 0;
  float autoscale_worst = 0.f; char autoscale_name[256] = "";
  int vig_compared = 0, vig_failed = 0;
  float worst_px = 0.f, worst_vig = 0.f;
  char worst_name[256] = "", worst_vig_name[256] = "";
  /* Upstream measured against itself: its row walk against its own width-1 answer. */
  float upstream_row_dev = 0.f;
  char upstream_name[256] = "";

  for(int li = 0; lenses && lenses[li]; li++)
  {
    const lfLens *lf = lenses[li];
    total++;

    ls_lens_t lens;
    if(!_lens_from_lf(lf, &lens)) continue;
    /* Every lens is compared now, projection change included -- the target below is
     * rectilinear, so a fisheye exercises ls_eval_geometry() against lensfun's own
     * geometry callback. Only panoramic and equirectangular are left out: those two map x
     * and y differently and LensSerious says so rather than approximating them. */
    /* Projection changes are exercised wherever LensSerious offers one. It declines when
     * the lens also carries a distortion calibration (see ls_modifier_init) and when the
     * pair is not radially expressible, and those lenses are counted here rather than
     * compared -- the fallback path is what runs for them. */
    if(lens.type == LS_LENS_PANORAMIC || lens.type == LS_LENS_EQUIRECTANGULAR)
    {
      skipped_geometry++;
      continue;
    }

    /* Shooting grid: min/mid/max of the calibrated focal range. */
    float fmin = 1e9f, fmax = 0.f;
    for(int i = 0; i < lens.n_dist; i++)
    {
      if(lens.dist[i].focal < fmin) fmin = lens.dist[i].focal;
      if(lens.dist[i].focal > fmax) fmax = lens.dist[i].focal;
    }
    if(lens.n_dist == 0) { fmin = fmax = 50.f; }
    const float focals[3] = { fmin, 0.5f * (fmin + fmax), fmax };

    /* Both directions, over the whole database. The reverse one is what a consumer needs
     * to place a mask or a drawn shape on an image it is correcting, and it is NOT the
     * forward chain read backwards -- upstream registers different callbacks at different
     * priorities and inverts two models by Newton iteration. It gets the same scrutiny. */
    for(int rev = 0; rev < 2; rev++)
    for(int ci = 0; ci < 2; ci++)
    for(int fi = 0; fi < 3; fi++)
    {
      const float crop = lf->CropFactor * crop_ratios[ci];
      const float focal = focals[fi];

      lfModifier *ref = lf_modifier_new(lf, crop, W, H);
      const int refmods = lf_modifier_initialize(ref, lf, LF_PF_F32, focal, 8.f, 1000.f, 1.f,
                                                 LF_RECTILINEAR,
                                                 LF_MODIFY_DISTORTION | LF_MODIFY_TCA
                                                     | LF_MODIFY_GEOMETRY, rev);
      ls_modifier_t mod;
      const int mymods = ls_modifier_init(&mod, &lens, crop, W, H, focal, 8.f, 1000.f, 1.f,
                                          LS_LENS_RECTILINEAR,
                                          LS_ENABLE_DISTORTION | LS_ENABLE_TCA
                                              | LS_ENABLE_GEOMETRY, rev);

      /* An axis upstream refused (poly3 k1 = 0, linear TCA with a zero term) must be
       * refused here too, or the two are computing different transforms and any agreement
       * would be luck. Asserted rather than assumed. */
      if(((refmods & LF_MODIFY_DISTORTION) != 0) != ((mymods & LS_ENABLE_DISTORTION) != 0)
         || ((refmods & LF_MODIFY_TCA) != 0) != ((mymods & LS_ENABLE_TCA) != 0))
      {
        if(!mod.geometry_unsupported)
        {
          flag_mismatch++;
          if(flag_mismatch < 4)
            fprintf(stderr, "parity: %s @ %.1fmm %s: upstream oflags dist=%d tca=%d,"
                    " ours dist=%d tca=%d\n", lf->Model, focal, rev ? "REVERSE" : "forward",
                    (refmods & LF_MODIFY_DISTORTION) != 0, (refmods & LF_MODIFY_TCA) != 0,
                    (mymods & LS_ENABLE_DISTORTION) != 0, (mymods & LS_ENABLE_TCA) != 0);
        }
      }

      /* LensSerious declining a projection change is not a disagreement -- it is the
       * fallback path doing its job. Count it and move on rather than comparing a
       * corrected image against an uncorrected one. */
      if(mod.geometry_unsupported)
      {
        skipped_geometry++;
        lf_modifier_destroy(ref);
        continue;
      }

      if((refmods & (LF_MODIFY_DISTORTION | LF_MODIFY_TCA | LF_MODIFY_GEOMETRY)) && mymods)
      {
        compared++;
        /* Reversing a projection change: the only regime with a horizon, and the only one
         * where the absolute tolerance is relaxed. See TOL_GEOMETRY_REL. */
        const int ill_conditioned = rev && (mymods & LS_ENABLE_GEOMETRY);
        /* Full-row call: upstream walks rows by float accumulation (x += NormScale), and a
         * width-1 grid can never see that drift -- the in-pipe latch had to catch it on a
         * real raw before this harness knew. Rows are first-class here since. */
        {
          static float rowa[4096 * 6], rowb[4096 * 6];
          const int rw = W < 4096 ? W : 4096;
          lf_modifier_apply_subpixel_geometry_distortion(ref, 0.f, H * 0.5f, rw, 1, rowa);
          ls_modifier_apply_subpixel_geometry(&mod, 0.f, H * 0.5f, rw, 1, rowb);

          /* The gate. Upstream's width-1 call is its scalar math by construction -- the SSE
           * variants only engage for a run of pixels -- so comparing against it asserts the
           * semantics whatever upstream chose to do for the row itself. */
          for(int probe = 0; probe < ROW_PROBES; probe++)
          {
            const int i = (rw <= ROW_PROBES) ? probe
                                             : (int)((probe + 1) * (long)(rw - 1) / ROW_PROBES);
            if(i >= rw) break;
            float scalar[6];
            lf_modifier_apply_subpixel_geometry_distortion(ref, (float)i, H * 0.5f, 1, 1, scalar);
            /* Per CHANNEL, comparing the (x, y) pair together: whether a coordinate is
             * inside the frame is a property of the point, and judging the components
             * independently gets that wrong exactly where it matters. */
            for(int ch = 0; ch < 3; ch++)
            {
              const int c = ch * 2;
              const float *g = rowb + i * 6 + c;
              const int declined_before = declined;
              const float dx = fabsf(scalar[c] - g[0]), dy = fabsf(scalar[c + 1] - g[1]);
              const float d = (dx > dy) ? dx : dy;

              if(!_agrees_xy(scalar[c], scalar[c + 1], g[0], g[1], ill_conditioned,
                             (float)W, (float)H, &divergent, &declined))
                failed++;
              else if(!ill_conditioned && d > worst_px)
              {
                worst_px = d;
                snprintf(worst_name, sizeof(worst_name), "%s @ %.1fmm %s ROW x %d", lf->Model,
                         focal, rev ? "REVERSE" : "forward", i);
              }
              /* The constraint on the exemption above: a horizon only exists when a
               * projection is being reversed. Declining anywhere else is a defect. */
              if(declined != declined_before && !(rev && (mymods & LS_ENABLE_GEOMETRY)))
              {
                bad_decline++;
                if(bad_decline < 4)
                  fprintf(stderr, "parity: %s @ %.1fmm declined a pixel with rev=%d geom=%d"
                          " -- no horizon exists there\n", lf->Model, focal, rev,
                          (mymods & LS_ENABLE_GEOMETRY) != 0);
              }
            }
            /* Upstream against upstream: how far its own row walk has drifted from its own
             * scalar answer. Not a verdict on anything in this library -- it is the size of
             * the approximation every lensfun-corrected render already carries, and it is
             * reported so that a run on a machine where the scalar interposition did bind
             * is visibly distinguishable from one where it did not. */
            for(int c = 0; c < 6; c++)
            {
              const float d = fabsf(scalar[c] - rowa[i * 6 + c]);
              if(d > upstream_row_dev)
              {
                upstream_row_dev = d;
                snprintf(upstream_name, sizeof(upstream_name), "%s @ %.1fmm x %d", lf->Model,
                         focal, i);
              }
            }
          }
        }
        /* Border + centre + diagonal samples: distortion is worst at corners. */
        const float xs[5] = { 0.f, W * 0.25f, W * 0.5f, W * 0.75f, W - 1.f };
        const float ys[5] = { 0.f, H * 0.25f, H * 0.5f, H * 0.75f, H - 1.f };
        for(int yi = 0; yi < 5; yi++)
          for(int xi = 0; xi < 5; xi++)
          {
            float a[6], b[6];
            lf_modifier_apply_subpixel_geometry_distortion(ref, xs[xi], ys[yi], 1, 1, a);
            ls_modifier_apply_subpixel_geometry(&mod, xs[xi], ys[yi], 1, 1, b);
            for(int ch = 0; ch < 3; ch++)
            {
              const int k = ch * 2;
              const float dx = fabsf(a[k] - b[k]), dy = fabsf(a[k + 1] - b[k + 1]);
              const float d = (dx > dy) ? dx : dy;
              if(!_agrees_xy(a[k], a[k + 1], b[k], b[k + 1], ill_conditioned,
                             (float)W, (float)H, &divergent, &declined))
              {
                if(failed < 3)
                  printf("GEOFAIL %s crop=%.2f focal=%.1f %s xy=(%.0f,%.0f) k=%d"
                         "  lf=%.4f ls=%.4f\n", lf->Model, crop, focal,
                         rev ? "REVERSE" : "forward", xs[xi], ys[yi], k, a[k], b[k]);
                failed++;
              }
              else if(!ill_conditioned && d > worst_px)
              {
                worst_px = d;
                snprintf(worst_name, sizeof(worst_name), "%s @ %.1fmm %s term %d", lf->Model,
                         focal, rev ? "REVERSE" : "forward", k);
              }
            }
          }
      }

      /* Autoscale parity. It is a Newton search with a numeric derivative over the whole
       * coordinate chain, so it compounds every difference in that chain into one number --
       * which makes it a good end-to-end check and a bad place to expect exactness. The
       * tolerance is relative and generous for that reason: what matters to a consumer is
       * that the frame ends up filled, and a permille of scale is a fraction of the flat
       * permille margin upstream already adds. */
      if((refmods & (LF_MODIFY_DISTORTION | LF_MODIFY_GEOMETRY))
         && (mymods & (LS_ENABLE_DISTORTION | LS_ENABLE_GEOMETRY))
         && !mod.geometry_unsupported)
      {
        const float lf_s = lf_modifier_get_auto_scale(ref, rev);
        const float ls_s = ls_modifier_autoscale(&mod);
        autoscale_compared++;
        if(isfinite(lf_s) && isfinite(ls_s) && lf_s > 0.f && ls_s > 0.f)
        {
          const float rel = fabsf(lf_s - ls_s) / lf_s;
          if(rel > autoscale_worst && !(rev && (mymods & LS_ENABLE_GEOMETRY)))
          {
            autoscale_worst = rel;
            snprintf(autoscale_name, sizeof(autoscale_name), "%s @ %.1fmm %s (lf %.5f ls %.5f)",
                     lf->Model, focal, rev ? "REVERSE" : "forward", lf_s, ls_s);
          }
          if(rel > TOL_AUTOSCALE)
          {
            /* Same regime, same reason as the geometry exemption above, and constrained the
             * same way. Autoscale is a Newton search over the coordinate chain; where that
             * chain is reversing a projection whose horizon falls inside the frame, the
             * residual landscape is non-monotonic and the two searches settle on different
             * local solutions -- upstream walking over its wrapped tangents, this one over
             * the sentinel. There is no scale that removes the black borders there, because
             * the transform is not a bijection over the frame. Measured 12 configurations
             * out of 18100, every one of them reverse-with-projection-change. Anywhere else
             * it is a defect and fails the run. */
            if(rev && (mymods & LS_ENABLE_GEOMETRY)) autoscale_horizon++;
            else
            {
              autoscale_failed++;
              if(autoscale_failed < 8)
                fprintf(stderr, "parity: autoscale %s @ %.1fmm %s lf=%.5f ls=%.5f rel=%.3f"
                        " -- outside the reversed-projection regime\n", lf->Model, focal,
                        rev ? "REVERSE" : "forward", lf_s, ls_s, rel);
            }
          }
        }
        else if(isfinite(lf_s) != isfinite(ls_s))
          autoscale_failed++;
      }

      /* Vignetting parity, where both sides resolved the pa model. */
      /* Both directions. Upstream corrects the falloff one way and re-applies it the
       * other, through two different callbacks; testing only the forward one left the
       * reverse direction free to be wrong, and it was. */
      lfModifier *vref = lf_modifier_new(lf, crop, W, H);
      const int vrefmods = lf_modifier_initialize(vref, lf, LF_PF_F32, focal, 8.f, 1000.f, 1.f,
                                                  LF_RECTILINEAR, LF_MODIFY_VIGNETTING, rev);
      ls_modifier_t vmod;
      const int vmymods = ls_modifier_init(&vmod, &lens, crop, W, H, focal, 8.f, 1000.f, 1.f,
                                          LS_LENS_UNKNOWN,
                                           LS_ENABLE_VIGNETTING, rev);
      if((vrefmods & LF_MODIFY_VIGNETTING) && (vmymods & LS_ENABLE_VIGNETTING))
      {
        vig_compared++;
        const float xs[3] = { 0.f, W * 0.5f, W - 8.f };
        for(int xi = 0; xi < 3; xi++)
        {
          float rowa[8 * 4], rowb[8 * 4];
          for(int k = 0; k < 8 * 4; k++) rowa[k] = rowb[k] = 1.f;
          lf_modifier_apply_color_modification(vref, rowa, xs[xi], H * 0.25f, 8, 1,
                                               LF_CR_4(RED, GREEN, BLUE, UNKNOWN), 8 * 4 * sizeof(float));
          ls_modifier_apply_vignetting(&vmod, xs[xi], H * 0.25f, 8, 1, rowb, 8 * 4 * sizeof(float));
          for(int k = 0; k < 8 * 4; k++)
          {
            /* Alpha is not compared against upstream, because upstream's own two paths
             * disagree about it: the SSE2 DeVignetting multiplies all four components, the
             * scalar one -- which this harness forces via the CPU-features interposition --
             * leaves the fourth alone. There is no single upstream answer to be faithful to.
             *
             * LensSerious leaves it alone, deliberately: the fourth component is not a
             * colour, so there is no falloff in it to remove, and a consumer keeping a mask
             * or premultiplied coverage there would have it silently corrupted. That is a
             * decision this library makes rather than inherits, so it is ASSERTED here
             * rather than merely excluded -- the row went in at 1.0 and must come back at
             * exactly 1.0. */
            if((k & 3) == 3)
            {
              if(rowb[k] != 1.f)
              {
                if(vig_failed < 3)
                  printf("VIGFAIL %s crop=%.2f focal=%.1f x=%.0f k=%d  alpha scaled to %.6f,"
                         " must be untouched\n", lf->Model, crop, focal, xs[xi], k, rowb[k]);
                vig_failed++;
              }
              continue;
            }

            const float d = fabsf(rowa[k] - rowb[k]);
            if(d > worst_vig)
            {
              worst_vig = d;
              snprintf(worst_vig_name, sizeof(worst_vig_name), "%s @ %.1fmm", lf->Model, focal);
            }
            if(d > TOL_VIGNETTING)
            {
              if(vig_failed < 3)
                printf("VIGFAIL %s crop=%.2f focal=%.1f x=%.0f k=%d  lf=%.6f ls=%.6f  (lfmods=%d lsmods=%d)\n",
                       lf->Model, crop, focal, xs[xi], k, rowa[k], rowb[k], vrefmods, vmymods);
              vig_failed++;
            }
          }
        }
      }
      lf_modifier_destroy(vref);
      lf_modifier_destroy(ref);
    }
  }
  /* GetLenses() returns the database-owned array: not ours to free. */
  lf_db_destroy(ldb);

  printf("parity: %d lenses total, %d geometry-compared, %d skipped (projection not yet"
         " offered), %d vignetting-compared\n", total, compared, skipped_geometry, vig_compared);
  printf("parity: worst geometry delta %.6f px (%s); everything held to upstream's SCALAR"
         " answer at %.3f px -> %s\n",
         worst_px, worst_name, TOL_GEOMETRY_PX, failed ? "FAIL" : "pass");
  printf("parity: upstream's own row walk deviates %.6f px from its own width-1 answer (%s)"
         " -- %s\n", upstream_row_dev, upstream_name[0] ? upstream_name : "none",
         upstream_row_dev > TOL_GEOMETRY_PX
             ? "liblensfun ran its SSE path here; that figure is the approximation every"
               " lensfun-corrected render carries, and is not a verdict on this library"
             : "liblensfun ran its scalar path here");
  printf("parity: %d pixels declined as past a reversed projection's horizon (upstream"
         " answers with a wrapped tangent there); %d of them outside that case -> %s\n",
         declined, bad_decline, bad_decline ? "FAIL" : "ok");
  printf("parity: %d samples judged on the relative rule at %g (reversing a projection near"
         " its horizon, where dr/dtheta = f*sec^2 amplifies a one-ulp theta)\n",
         divergent, TOL_GEOMETRY_REL);
  printf("parity: %d oflag mismatches against upstream (axes refused on one side only)\n",
         flag_mismatch);
  printf("parity: autoscale over %d configurations, worst %.2e relative (%s), tolerance %g"
         " -> %s\n", autoscale_compared, autoscale_worst,
         autoscale_name[0] ? autoscale_name : "none", TOL_AUTOSCALE,
         autoscale_failed ? "FAIL" : "pass");
  printf("parity: %d autoscale configurations differ where a reversed projection's horizon"
         " falls inside the frame (no scale removes those borders); %d elsewhere\n",
         autoscale_horizon, autoscale_failed);
  printf("parity: worst vignetting delta %.6f (%s), tolerance %g -> %s\n",
         worst_vig, worst_vig_name, TOL_VIGNETTING, vig_failed ? "FAIL" : "pass");
  if(flag_mismatch) failed += flag_mismatch;
  if(bad_decline) failed += bad_decline;
  if(autoscale_failed) failed += autoscale_failed;
  if(failed || vig_failed)
  {
    printf("parity: %d geometry samples, %d vignetting samples out of tolerance\n", failed, vig_failed);
    return 1;
  }
  return 0;
}
