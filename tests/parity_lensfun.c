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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TOL_GEOMETRY_PX 0.01f   /* max |Δ| in pixels over the sample grid */
#define TOL_VIGNETTING  1e-4f   /* max |Δ| on the multiplier */

static int _lens_from_lf(const lfLens *lf, ls_lens_t *out)
{
  memset(out, 0, sizeof(*out));
  switch(lf->Type)
  {
    case LF_RECTILINEAR: out->type = LS_LENS_RECTILINEAR; break;
    case LF_UNKNOWN:     out->type = LS_LENS_UNKNOWN; break;
    default:             out->type = LS_LENS_FISHEYE; break; /* any non-rectilinear: skip */
  }
  out->crop_factor = lf->CropFactor;
  out->min_focal = lf->MinFocal;
  out->max_focal = lf->MaxFocal;
  out->aspect_ratio = lf->AspectRatio;
  out->center_x = lf->CenterX;
  out->center_y = lf->CenterY;

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

  const lfLens **lenses = lf_db_get_lenses(ldb);
  int total = 0, compared = 0, skipped_geometry = 0, failed = 0;
  int vig_compared = 0, vig_failed = 0;
  float worst_px = 0.f, worst_vig = 0.f;
  char worst_name[256] = "", worst_vig_name[256] = "";

  for(int li = 0; lenses && lenses[li]; li++)
  {
    const lfLens *lf = lenses[li];
    total++;

    ls_lens_t lens;
    if(!_lens_from_lf(lf, &lens)) continue;
    if(lens.type != LS_LENS_RECTILINEAR && lens.type != LS_LENS_UNKNOWN)
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

    for(int ci = 0; ci < 2; ci++)
    for(int fi = 0; fi < 3; fi++)
    {
      const float crop = lf->CropFactor * crop_ratios[ci];
      const float focal = focals[fi];

      lfModifier *ref = lf_modifier_new(lf, crop, W, H);
      const int refmods = lf_modifier_initialize(ref, lf, LF_PF_F32, focal, 8.f, 1000.f, 1.f,
                                                 LF_RECTILINEAR,
                                                 LF_MODIFY_DISTORTION | LF_MODIFY_TCA, 0);
      ls_modifier_t mod;
      const int mymods = ls_modifier_init(&mod, &lens, crop, W, H, focal, 8.f, 1000.f, 1.f,
                                          LS_ENABLE_DISTORTION | LS_ENABLE_TCA);

      if((refmods & (LF_MODIFY_DISTORTION | LF_MODIFY_TCA)) && mymods)
      {
        compared++;
        /* Border + centre + diagonal samples: distortion is worst at corners. */
        const float xs[5] = { 0.f, W * 0.25f, W * 0.5f, W * 0.75f, W - 1.f };
        const float ys[5] = { 0.f, H * 0.25f, H * 0.5f, H * 0.75f, H - 1.f };
        for(int yi = 0; yi < 5; yi++)
          for(int xi = 0; xi < 5; xi++)
          {
            float a[6], b[6];
            lf_modifier_apply_subpixel_geometry_distortion(ref, xs[xi], ys[yi], 1, 1, a);
            ls_modifier_apply_subpixel_geometry(&mod, xs[xi], ys[yi], 1, 1, b);
            for(int k = 0; k < 6; k++)
            {
              const float d = fabsf(a[k] - b[k]);
              if(d > worst_px)
              {
                worst_px = d;
                snprintf(worst_name, sizeof(worst_name), "%s @ %.1fmm term %d", lf->Model, focal, k);
              }
              if(d > TOL_GEOMETRY_PX)
              {
                if(failed < 3)
                  printf("GEOFAIL %s crop=%.2f focal=%.1f xy=(%.0f,%.0f) k=%d  lf=%.4f ls=%.4f\n",
                         lf->Model, crop, focal, xs[xi], ys[yi], k, a[k], b[k]);
                failed++;
              }
            }
          }
      }

      /* Vignetting parity, where both sides resolved the pa model. */
      lfModifier *vref = lf_modifier_new(lf, crop, W, H);
      const int vrefmods = lf_modifier_initialize(vref, lf, LF_PF_F32, focal, 8.f, 1000.f, 1.f,
                                                  LF_RECTILINEAR, LF_MODIFY_VIGNETTING, 0);
      ls_modifier_t vmod;
      const int vmymods = ls_modifier_init(&vmod, &lens, crop, W, H, focal, 8.f, 1000.f, 1.f,
                                           LS_ENABLE_VIGNETTING);
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
         " implemented), %d vignetting-compared\n", total, compared, skipped_geometry, vig_compared);
  printf("parity: worst geometry delta %.6f px (%s), tolerance %.3f px -> %s\n",
         worst_px, worst_name, TOL_GEOMETRY_PX, failed ? "FAIL" : "pass");
  printf("parity: worst vignetting delta %.6f (%s), tolerance %g -> %s\n",
         worst_vig, worst_vig_name, TOL_VIGNETTING, vig_failed ? "FAIL" : "pass");
  if(failed || vig_failed)
  {
    printf("parity: %d geometry samples, %d vignetting samples out of tolerance\n", failed, vig_failed);
    return 1;
  }
  return 0;
}
