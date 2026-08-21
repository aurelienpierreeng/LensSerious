/*
    LensSerious — the geometry map as an OpenCL kernel.

    Copyright (C) 2026 Aurélien PIERRE.  License: LGPL-3.0-or-later.

    The SAME closed forms as src/lensserious.c, evaluated per work-item. This is the
    point of the whole library: lensfun can only produce this map through single-threaded
    CPU callbacks (measured 278 ms / 24 Mpx frame) which a GPU pipeline must then upload;
    here the map costs a handful of FMAs inside the kernel that consumes it. The host
    fills ls_cl_params_t from a resolved ls_modifier_t — plain values, no handles.
*/

typedef struct ls_cl_params_t
{
  float norm_scale, norm_unscale;
  float center_x, center_y;
  float scale;              /* reciprocal, pre-resolved; 1.0 when disabled */
  int   dist_model;         /* ls_dist_model_t */
  float dist_terms[4];
  int   tca_model;          /* ls_tca_model_t */
  float tca_terms[6];
  int   enabled;            /* LS_ENABLE_* bits, matching lensserious.h */
} ls_cl_params_t;

static float2 ls_dist_eval(const int model, __constant const float *terms, float2 p)
{
  const float ru2 = p.x * p.x + p.y * p.y;
  float m = 1.0f;
  if(model == 1)        /* POLY3: Rd = Ru · (1 − k1 + k1·Ru²) */
    m = (1.0f - terms[0]) + terms[0] * ru2;
  else if(model == 2)   /* POLY5 */
    m = 1.0f + terms[0] * ru2 + terms[1] * ru2 * ru2;
  else if(model == 3)   /* PTLENS */
  {
    const float r = native_sqrt(ru2);
    m = terms[0] * ru2 * r + terms[1] * ru2 + terms[2] * r
        + (1.0f - terms[0] - terms[1] - terms[2]);
  }
  return p * m;
}

/* One work-item = one output pixel; writes the 6-float R/G/B source coordinates,
 * byte-compatible with lfModifier::ApplySubpixelGeometryDistortion's buffer. */
kernel void ls_subpixel_geometry(global float *res,
                                 const int width, const int height,
                                 const float xu, const float yu,
                                 constant ls_cl_params_t *p)
{
  const int col = get_global_id(0);
  const int row = get_global_id(1);
  if(col >= width || row >= height) return;

  float2 c = (float2)((xu + col) * p->norm_scale - p->center_x,
                      (yu + row) * p->norm_scale - p->center_y);

  if(p->enabled & (1 << 3)) c *= p->scale;
  if(p->enabled & (1 << 0)) c = ls_dist_eval(p->dist_model, p->dist_terms, c);

  float2 r = c, b = c;
  if(p->enabled & (1 << 1))
  {
    if(p->tca_model == 1)            /* LINEAR */
    {
      r *= p->tca_terms[0];
      b *= p->tca_terms[1];
    }
    else if(p->tca_model == 2)       /* POLY3: Rd = Ru·(b·Ru² + c·Ru + v) */
    {
      float ru2 = r.x * r.x + r.y * r.y;
      r *= p->tca_terms[4] * ru2 + p->tca_terms[0]
           + (p->tca_terms[2] != 0.0f ? p->tca_terms[2] * native_sqrt(ru2) : 0.0f);
      ru2 = b.x * b.x + b.y * b.y;
      b *= p->tca_terms[5] * ru2 + p->tca_terms[1]
           + (p->tca_terms[3] != 0.0f ? p->tca_terms[3] * native_sqrt(ru2) : 0.0f);
    }
  }

  global float *out = res + ((size_t)row * width + col) * 6;
  out[0] = (r.x + p->center_x) * p->norm_unscale;
  out[1] = (r.y + p->center_y) * p->norm_unscale;
  out[2] = (c.x + p->center_x) * p->norm_unscale;
  out[3] = (c.y + p->center_y) * p->norm_unscale;
  out[4] = (b.x + p->center_x) * p->norm_unscale;
  out[5] = (b.y + p->center_y) * p->norm_unscale;
}
