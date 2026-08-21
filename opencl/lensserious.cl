/*
    LensSerious — lens correction evaluated inside the kernel that consumes it.

    Copyright (C) 2026 Aurélien PIERRE.  License: LGPL-3.0-or-later.

    Lensfun can only produce its displacement map through single-threaded CPU callbacks
    (measured 278 ms per 24 Mpx frame), after which a GPU pipeline uploads six floats per
    pixel -- 576 MB for that same frame -- purely so a resampler can read them back. Here
    the correction crosses as an ls_eval_t of scalars and each work-item evaluates its own
    coordinates in a handful of FMAs.

    THE MATH IS NOT IN THIS FILE. It is in include/lensserious_eval.h, which the CPU
    library includes too, so there is one source text and no second copy to drift. What
    lives here is the OpenCL-only part: the address-space plumbing and the entry points.

    A consuming kernel should not call ls_map_pixel_buffer() below -- writing a map to
    memory is the cost this library exists to remove. It should include this file and call
    ls_eval_map()/ls_eval_vignette_factor() directly at the point it needs a source
    coordinate. The kernels here exist so a host can verify the device against the CPU.
*/

#include "lensserious_eval.h"

/**
 * Reference entry point: fill the 6-float-per-pixel map, byte-compatible with
 * lfModifier::ApplySubpixelGeometryDistortion()'s buffer, so a host can compare the two
 * element-wise. Production paths evaluate in place instead; see the note above.
 */
kernel void ls_subpixel_geometry(global float *res,
                                 const int width, const int height,
                                 const float xu, const float yu,
                                 const ls_eval_t p)
{
  const int col = get_global_id(0);
  const int row = get_global_id(1);
  if(col >= width || row >= height) return;

  float out[6];
  ls_eval_map(&p, xu + (float)col, yu + (float)row, out);

  global float *dst = res + ((size_t)row * width + col) * 6;
  for(int k = 0; k < 6; k++) dst[k] = out[k];
}

/** Reference entry point for vignetting: the multiplier per pixel, one float each. */
kernel void ls_vignette_map(global float *res,
                            const int width, const int height,
                            const float xu, const float yu,
                            const ls_eval_t p)
{
  const int col = get_global_id(0);
  const int row = get_global_id(1);
  if(col >= width || row >= height) return;

  res[(size_t)row * width + col] =
      ls_eval_vignette_factor(&p, xu + (float)col, yu + (float)row);
}
