Fused evaluation: why the map should usually not exist  {#fused-evaluation}
======================================================

[TOC]

A lens correction is only ever wanted so that something can **resample** with it. That
sentence is the whole design, and it is worth stating before any code: the displacement map
is not the product, it is an intermediate that one particular implementation happens to
need.

lensfun needs it. Its corrections are delivered through callbacks that fill a buffer —
`lfModifier::ApplySubpixelGeometryDistortion()` writes six floats for every output pixel, and
the consumer reads them back to sample with. There is no way to ask lensfun for "the source
coordinate of this one pixel" cheaply, so the buffer is mandatory.

LensSerious has no such constraint. @ref ls_eval_map() is a pure function of an @ref
ls_eval_t and a coordinate, so a consumer can call it *where it needs the answer* and never
write a map at all. On a GPU that is exactly what happens: the evaluator is `#include`d into
the kernel that consumes the coordinates, and the correction crosses the bus as ~80 bytes of
coefficients instead of six floats per pixel.

@section fused-cost What the map actually costs

For a 24 Mpx frame, six floats per pixel is **549 MB**. A two-pass pipeline writes all of it
and then reads all of it back: 1.1 GB of traffic that a fused pass simply does not generate.

On a GPU that is the difference between 558 ms and 3.4 ms — the map has to be built on the
CPU (482 ms, single-threaded, because lensfun's callbacks are) and then uploaded (76 ms).
Evaluating in the consuming kernel is 165× cheaper, and the ratio is not subtle because
neither term survives.

@section fused-cpu On the CPU it depends on the thread count

This is the part that is easy to get wrong, and the reason this page exists.

The same fusion on the CPU is **not** an unconditional win. Measured with a bilinear
resampler over a 24 Mpx frame (`tests/bench_lensfun.c`), all three producing the same image:

| correct + resample | 1 thread | 2 | 4 | 8 |
|---|---|---|---|---|
| lensfun: map, then resample | 765 ms | 355 | 237 | 236 |
| LensSerious: map, then resample | 744 ms | 347 | 229 | 218 |
| LensSerious: **fused, no map** | 874 ms | 438 | 227 | **174** |
| fused vs. lensfun | 0.88× | 0.81× | 1.04× | **1.35×** |

@note Best of three passes with the thread counts **interleaved**, not swept in order. Swept,
the same binary reports 1.1x at eight threads: the low-thread passes heat the machine, and
the two-pass row -- being memory-bound -- loses more to the lower clock than the fused one
does. Measuring a ratio at one thread count in isolation is not enough here.

Fusing **loses** below four threads and wins above, crossing over at four:

- With one thread the work is latency-bound. Two tight loops each vectorise and prefetch
  well; one fused loop interleaves coordinate evaluation with scattered source reads, and
  neither half runs at its best.
- With eight, the two-pass form saturates the memory bus and stops scaling — lensfun goes
  237 ms to 236 ms from four threads to eight, gaining *nothing* — while the fused pass has
  no map traffic to saturate on and keeps scaling.

A GPU sits so far past that crossover that the question does not arise. A CPU sits on either
side of it depending on how many cores the caller gave the job.

@section fused-guidance What this means for a consumer

- **On a GPU, always fuse.** Include `lensserious_eval.h` in the kernel that consumes the
  coordinates and call @ref ls_eval_map() per work-item. `opencl/lensserious.cl` keeps
  map-writing kernels only so a host can verify the device against the CPU; they are not the
  way to use this.
- **On a CPU, measure at the thread count you will actually run at.** Do not assume the GPU
  result transfers. Below the crossover the map is faster, and the crossover depends on the
  machine's memory bandwidth and on how expensive the resampling filter is — a costlier
  filter moves it left, since it raises the compute per byte.
- There is deliberately **no `ls_correct_and_resample()`** in the API. A convenience function
  that fused for you would look like a free win and would be a loss on a two-core machine.
  The library gives you a per-pixel evaluator; the loop structure is the caller's decision,
  because only the caller knows the thread count and the filter.

@section fused-exactness Fusing does not change the numbers

@ref ls_eval_map() evaluates one pixel from absolute coordinates and has no state, so a
fused loop, a row walker and a GPU work-item all pass identical operands to identical
expressions. That is why `ls_modifier_apply_subpixel_geometry()` takes absolute `xu`/`yu`
and sums `+ col` at the call site rather than stepping an accumulator: an incremental step
would be cheaper and would round differently, and the CPU and the kernel would stop agreeing.

The one caveat is not in this library: a consumer built with `-ffast-math` (Ansel is) lets
the compiler reassociate and contract the CPU side, so its results are no longer bit-identical
to the GPU's. The parity harness is built without it, and that is where the bit-exactness
claim is asserted.
