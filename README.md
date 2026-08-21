# LensSerious

Lens-correction mathematics as **data**, not as a library of callbacks.

## Why

[Lensfun](https://lensfun.github.io/) carries an excellent community database inside a
runtime that fights a modern raw pipeline at every step:

- the database is 8.4 MB of XML parsed into flat lists on every startup (measured
  89–102 ms), searched by linear fuzzy scans, with `setlocale()` flipped during parsing
  (process-global, a real concurrency hazard);
- the displacement map can only be produced by single-threaded C++ callbacks — measured
  **278 ms per 24 Mpx frame** — after which a GPU pipeline must upload six floats per
  pixel just to resample them.

Yet the actual mathematics is six small closed forms. A census of the complete
`version_1` database (2026): distortion **ptlens/poly3/poly5** (4810/875/5 entries),
TCA **poly3/linear** (3355/6), vignetting **pa** (25269 — the only model). LensSerious
expresses those forms as plain C over a plain struct of coefficients, so the same
evaluation runs vectorised on the CPU, inside an OpenCL kernel, or anywhere a float goes.

Not "the same math, textually" — **the same text**.
[`include/lensserious_eval.h`](include/lensserious_eval.h) holds the evaluators, and it
is compiled twice: as C99 into the library, and as OpenCL C by the host's driver. There
is no second copy to drift, which matters because the two copies this replaced already
had: the kernel computed its square roots with `native_sqrt` — a 12-bit approximation,
against the library's `sqrtf` — and never grew a vignetting evaluator at all. Neither
was visible to a harness that only tests the CPU side. `ctest` now compiles the kernel
offline through clang's OpenCL mode, so the second compilation is checked like the first.

A correction therefore crosses to the GPU as an
[`ls_eval_t`](include/lensserious_eval.h) — about 80 bytes of scalars, passed by value as
a kernel argument — and each work-item evaluates its own coordinates. Not a
six-float-per-pixel map: that is 576 MB of transfer for a 24 Mpx frame, on top of the
278 ms spent building it. A consuming kernel `#include`s the header and calls
`ls_eval_map()` where it needs a source coordinate;
[`opencl/lensserious.cl`](opencl/lensserious.cl) keeps map-writing kernels only so a host
can verify the device against the CPU.

**The lensfun project is not being forked.** Its XML database remains the interchange
format and its community remains the source of calibrations. Only the runtime is
replaced.

## Parity is measured, never assumed

Every convention is ported from lensfun 0.3.4 source *with its quirks kept*: the
`(width − 1)` sizing, the calibration-sensor crop/aspect correction chain, the Hermite
interpolation with one-sided tangents. A "cleaner" convention would silently shift every
corrected render ever made.

Two instruments hold that line:

- [`tests/parity_lensfun.c`](tests/parity_lensfun.c) links **both** libraries and walks
  the entire installed database: every lens × three focal lengths × a 5×5 sample grid,
  geometry asserted to 0.01 px and vignetting to 1e-4. Lenses whose projection is not
  implemented yet are **skipped and counted**, never silently passed.
- In [Ansel](https://github.com/aurelienpierreeng/ansel), `iop/lens.cc` latches against
  lensfun **and** LensSerious simultaneously and logs the live per-frame deviation, so
  the comparison runs on real raws in the real pipeline until the day the two part ways.

## What the parity work found in upstream

Two genuine discoveries, both invisible to any width-1 test and caught only by comparing
full rows in the real pipeline:

- **Lensfun's SSE row path is ~0.2–0.8 px approximate.** Its SIMD variants compute
  `sqrt(r²)` as `_mm_rcp_ps(_mm_rsqrt_ps(r²))` — two chained 12-bit approximations with no
  Newton step (`mod-coord-sse.cpp`) — and disagree with lensfun's *own scalar math* by up
  to 0.8 px at the end of wide rows on strong wide-angle lenses. Every lensfun-corrected
  render carries that error. LensSerious evaluates exactly (identically on CPU and GPU) and
  matches upstream's **scalar** semantics to < 0.01 px over the whole database; the parity
  harness forces upstream onto its scalar path by interposing its CPU detection.
- **Lensfun's two vignetting paths disagree about alpha.** The SSE2 `DeVignetting`
  multiplies all four components; the scalar one leaves the fourth alone. LensSerious
  matches the SSE2 behaviour, which is what production runs on every x86-64.

## Status

| piece | state |
|---|---|
| distortion (poly3/poly5/ptlens), TCA (linear/poly3), scaling | ported, harness-gated |
| vignetting (pa) | ported; the (focal, aperture, distance) interpolation is inverse-distance weighted and flagged by the harness where it drifts from upstream |
| projection conversion (fisheye ↔ rectilinear, …) | **not yet** — reported via `geometry_unsupported`, callers fall back |
| database (XML → indexed store) & fuzzy matching | phase 2 |
| upstream sync (`version_2` conversion tooling) | phase 3 |

## License

LGPL-3.0-or-later. Contains code ported from Lensfun 0.3.4,
Copyright (C) 2007 Andrew Zabolotny and the Lensfun contributors.
The lens calibration data consumed at runtime is the lensfun database,
distributed under its own terms (CC-BY-SA), and is not part of this repository.
