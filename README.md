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

Four instruments hold that line:

- [`tests/parity_lensfun.c`](tests/parity_lensfun.c) links **both** libraries and walks
  the entire installed database: every lens × three focal lengths × a 5×5 sample grid,
  geometry asserted to 0.01 px and vignetting to 1e-4. Configurations LensSerious declines
  to serve are **skipped and counted**, never silently passed.
- [`tests/parity_db.c`](tests/parity_db.c) checks every field of every lens the importer
  wrote against liblensfun, *exactly*, and runs eight threads with their own handles over
  the same file to exercise the no-mutex contract.
- [`tests/match_lensfun.c`](tests/match_lensfun.c) uses liblensfun as an **oracle** for the
  fuzzy matcher — an installed lensfun ships no sources, so what is asserted is that the
  two reach the same answer, not that they run the same algorithm.
- In [Ansel](https://github.com/aurelienpierreeng/ansel), `iop/lens.cc` latches against
  lensfun **and** LensSerious simultaneously and logs the live per-frame deviation, so
  the comparison runs on real raws in the real pipeline until the day the two part ways.

## What the parity work found in upstream

Three genuine findings, none of them visible to a width-1 test — the first two were caught
by comparing full rows in the real pipeline, the third by including fisheye lenses in the
harness for the first time:

- **Lensfun's SSE row path is ~0.2–0.8 px approximate.** Its SIMD variants compute
  `sqrt(r²)` as `_mm_rcp_ps(_mm_rsqrt_ps(r²))` — two chained 12-bit approximations with no
  Newton step (`mod-coord-sse.cpp`) — and disagree with lensfun's *own scalar math* by up
  to 0.8 px at the end of wide rows on strong wide-angle lenses. Every lensfun-corrected
  render carries that error. LensSerious evaluates exactly (identically on CPU and GPU) and
  matches upstream's **scalar** semantics to < 0.01 px over the whole database; the parity
  harness forces upstream onto its scalar path by interposing its CPU detection.
- **The `pa` vignetting polynomial is not constrained to stay positive**, and for some
  lenses it crosses zero *inside the frame* — the Canon EF 8-15mm Fisheye at 8mm has
  k = (−0.625, 5.648, −19.330), whose root sits near r = 0.65. Past it the correction 1/c
  is negative, which is not a brightness. Upstream clamps it at `apply_multiplier()`;
  LensSerious clamps the factor, which for non-negative pixels is the same thing one step
  earlier.
- **Lensfun's two vignetting paths disagree about alpha.** The SSE2 `DeVignetting`
  multiplies all four components; the scalar one leaves the fourth alone. LensSerious
  matches the SSE2 behaviour, which is what production runs on every x86-64.

## Design

Three layers, and only the first ends up inside a pixel pipeline.

**The evaluators** — [`include/lensserious_eval.h`](include/lensserious_eval.h). Six closed
forms plus a projection stage, over `ls_eval_t`: a flat block of ~30 scalars, no pointers,
no vectors, laid out identically under every host C ABI and under OpenCL C. It is compiled
**twice** — as C99 into `liblensserious`, and as OpenCL C by the host's driver — so there is
one source text and no second copy to drift. That is not a style preference: the two
hand-maintained copies this replaced had already drifted, one computing its square roots
with `native_sqrt` (~12 bits, against the library's `sqrtf`) and never growing a vignetting
evaluator at all, with nothing in the harness able to see either. `ctest` compiles the
kernel offline through clang's OpenCL mode so the second compilation is checked like the
first. Everything transcendental goes through `LS_SQRT`/`LS_ATAN`/… so the host branch gets
the `f`-suffixed names: bare `atan(float)` in C promotes to double, evaluates in double and
truncates back, which is both slower and enough to stop the loops vectorising.

**The resolver** — [`src/lensserious.c`](src/lensserious.c). `ls_modifier_init()` turns a
lens plus a shooting configuration into that scalar block: interpolating the calibration
across focal, aperture and distance, building the coordinate system, and deciding which
axes it can actually serve. Plain values in, plain values out; nothing owned, nothing
locked.

**The database** — [`include/lensserious_db.h`](include/lensserious_db.h) over SQLite,
built offline by [`tools/import_lensfun_xml.c`](tools/import_lensfun_xml.c). Stateless: no
initialisation, no shutdown, no global, nothing cached between calls, and an `ls_lens_t`
that stays valid after the handle is closed. Read-only and genuinely lock-free — opened
`mode=ro&immutable=1` with `SQLITE_OPEN_NOMUTEX`, so SQLite takes no file locks, no
shared-memory segment and no mutex, and neither does this library. The price is one rule:
**one handle per thread**, and a database file is **replaced by rename, never edited**.

The importer reads the XML *through liblensfun*, deliberately. It is an offline tool run
when upstream publishes calibrations, so its dependencies cost a rendered frame nothing —
and it buys `GuessParameters()`, which infers the focal and aperture ranges that are absent
from the XML for ~95% of lenses and that feed the vignetting interpolation's distance
metric. Nothing in the schema or the reader knows liblensfun exists, so a native XML reader
later replaces one file.

## Measurements

Head-to-head against liblensfun 0.3.4, 24 Mpx, single-threaded, both sides on the CPU
([`tests/bench_lensfun.c`](tests/bench_lensfun.c), Nikon AF-S 16-35mm f/4G ED VR):

| | lensfun | LensSerious | |
|---|---|---|---|
| open the database | 99.2 ms | **0.261 ms** | 380× faster |
| find the lens (fuzzy) | 0.038 ms | 0.28 ms | 7× slower |
| resolve at focal/aperture | 0.5 µs | 0.04 µs | 13× faster |
| build the whole geometry map | 288.5 ms | 293.4 ms | **1.0× — no faster** |
| apply vignetting, whole frame | 71.4 ms | **49.3 ms** | 1.4× faster |
| one image, cold start | 459.2 ms | **343.2 ms** | 1.3× |
| one more image, warm | 359.9 ms | **342.7 ms** | 1.1× |

**The geometry map is not faster on the CPU**, and that is worth stating plainly: the 278 ms
quoted above was always *lensfun's* cost, and closed-form C pays essentially the same to
push six floats per pixel through memory. That loop does not vectorise either — the model
dispatch inside the per-pixel evaluator is loop-invariant but the compiler will not unswitch
it, so the generated code is 10 scalar divides and 10 scalar square roots. Real work left on
the table, and the row that matters least, because on a GPU that map is never built.

The **fuzzy lookup** was 2.89 ms — 96× slower — until it was profiled. It was not the
database losing to an in-memory search, it was two mistakes. It scored every one of the
~4700 catalogue names on every lookup, where lensfun rejects almost every candidate on two
float compares before touching a string; and `lens_name` was indexed on `norm` and nothing
else, so the query gathering candidates scanned the whole table anyway.

The fix is to score fewer names, not to score faster: `lens_token` and `token_df` let a
lookup ask which of the *query's* tokens is rarest — `16` or `nikkor` reaches tens of lenses
where `mm` or `f` reaches thousands — and gather only those, falling back to the full scan
if that finds nothing, so the answer never depends on the pruning. 2.89 ms → 0.19 ms, with
agreement unchanged at 99.0%.

Two things were tried, measured, and are deliberately *not* in the code: gathering on **all**
the query's tokens (no better than the scan it replaces), and hashing the tokens to compare
them faster (moved nothing — the per-comparison cost was already ~5 ns; there were simply
half a million comparisons).

**Vignetting** was 1.8× *slower* until the divides were counted. Three causes, worth
recording because two of them are not about arithmetic: the half-diagonal coordinate system
was computed with a divide per axis per pixel (now a folded multiply-and-subtract carried in
`ls_eval_t`); the row-invariant y term was recomputed per pixel; and a `c == 0` guard —
which was not upstream behaviour either — was a second conditional that stopped the caller's
loop vectorising **entirely**, leaving `divss` where there should be `divps`. Removing it and
splitting the row into fill-then-scale gives four divides per instruction. 117 ms → 43 ms.

## What it costs to reach a GPU

This is where the design claim actually lives, and it is not about arithmetic
([`tests/bench_opencl.c`](tests/bench_opencl.c), Quadro M2200, 24 Mpx):

| geometry map, 549 MB of coordinates | |
|---|---|
| lensfun: build on the CPU | 555.5 ms |
| lensfun: upload it | 107.8 ms |
| **lensfun: total to reach the device** | **663.3 ms** |
| LensSerious: the same map, on the device | 37.0 ms — **17.9×** |
| LensSerious: evaluated, never stored | 3.4 ms — **196×** |

| vignetting, 366 MB of multipliers | |
|---|---|
| **lensfun: total to reach the device** | **317.1 ms** |
| LensSerious: on the device | 1.7 ms — **187×** |

The CPU build is 555 ms here against 269 ms in the table above. Same single-threaded loop —
the difference is that this one materialises the whole 549 MB buffer, while the other reuses
one row. 555 ms is what a GPU pipeline actually pays, because it needs the entire map before
it can upload any of it.

The middle row is the apples-to-apples comparison: the same buffer, produced the two
available ways. The last row is the shape production uses, and even that is pessimistic — it
still launches a kernel over the frame, where inlining the evaluator into the kernel that
*consumes* the coordinates makes the map disappear altogether.

Where the design wins is not arithmetic, it is architecture:

- the **88 ms XML parse stops being paid at all**, in every session, including every
  session that never corrects a lens;
- the map **need not exist**. The same evaluator inside the consuming kernel costs a
  handful of FMAs per work-item and no transfer. Measured in Ansel on a Quadro M2200, the
  whole distort + vignette module went **0.501 s → 0.209 s** with ~900 MB of per-frame
  upload removed.

## Status

| piece | state |
|---|---|
| distortion (poly3/poly5/ptlens), TCA (linear/poly3), scaling | ported; 9114 configurations at ≤ 0.0067 px |
| vignetting (pa) | ported; 4054 configurations at ≤ 3e-6 |
| projection conversion (fisheye ↔ rectilinear, orthographic, stereographic, equisolid, thoby) | ported and verified; declines for lenses carrying `<real-focal-length>` data, and for panoramic/equirectangular, which are not radially expressible — 42 configurations |
| database (XML → SQLite) + stateless lock-free reader | done; whole database verified field-by-field against liblensfun |
| fuzzy matching | done; **99.0%** agreement with liblensfun's own top pick over the whole database (100% verbatim, 99.5% maker-stripped), at 0.19 ms a lookup |
| OpenCL | evaluators compile as OpenCL C from the same header; Ansel consumes them per work-item |
| upstream sync (`version_2` conversion tooling) | not started |
| native XML reader (dropping liblensfun from the importer) | not started |
| CPU geometry map vectorisation | **open** — the loop is scalar because the model dispatch inside the per-pixel evaluator is loop-invariant and the compiler will not unswitch it |
| `<real-focal-length>` interpolation | **open** — would let the last projection cases off the fallback path; the schema already has the table, the importer does not fill it |

## License

LGPL-3.0-or-later. Contains code ported from Lensfun 0.3.4,
Copyright (C) 2007 Andrew Zabolotny and the Lensfun contributors.
The lens calibration data consumed at runtime is the lensfun database,
distributed under its own terms (CC-BY-SA), and is not part of this repository.
