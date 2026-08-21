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

Against liblensfun 0.3.4, 24 Mpx, Nikon AF-S 16-35mm f/4G ED VR, on an idle 8-core Xeon
E3-1505M ([`tests/bench_lensfun.c`](tests/bench_lensfun.c)). Each stage runs five times and
the best is kept — a shared machine only ever makes a measurement slower, so the minimum is
the closest thing to the cost of the code. Absolute numbers still drift a few percent
between sessions; the ratios are what to read.

The two halves are separated below because they behave completely differently, and averaging
them into one "how fast is it" number hides the only result that matters. Getting the
calibration data is where LensSerious wins by two orders of magnitude and always will.
Processing the pixels is where it is *level with lensfun* on a CPU — the win there is
structural, and it is [further down](#readme-fused).

### Getting the calibration data

Paid on the CPU whatever the pixels are processed on: open a database, resolve a free-text
lens name, then resolve that lens at the shot's focal and aperture.

| | lensfun | LensSerious | |
|---|---|---|---|
| open the database | 101.7 ms | **0.18 ms** | 576× faster |
| find the lens (fuzzy) | 0.031 ms | 0.174 ms | **5.5× slower** |
| resolve at focal/aperture | 0.44 µs | ~0.03 µs | >10× faster |
| **= open + find + resolve** | **101.75 ms** | **0.35 ms** | **289× faster** |

The total is the row worth reading, and it is 289× because of one term. lensfun's 101.7 ms
is almost entirely the XML parse: it reads and DOM-parses the whole upstream database on
open, so the first lens costs a tenth of a second before any lens has been looked at.
LensSerious opens an SQLite file that was already parsed offline, `mmap`s it read-only and
immutable, and reads what the query touches — 0.18 ms.

That is also why the fuzzy match being **5.5× slower is not a defect worth fixing**. It is
0.174 ms against a 101.7 ms parse: even 20× slower matching would be invisible next to the
term it replaced. It is the honest number and it stays in the table, but the useful
comparison is the total, not the row. (It got there from 2.89 ms — the three approaches that
did *not* work are in [the maintainer's log](doc/maintainers.md#log-matcher).)

Two structural consequences follow from that row, and they are the actual reason for the
rewrite:

- lensfun's cost is per *process*, so it is paid again by every worker, and it scales with
  the size of the upstream database rather than with the number of lenses used. LensSerious'
  is per *query*.
- Because the file is immutable and lock-free, that 0.35 ms is also what it costs from
  **any** thread concurrently, with no mutex and no shared handle.

### Processing the pixels

Paid per image, and the part a GPU can take.

| single-threaded, 24 Mpx | lensfun | LensSerious | |
|---|---|---|---|
| build the whole geometry map | 296.0 ms | 293.0 ms | **1.0× — no faster** |
| apply vignetting, whole frame | 71.8 ms | **47.6 ms** | 1.5× faster |

The map row is **1.0×**, and that is worth stating plainly rather than burying: the ~300 ms
was always *lensfun's* cost, and closed-form C pays essentially the same to push six floats
per pixel through memory. Vignetting is 1.5× only because upstream leaves three divides and
a vectorisation-blocking branch per pixel on the table; the arithmetic is the same
arithmetic.

**A consumer's compiler flags matter more than anything in this library.** The same map is
286 ms at plain `-O3`, 217 ms with `-march=native`, and 103 ms with `-ffast-math`. Ansel
builds with the last, so its map costs ~100 ms, not ~300. (That also means its CPU results
are no longer bit-identical to its GPU ones; see
[Fused evaluation](doc/fused-evaluation.md#fused-exactness).)

#### In parallel

A pixel pipeline runs the per-pixel stages across threads, so they are measured that way
too. Rows are independent on both sides — lensfun's modifier is const during apply and
writes only the caller's buffer, which is what lets Ansel parallelise it — so this is the
same work, not one library's threading model against another's.

| 8 threads, 24 Mpx | lensfun | LensSerious | |
|---|---|---|---|
| build the whole geometry map | 91.7 ms | **80.7 ms** | 1.1× |
| apply vignetting, whole frame | 32.2 ms | 33.3 ms | 1.0× |
| scaling vs. its own 1 thread | 3.23× | 3.63× | |

Neither scales linearly, and vignetting converges to ~33 ms for both: 732 MB read and
written in 0.033 s is about 22 GB/s, which is the memory bus rather than either library. It
is also why the 1.5× single-threaded advantage disappears — LensSerious is already close to
the ceiling with one thread.

**This is the finding, not a disappointment.** Once both libraries are bandwidth-bound
moving the same map through the same memory, no amount of arithmetic wins. The way past it
is to stop moving the map.

### One image, end to end

| | lensfun | LensSerious | |
|---|---|---|---|
| one image, cold start | 473.3 ms | **343.6 ms** | 1.4× |
| one more image, warm | 369.8 ms | **343.1 ms** | 1.1× |

The gap between the two rows is the whole database story: 103 ms of it is lensfun's XML
parse, and it disappears from the second image because that cost is per-process. The warm
row — 1.1× — is what a batch export actually converges to on a CPU with the map, and it is
the honest steady-state number for this library as a drop-in replacement.

## The map, versus never building one  {#readme-fused}

A correction is only ever wanted so that something can **resample** with it, and that is
where the two libraries actually differ. lensfun can only deliver a correction as a buffer:
its callbacks fill six floats per pixel and the consumer reads them back. LensSerious can be
called per pixel, so a consumer can evaluate the coordinates where it needs them and never
materialise anything — which is what the OpenCL kernels do, and what the CPU can do equally
well.

All three produce the same image. The difference is 1.1 GB of memory traffic per frame: 549
MB written by the map pass, 549 MB read back by the resampler.

| correct + resample, 24 Mpx | 1 thread | 2 | 4 | 8 |
|---|---|---|---|---|
| lensfun: map, then resample | 765 ms | 355 | 237 | 236 |
| LensSerious: map, then resample | 744 ms | 347 | 229 | 218 |
| LensSerious: **fused, no map** | 874 ms | 438 | 227 | **174** |
| fused vs. lensfun | 0.88× | 0.81× | 1.04× | **1.35×** |

**Fusing is a loss until four threads and a win after**, and the crossover is the point. With
one thread the work is latency-bound: two tight loops each optimise and prefetch better than
one fused loop that interleaves evaluation with scattered reads. With eight, the two-pass
form saturates the memory bus and stops scaling — lensfun goes 237 → 236 ms from four threads
to eight, gaining *nothing* — while the fused pass has no map traffic to saturate on and keeps
scaling.

(That table is best-of-three with the thread counts interleaved rather than swept in order.
Run as a plain sweep it reports 1.1× at eight threads instead of 1.35×, because the low-thread
runs heat the machine first and the memory-bound row degrades more than the compute-bound one.
The ordering is part of the measurement.)

That is the same effect the GPU shows, at the other end of the scale: a device with thousands
of threads is so far past the crossover that not building the map is worth 165× (below).

The consequences for a consumer — always fuse on a GPU, measure before fusing on a CPU, and
why there is deliberately no `ls_correct_and_resample()` in the API — are in
[Fused evaluation](doc/fused-evaluation.md).

## What it costs to reach a GPU

This is where the design claim actually lives, and it is not about arithmetic
([`tests/bench_opencl.c`](tests/bench_opencl.c), Quadro M2200, 24 Mpx):

| geometry map, 549 MB of coordinates | |
|---|---|
| lensfun: build on the CPU | 481.7 ms |
| lensfun: upload it | 76.5 ms |
| **lensfun: total to reach the device** | **558.2 ms** |
| LensSerious: the same map, on the device | 35.9 ms — **15.5×** |
| LensSerious: evaluated, never stored | 3.4 ms — **165×** |

| vignetting, 366 MB of multipliers | |
|---|---|
| **lensfun: total to reach the device** | **258.7 ms** |
| LensSerious: on the device | 1.7 ms — **153×** |

The CPU build is 482 ms here against 300 ms in the table above. Same single-threaded loop —
the difference is that this one materialises the whole 549 MB buffer, while the other reuses
one row. 482 ms is what a GPU pipeline actually pays, because it needs the entire map before
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
| map-less (fused) evaluation on the CPU | measured; a win from four threads up, a loss below — see above. Not yet offered as an API: consumers call `ls_eval_map()` per pixel themselves, as `tests/bench_lensfun.c` does |
| upstream sync (`version_2` conversion tooling) | not started |
| native XML reader (dropping liblensfun from the importer) | not started |
| CPU geometry map vectorisation | **open** — the loop is scalar because the model dispatch inside the per-pixel evaluator is loop-invariant and the compiler will not unswitch it |
| `<real-focal-length>` interpolation | **open** — would let the last projection cases off the fallback path; the schema already has the table, the importer does not fill it |

## Documentation

The API reference and the design notes are generated with Doxygen and published to GitHub
Pages by `.github/workflows/docs.yml`. To build them locally:

    git submodule update --init doc/doxygen-awesome-css
    doxygen doc/Doxyfile          # writes doc/api/html

Two pages there are worth reading before changing anything:

- **[Fused evaluation](doc/fused-evaluation.md)** — why the displacement map should usually
  not exist, and the thread count at which that stops being true.
- **[Maintainer's log](doc/maintainers.md)** — the experiments that were measured and
  reverted, each with the number that killed it. Read it before optimising anything here;
  most of the obvious ideas are in it already.

## License

LGPL-3.0-or-later. Contains code ported from Lensfun 0.3.4,
Copyright (C) 2007 Andrew Zabolotny and the Lensfun contributors.
The lens calibration data consumed at runtime is the lensfun database,
distributed under its own terms (CC-BY-SA), and is not part of this repository.
