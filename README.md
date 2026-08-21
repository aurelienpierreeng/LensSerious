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
| open the database | 88.2 ms | **0.166 ms** | 531× faster |
| find the lens (fuzzy) | 0.028 ms | 3.15 ms | **112× slower** |
| resolve at focal/aperture | 0.4 µs | 0.04 µs | 11× faster |
| build the whole geometry map | 271.8 ms | 268.5 ms | **1.0× — no faster** |
| apply vignetting, whole frame | 66.1 ms | 117.3 ms | **1.8× slower** |
| one image, cold start | 426 ms | 389 ms | 1.1× |
| one more image, warm | 338 ms | 386 ms | 0.9× |

**Read that honestly: the closed forms are not faster than lensfun on the CPU.** The 278 ms
quoted above was always *lensfun's* cost, and plain C pays essentially the same to push six
floats per pixel through memory. Two rows are slower on purpose:

- **vignetting**, because `ls_eval_vignette_factor()` is written to be independent per
  work-item — it recomputes the row term and divides for every pixel — which is what a GPU
  wants and what keeps one source text for both. On the CPU that trade costs ~51 ms/frame.
- **the fuzzy lookup**, which is simply a database losing to an in-memory search: ~4700
  rows through a b-tree against structs already in RAM. It is once per image against 268 ms
  of map building. An inverted token index was built for it and *measured*, and is not in
  the code because it made things worse (5.13 ms against 3.15 ms, identical agreement):
  model tokens — `mm`, `f`, `ed`, `vr`, the focal digits — are not selective.

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
| fuzzy matching | done; **99.0%** agreement with liblensfun's own top pick over the whole database (100% verbatim, 99.5% maker-stripped) |
| OpenCL | evaluators compile as OpenCL C from the same header; Ansel consumes them per work-item |
| upstream sync (`version_2` conversion tooling) | not started |
| native XML reader (dropping liblensfun from the importer) | not started |

## License

LGPL-3.0-or-later. Contains code ported from Lensfun 0.3.4,
Copyright (C) 2007 Andrew Zabolotny and the Lensfun contributors.
The lens calibration data consumed at runtime is the lensfun database,
distributed under its own terms (CC-BY-SA), and is not part of this repository.
