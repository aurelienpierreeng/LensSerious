Maintainer's log: what was tried and did not work  {#maintainers}
=================================================

[TOC]

This page is not for users of the library. It is the record of experiments that were run,
measured, and **reverted** — the ones whose absence from the code would otherwise look like
an oversight.

Every entry has the number that killed it. That is the point: each of these is an idea a
competent reader would have on first contact with the code, and rediscovering that it does
not work costs an afternoon. The rule this project runs on is that a plausible explanation
is not a result, and several entries below are cases where a plausible explanation was
confidently wrong.

@section log-matcher The fuzzy matcher

The matcher went from 2.89 ms to 0.19 ms a lookup. Three of the things tried along the way
are not in the code.

**Gathering candidates on all of the query's tokens.** An inverted index over
`lens_token`, unioning every lens that shares any token with the query. Measured 5.13 ms
against 3.15 ms for the plain scan it replaced — *worse*. Model tokens are not selective:
`mm` appears in 1252 of 1562 lenses, `f` in 1240, and the focal digits in hundreds. The
candidate set stayed near the full table and the subquery was pure added cost. What works is
the opposite: gather on the **rarest** token only (`16` reaches 55 lenses), which is what the
code does, with a full-scan fallback so the answer never depends on the pruning.

**Hashing the tokens to compare them faster.** Comparing a 32-bit FNV hash and a length
before `memcmp`. Measured: no change at all. The per-comparison cost was already ~5 ns; there
were simply half a million comparisons. This is the clearest instance in the project of
optimising the inner operation when the problem was the loop count. (The hashes *are* in the
code now, but for a different reason — they let the candidate query read a digest blob out of
a covering index instead of decoding and re-tokenising text.)

**Resolving the brand first**, to split the search space before matching the model. The
database says the opposite: the query's rarest model token reaches 55 lenses where the maker
`nikon` reaches 238, so the token already prunes 4.3× harder. As a *second* filter it is
worse still — an `EXISTS` on the maker returns 12× fewer rows and takes twice as long
(0.144 ms against 0.071 ms), because the correlated subquery costs more per candidate than
the candidates it removes. Requiring the maker would also cost accuracy: it is scored rather
than required precisely because vendors and upstream disagree about vendors' own names
("Nikon Corporation" against "Nikon").

**The prefix-match rule** was removed rather than kept: it scored `0.5 × lmin/lmax` for
tokens sharing three or more leading characters, and deleting it entirely left agreement with
liblensfun at 99.0%, shape for shape. It was costing a `memcmp` per token pair for nothing.

@section log-storage Storage layout

**`WITHOUT ROWID` for the calibration tables was rejected, then adopted.** The first
measurement said it cost 22% of the file size for 0.002 ms — but that comparison was against
a baseline taken *before* `lens_token` and `token_df` existed. Comparing a size across two
different schemas measures the schema difference, not the change. Measured properly, all
three against the same schema:

| | 70-row fetch | file |
|---|---|---|
| `WITHOUT ROWID (lens_id, ord)` | **0.0385 ms** | **3.76 MB** |
| plain table + `index(lens_id)` | 0.0436 ms | 4.19 MB |
| plain table + `index(lens_id, ord)` | 0.0474 ms | 4.35 MB |

Faster *and* smaller — the key is repeated in every row, but dropping three `lens_id` indexes
more than pays for it.

**`SEARCH … USING INDEX` is not the end of the story.** Both matcher queries reported it and
still spent 9% of their time in `sqlite3BtreeTableMoveto`: the row an index lands on is an
index *entry*, and every column outside the index is another b-tree descent into the table.
Making them covering (`lens_name(lens_id, kind, tokens)`, `lens_token(kind, token, lens_id)`)
removed those descents. Always read the plan, and read it for the word COVERING.

@section log-map Optimising the geometry map: four attempts, no result

The CPU map is ~1.0× against lensfun and the loop does not vectorise. The obvious cause is
the model dispatch inside the per-pixel evaluator, which is loop-invariant but sits inside a
branch. Four ways of removing it were tried; the code is unchanged because none worked.

| attempt | result |
|---|---|
| specialise all twelve (distortion, TCA) pairs | **103 → 263 ms** — twelve call sites exhaust the inliner's budget, so nothing folds at all |
| specialise only the dominant pair (ptlens + poly3: 4810/5690 and 3355/3361 of the database) | no change |
| `restrict` on the evaluator arguments | ~5%, within run-to-run noise |
| make the row the primitive and the pixel `count == 1`, hoisting `y·scale − centre` and `1 − a − b − c` | no change — 100.7 ms against 101.1, best of five |

The premise was wrong. A hand-written loop with **every** branch resolved at compile time is
292 ms against the generic 286 under plain `-O3`: the branches were never the cost. What
remains is one square root, ~20 flops, and a 24-byte scatter into an interleaved output whose
layout is fixed by compatibility with lensfun's buffer. Closing it needs SIMD across pixels
with SoA staging — a different function, and one that would put the CPU/GPU bit-exactness at
risk for @ref fused-evaluation "the stage that a GPU does not build at all".

@section log-flags What the same experiment did turn up

The consumer's compiler flags dominate everything in this library. Same code, same 24 Mpx
map:

| | |
|---|---|
| `-O3` (baseline x86-64) | 286 ms |
| `-O3 -march=native` | 217 ms |
| `-O3 -march=native -ffast-math` | **103 ms** |

Ansel builds with the last, so its map costs ~100 ms and not the ~300 this project's own
benchmark reports with default flags. Any CPU figure quoted without its flags is close to
meaningless.

It also means a consumer built that way is **not** bit-identical to its own GPU path.
`-ffast-math` lets the compiler reassociate and contract; the parity harness is built without
it, and that is where the bit-exactness claim is asserted.

@section log-vignetting Vignetting: three causes, only one of them arithmetic

Vignetting was 1.8× *slower* than lensfun before it was measured properly. It is now 1.5×
faster, and only one of the three causes was about the maths:

1. **Three divides per pixel against lensfun's one.** The half-diagonal coordinate system was
   computed as `(xu · norm_scale − center) / aspect_ratio_correction`, per axis, per pixel.
   Folded into `vig_scale`/`vig_center_{x,y}` at resolve time.
2. **The row-invariant `y` term recomputed per pixel.** Hoisted.
3. **A `c == 0` guard that stopped the loop vectorising entirely** — the generated code had
   two `divss` and no `divps`. It was not upstream behaviour either: lensfun divides with no
   special case and clamps the product afterwards. Removing it leaves one divide and one
   select, and the loop vectorises across pixels: four divides per instruction.

Structuring the row as two passes — fill a block of multipliers, then scale the pixels — is
what lets that happen. A single fused loop gives the compiler a four-component store to
vectorise, and it takes that instead, leaving the divide scalar. On its own that change was
worth nothing; with the branch removed, the two together are the whole 2.7×.

@section log-projection Projection: two wrong theories before reading the source

The projection stage was gated off for a while behind a wrong diagnosis, and the fix was to
stop reasoning and clone lensfun.

The composition order was blamed first. The public header says the geometry callback has
priority 500 and distortion 750, so geometry-then-distortion — which was 28 px out at the
centre of the frame. The reverse order agreed to 0.03 px near the axis and was worse overall
(551k samples out of tolerance against 140k). Neither was right, and the conclusion drawn —
that the composition was not understood — was itself wrong.

Reading lensfun 0.3.4 settled it in minutes. The **order was right the first time**: the
header's "distortion 750" refers to `ModifyCoord_Dist_*`, the correction direction; the 250
`UnDist_*` variants are the reverse direction and are a Newton inversion. What was actually
wrong was the geometry **focal**: lensfun runs the callback on
`GetRealFocalLength(focal) / get_hugin_focal_correction(focal)`, and for a lens with no
`<real-focal-length>` data those cancel to exactly the nominal focal — measured ratio 1.0000
across every fisheye — while for the handful that carry it they do not. On the Sigma 4.5mm
circular fisheye the geometry focal is 0.4727× nominal, and using the nominal is 28 px out at
the centre and 135 px at the edge.

Two lessons, both cheap: **clone the sources** (`git clone --depth 1 --branch v0.3.4
https://github.com/lensfun/lensfun` — note HEAD is 0.3.99 and differs), and an installed
lensfun shipping only headers and a shared object is not a reason to guess.

@section log-benchmarking How to measure here

A shared machine only ever makes a measurement slower, never faster, so the **minimum** of
several runs is the closest thing to the cost of the code; the mean measures the machine's
other tenants. `tests/bench_lensfun.c` runs each stage five times and keeps the best.

That is not fussiness. With single samples, lensfun's single-threaded map read 333 ms in one
run and 532 ms in the next, and a thread sweep appeared to show lensfun scaling at 1.31×
against LensSerious at 3.44× — a 2.8× parallel advantage that does not exist. Best-of-five
puts both at ~3× and the ratios where they belong.

@section log-traps Traps that cost a round each

- **A blanket rename can rewrite its own target.** `sed s/_exif_read_exif_tag/dt_exif_read_exif_tag/`
  turned the *definition* into `dtdt_…`; it compiled, because the macro that called it was
  renamed consistently, and failed only at run time as an undefined symbol. The same shape
  turned `#define LS_ATAN(x) atan(x)` into a self-referential macro, which only the OpenCL
  syntax test caught.
- **SQLite loses negative zero.** It stores a REAL with no fractional part as an integer, so
  upstream's 40 terms of `-0.0` read back as `+0.0`. Inert here — every evaluator only
  multiplies or adds these, and the one divide is guarded by `c != 0.f`, false for both zeros
  — but `tests/parity_db.c` counts them rather than tolerating a range.
- **Guard on all three OpenCL macros.** OpenCL mandates `__OPENCL_VERSION__` for
  `clBuildProgram`, but clang's offline mode — the only way to syntax-check a kernel without
  a device — defines `__OPENCL_C_VERSION__` and not the other. Guarding on one alone makes
  the shared header take the host branch and `#include <math.h>` into a kernel.
- **Fit a constant across the range it will be used in.** The projection focal
  (`focal · crop · √(ar²+1) / 21.633307`) was first fitted on 3:2 lenses alone, where it
  reduces to a bare `focal/12`, and was **284 px out** on a 4:3 compact.
