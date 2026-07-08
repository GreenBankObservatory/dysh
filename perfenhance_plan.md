# Performance Enhancement Plan 

## Context

dysh is slower than GBTIDL for common operations (`GBTFITSLoad.get*`, `Spectrum` operations),
which blocks user migration from GBTIDL. PR #1049 ("Performance Enhancements", closed 2026-04-23,
53 commits, based on release-1.1) achieved 8–10x speedups on real workflows but was closed
unmerged because it grew unmanageable. This plan recovers the still-applicable optimizations from
that PR as a series of small, independently reviewable PRs, each with before/after benchmarks.

**Already merged into main (do not redo):** `.index` file caching (`index_file.py`, optional Rust
parser), fitsio selective row reads (`SDFITSLoad.rawspectra` backend selection), lazy metadata
loading core (`_load_full_rows_if_needed`), lazy flags (`lazyflag.py`). `timeaverage`, `scale`,
`average_spectra`, and `subtract_baseline` are already vectorized.

**Method:** reimplement on top of current main using the PR-1049 commit diffs as a guide
(`git diff <sha>^ <sha> -- <file>` against the locally available `origin/perf-enhancements`),
not cherry-pick — 14 src commits have landed since the PR's base, and the branch contains
interleaved fixes and reverts.

**Decisions made with Marc (2026-07-07):**
- Vendor the PR-1049 workflow benchmark drivers into `benchmark/` (PR 0)
- IERS: preload only; do NOT set `iers_conf.auto_download = False`
- Pure perf only: no API additions. `Spectrum.stats(brange/erange)` (commit f3d4ae97) and new
  `Scan/ScanBlock.baseline()` methods (469be8fc) are out of scope — file as separate issues
- fitsio handle reuse: persistent handles (one per file), not an LRU pool
  (rationale/trade-offs recorded in `fitsio_handles.md`)

---

## PR 0 — Benchmark harness and protocol (merge first)

**Scope:**
- Port PR-1049 `dysh-gbtidl-bench/` workflow drivers into `benchmark/workflows/`:
  `hi_survey`, `argus_vanecal` (and `nod_kfpa` if data available), including the `verify.py`
  result checkers. Retrieve via `git show origin/perf-enhancements:dysh-gbtidl-bench/...`.
  The checkers double as correctness gates: capture golden values (Tsys arrays,
  `stats(qac=True)` strings) once on main.
- Add/extend a `DTime`-tagged driver (extend `bench_getps.py` or add `bench_calibration.py`)
  timing: load, first `getps` (lazy-load penalty), repeated `getps` (warm), `timeaverage`,
  and a `getspec` loop (exposes WCS/SkyCoord/deepcopy per-spectrum costs).
- Document the protocol in `benchmark/README.md`.

**Protocol (used verbatim in every PR description):**
1. `uv sync --frozen`; record Python version and hostname; `OMP_NUM_THREADS=1`;
   data via `dysh_data` aliases (`DYSH_DATA=/bigdisk/data/gbt/dysh_data` locally).
2. Warm runs (most PRs): 5 runs, discard first, report median ± min/max of DTime tags.
3. Cold runs (I/O PRs 2 and 7 only): `sync; echo 1 | sudo tee /proc/sys/vm/drop_caches`
   before each of 3 runs; additionally one NFS-path run if available (handle reuse is
   disproportionately an NFS win).
4. One `python -m cProfile` run before/after per PR to confirm the targeted function left the
   top-10 cumulative list. Never mix profiled and unprofiled timings (1.5–2.5x overhead).
5. Standard datasets: `example="getps"` (AGBT05B_047_01, 45 MB, warm), `example="getpslarge"`
   (TGBT21A_501_11, 7.5 GB, cold/lazy-load), plus workflow data.
6. Gate: `uv run pytest -m "not gbo_only and not requires_internet"` and workflow verifiers
   pass before numbers are reported.

---

## PR sequence (ordered by value/risk)

### PR 1 — `Selection._lightweight_copy()` for calibration selection
- **Files:** `src/dysh/util/selection.py`; `src/dysh/fits/gbtfitsload.py:1969`
- **Change:** add `_lightweight_copy()` (shallow DataFrame copy via `DataFrame.__init__(result,
  self)`, `dict(self._selection_rules)`, shared `_aliases`/`_defkeys`, `self._table.copy()`);
  use it in `_common_selection`, replacing `copy.deepcopy(self._selection)` — the single site
  feeding gettp/getsigref/getps/getfs/getnod/subbeamnod.
- **Informed by:** 4c90917d (impl), db3a6e9a (tests). PR-1049 called this THE per-call bottleneck.
- **Risk:** medium-low. Tests: port the db3a6e9a isolation guarantees (mutate rules/channel
  selections on copy, assert original unchanged, both directions); full selection + flag suites.
- **Benchmarks:** `bench_getps.py` warm (both datasets), `hi_survey`.

### PR 2 — Persistent fitsio handles; cheap `nchan`
- **Files:** `src/dysh/fits/sdfitsload.py` (fresh opens at lines 145, 324, 672, 723, 936)
- **Change:** cache one open `fitsio.FITS` per `SDFITSLoad`; explicit `close()` + `__del__`
  fallback; invalidate on write/online-refresh paths; reopen on error (ESTALE);
  `__getstate__`/`__setstate__` drop handles for pickling. See `fitsio_handles.md`.
- **Informed by:** e7319925.
- **Risk:** low-medium. Tests: loader suite; multi-file session load + re-read after `close()`;
  write-then-read invalidation.
- **Benchmarks:** `bench_gbtfitsload.py` + `bench_getps.py` cold (drop_caches ×3); NFS run.

### PR 3 — Vectorize calibration inner loops
- **Files:** `src/dysh/spectra/core.py` (add `mean_tsys_vectorized` beside `mean_tsys`);
  `src/dysh/spectra/scan.py` — `TPScan._calc_tsys` (~1787), `PSScan.calibrate` (three loops at
  2143/2157/2169), `NodScan` (2494/2506), `FSScan` (2920/2978), `SubBeamNodScan` (3207).
- **Change:** vectorize the common path (`smoothref == 1`, no vane) as
  `tsys[:, np.newaxis] * (sig - ref) / ref`; keep per-i loop as fallback for smoothref>1/vane
  (exactly PR-1049's shape).
- **Informed by:** 5f5162cb (measured hi_survey warm 6.44→5.80 s).
- **Risk:** medium (masked-array semantics, NaN-Tsys rows, dtype drift). Tests: assert
  vectorized ≡ stacked `mean_tsys` incl. NaN rows; existing GBTIDL-golden calibration tests;
  PR-0 verifiers.
- **Benchmarks:** `hi_survey`, `bench_getps.py` warm, `bench_otf.py`.

### PR 4 — Scan metadata fast path + `timeaverage` deepcopy removal
- **Files:** `src/dysh/spectra/scan.py` (`_make_meta`:894, `_set_all_meta`:887,
  `_add_calibration_meta`:931, `_update_scale_meta`:946, `timeaverage`:980);
  optionally `src/dysh/spectra/spectrum.py:1644`.
- **Change:** compute derived columns on the DataFrame before `.to_dict("records")` instead of
  per-row loops; NumPy indexing for meta arrays; replace `timeaverage`'s
  `deepcopy(self.getspec(0, use_wcs=use_wcs))` with direct `make_spectrum` on averaged data.
  Optionally a private `make_spectrum` fast path (skip `deepcopy(meta)`) for Scan callers that
  hand over a throwaway dict — keep the public API's deepcopy isolation.
- **Informed by:** b2e20d4a, 3bad2cab, 248fe0dd, d1718e1d, 4c90917d item 6.
- **Risk:** medium (meta aliasing — shared-dict bug; averaged-meta field drift). Tests:
  meta-isolation (mutate `scan._meta[0]`, assert `_meta[1]` unchanged; mutate returned
  `Spectrum.meta`, assert scan meta unchanged); compare timeaverage data+meta against main.
- **Benchmarks:** `bench_getps.py -t`, `hi_survey`.

### PR 5 — WCS template cache in `make_spectrum`
- **Files:** `src/dysh/spectra/spectrum.py:1656–1697`; small `scan.py` companions.
- **Change:** module-level `@lru_cache` WCS builder keyed on the invariant subset of
  `wcs_meta_keys` (CTYPE1–4, CUNIT1–3, CRPIX1 — match current key set exactly); per call,
  deepcopy the cached WCS and patch CRVAL1–4/DATE-OBS/obsgeo (call `wcs.wcs.set()` after
  patching). Include the two companions: skip WCS for getsigref TP reference spectra
  (a5421566) and skip redundant single-scan averages without WCS (41be4f1c).
- **Informed by:** 4c90917d item 5.
- **Risk:** medium (cache key must cover every non-patched varying value; unhashable keys).
  Tests: cached-and-patched vs fresh `WCS(header=...)` equivalence (`to_header()` equality +
  `pixel_to_world` spot checks) across test datasets; multi-IF spectral_axis regression.
- **Benchmarks:** getspec-loop micro-benchmark, `bench_getps.py -t`, `hi_survey`.

### PR 6 — Cache bintable index slices in calibration entry points (optional)
- **Files:** `src/dysh/fits/gbtfitsload.py` (per-call `select_from("FITSINDEX", ...)` /
  `_common_scan_list_selection` in getps at 2676–2701 and siblings).
- **Change:** cache keyed on (selection state, scans, proc key); invalidate on any
  `select*`/`flag*`/lazy-load mutation.
- **Informed by:** ac889ebd.
- **Risk:** medium — invalidation is the whole game. Tests: staleness (select → getps →
  select more → getps, row counts change; same for flags and lazy-load column arrival).
- **Benchmarks:** `hi_survey`, `bench_getps.py -l 10`.
- **Note:** worst value/risk of the top group; drop if the ac889ebd invalidation surface looks
  too broad after review. Do not develop concurrently with PR 7 (same call path).

### PR 7 — Lazy-load pipeline efficiency
- **Files:** `src/dysh/fits/sdfitsload.py:695–753` (`load_full_rows`);
  `src/dysh/fits/gbtfitsload.py:1668–1821` (`_load_full_rows_if_needed`,
  `_rebuild_merged_index`).
- **Change:** build the lazy-load DataFrame without the per-column decode loop (vectorized
  decode, single construction); batch the per-sdf × per-column full-length `pd.Series`
  allocations into one aligned frame write; track fully-loaded columns so repeated get* calls
  stop re-triggering load+rebuild (4c90917d item 3); reuse Selection storage in the
  single-file `_rebuild_merged_index` case instead of reconstructing `Selection(df)`+`Flag(df)`.
- **Informed by:** 376f8191, 820838dc, 54e487a9.
- **Risk:** medium-high — touches the `_index_source` hybrid-mode state machine and rule
  preservation. Tests: lazy-load/hybrid suites; selection- and flag-rule survival across a lazy
  load; `test_lazy_flags_bigfile.py`; getps identical with `.index` present vs absent.
- **Benchmarks:** `bench_partial_loading.py`, first-`getps` tag on getpslarge, cold runs.

### PR 8 — Startup hygiene: IERS preload, deferred GUI/CLI imports
- **Files:** `src/dysh/__init__.py`, `src/dysh/plot/plotbase.py:22–61`, `src/dysh/cli.py`.
- **Change:** preload the IERS table once at import (`IERS_B.open()`); do NOT change
  `iers_conf.auto_download` (per Marc — science-policy neutral). Move plot GUI backend
  selection + tk/ShellGUI warnings from import time to first-plot time (768f2f27); defer CLI
  imports until commands run (4cca0d99).
- **Risk:** low. Tests: plot suite headless/notebook modes; `python -c "import dysh"` smoke.
- **Benchmarks:** `python -X importtime`, first-`getspec` latency.

### PR 9 — Target/observer coordinate reuse in scan averages (HIGH RISK — last)
- **Files:** `src/dysh/spectra/scan.py`, `src/dysh/spectra/spectrum.py` (make_spectrum
  target/observer params), `src/dysh/coordinates/core.py`.
- **Change:** narrowly the surviving PR-1049 approach — compute `make_target(_meta)` once per
  scan average when CRVAL2/3/RADESYS/EQUINOX/DATE-OBS are identical across combined rows;
  reuse observer ITRS only when obstime is identical. Do NOT reintroduce
  `_precompute_observer`/`_precompute_target` on Spectrum (reverted in f8048549 — broke
  velocity frames) or `lru_cache` on `make_target` (reverted a5850bf0/09699b55).
- **Informed by:** 010c9401, ae5cd844, 45a6f589; revert history f8048549.
- **Risk:** high — this exact area was reverted once. **Write the velocity-frame regression
  suite against main BEFORE making the change:** for each frame (topo/LSRK/bary/helio…)
  assert frame-converted spectral axes identical to float tolerance, on data with
  per-integration DATE-OBS variation; heterogeneous-target averages must not take the reuse
  path.
- **Benchmarks:** getspec-loop micro, `bench_getps.py -t`, `nod_kfpa` workflow if ported.

---

## Grouping rationale

- PR 1 vs PR 6 both attack `_common_selection` cost but have independent failure modes
  (aliasing vs staleness); PR 1 is near-free, PR 6 needs invalidation machinery.
- PR 3 (numeric golden values) vs PR 4 (meta isolation) have different correctness gates;
  either can be reverted alone.
- WCS cache + skip-WCS companions share one regression test, so they combine (PR 5).
- The three lazy-load commits all mutate the same load→rebuild cycle, so they combine (PR 7).
- PR 9 is last because of its revert history; its regression tests must exist first.

## Verification (every PR)

1. `uv run pytest -m "not gbo_only and not requires_internet"` (parallel OK).
2. Workflow verifiers (`benchmark/workflows/*/verify.py`) — results byte/tolerance-identical
   to golden values captured on main.
3. Before/after benchmark table per the PR-0 protocol in the PR description.
4. `uv run ruff check .` and `uv run ruff format --check .`.

## Out of scope (file as issues)

- `Spectrum.stats(brange/erange)` region args (f3d4ae97) — API addition.
- Vectorized `Scan.baseline()`/`ScanBlock.baseline()` (469be8fc) — new feature; these methods
  do not exist on main and `subtract_baseline` is already vectorized.
