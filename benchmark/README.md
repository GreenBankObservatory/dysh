# Notes for benchmarking Dysh

The purpose of this benchmarking is to find and address the various
CPU and I/O bottlenecks in dysh.

- `run_bench.py` (this directory): the runner for the workflow and micro-benchmarks below,
  and the before/after protocol.
- `workflows/`: end-to-end user-workflow benchmarks (GBTIDL comparison, verifiers, goldens);
  see `workflows/README.md`.
- `development/`: earlier, ad hoc `DTime`/ECSV benchmarks from the Q8-Q12 development phase
  (`bench_getps.py`, `bench_calibration.py`, `bench_gbtfitsload.py`, and others); `getps` and
  `calibration` are also runnable through `run_bench.py`. See `development/README.md`.

## Before/after protocol for performance PRs

Every performance PR must include a before/after table in its description,
produced as follows. "Before" is the PR's base branch, "after" is the PR head,
both benchmarked the same way on the same machine.

1. **Pin the environment.** `uv sync --frozen`; record the Python version and
   hostname alongside the numbers. Set `OMP_NUM_THREADS=1` (see below). Data
   comes from the `dysh_data` aliases with `DYSH_DATA` set.
   **Pin the `.index` state too.** Whether an `.index` file sits next to the
   data changes which code path dysh takes: with one, metadata loads lazily
   (fast `GBTFITSLoad`, but the first `get*` call pays the load cost); without
   one, dysh does a full load up front. On the 45 MB `getps` example
   `GBTFITSLoad` took ~0.10-0.16 s and the first `getps` ~1.5-1.9 s with the
   index, versus ~0.25 s and no first-call penalty without it. So the "before"
   and "after" runs must have the same index state, or the comparison measures
   the index, not the change. Before each series, check for the `.index` file
   beside the data (e.g. `ls AGBT05B_047_01.raw.acs.index`) and record whether
   it is present. dysh may also write an `.index` on load (see
   `index_file_threshold`), so re-check between the "before" and "after" runs.
   `run_bench.py --mode cold` copies a data file together with its
   same-stem siblings (`.index`, `.flag`, ...) so that cold and warm runs take
   the same path. GBTIDL comparisons should be run with no index file (see
   Notes below).
2. **Warm runs** (default, for CPU-bound changes): run the relevant benchmark
   5 times, discard the first (import and lazy-load warm-up), report the
   median and min/max of the timing tags. `run_bench.py` prints
   these as the `median` and `min-max` columns of its stage tables and stores
   them in its JSON output; `--iterations 4` (it runs one untimed warm-up
   first) is equivalent to 5 runs with the first discarded.
3. **Cold runs** (only for I/O-bound changes, e.g. FITS reading): drop the
   page cache before each of 3 runs (`sync; echo 1 | sudo tee
   /proc/sys/vm/drop_caches`, or use `run_bench.py --mode cold`
   which evicts with `posix_fadvise`); report all 3. If NFS-mounted data is
   available, include one NFS run — I/O optimizations are disproportionately
   an NFS win.
4. **Confirm the hotspot moved.** One `python -m cProfile` run before and
   after, checking that the targeted function left the top of the cumulative
   list. Never mix profiled and unprofiled timings (profiling overhead is
   1.5–2.5x).
5. **Standard datasets:** `example="getps"` (45 MB, warm CPU benchmarks),
   `example="getpslarge"` (7.5 GB, cold I/O and lazy-load benchmarks), plus
   the workflow datasets (see `workflows/README.md`).
6. **Regression gate:** `uv run pytest -m "not gbo_only and not
   requires_internet"` passes, and `run_bench.py --verify` passes,
   before any numbers are reported.

The relevant drivers: `run_bench.py` (end-to-end user workflows, GBTIDL
comparison, and the `getps`/`calibration` micro-benchmarks below), plus the
rest of `development/` (per-spectrum overhead, SDFITS loading, and so on;
see `development/README.md`).
