# Workflow benchmarks

End-to-end benchmarks that reproduce real user reductions, with GBTIDL
equivalents where they exist. Ported from the benchmark harness developed in
PR #1049. These complement the micro-benchmarks in the parent `benchmark/`
directory: micro-benchmarks localize a cost, workflow benchmarks show whether
users will feel an improvement.

## Benchmarks

| Name            | Workflow                                                | GBTIDL | Verifier |
| --------------- | ------------------------------------------------------- | ------ | -------- |
| `hi_survey`     | gettp/getsigref survey reduction + smooth/baseline/stats | yes    | RMS of two line-free regions, atol 1 mK |
| `argus_vanecal` | Argus VANE/SKY Tsys calibration across 16 feeds          | yes    | mean Tsys per feed, atol 1 mK |
| `nod_kfpa`      | KFPA nodding data load                                   | no     | none |
| `exit`          | process-startup baseline (subtract from the others)      | yes    | none |


## Running

`run_bench.py` can be run from any directory (`uv run --project <repo> python
<path>/run_bench.py ...`). Script, golden and verifier paths are relative to `run_bench.py`,
and every child process runs with `benchmark/workflows/` as its working directory. dysh scripts
are launched with `uv run --frozen`, so `uv.lock` is used as is and never rewritten. From this
directory:

```bash
# time everything available on this host, 3 iterations, warm cache
uv run python run_bench.py

# one benchmark, cold cache, JSON output
uv run python run_bench.py --benchmarks hi_survey --mode cold --output results.json

# check correctness against GBTIDL (at GBO) or golden captures (elsewhere)
uv run python run_bench.py --verify
```

Data paths resolve in this order: `--data-path`, `$DYSH_BENCH_DATA_PATH`, the
canonical GBO path if it exists on this host, then the `dysh_data` alias
(which may download data on first use). Benchmarks whose data cannot be found
are skipped. Note for `argus_vanecal`: the canonical GBO dataset is a
vane/sky-only subset; the `otf4` alias fallback is the full session, which
inflates the `GBTFITSLoad` stage but leaves the vanecal stages comparable.

## Adding a benchmark

Add an entry to `BENCHMARKS` in `run_bench.py`. `dysh_script` (and `gbtidl_script`) may be a
path, or an argv list with the script first and its arguments after, and need not live under
`scripts/`:

```python
"getps": {"dysh_script": ["../bench_getps.py", "-t", "-l", "4"], "gbtidl_script": None, ...}
```

Only the first item is resolved (relative to `run_bench.py`); the rest are passed unchanged. A
script must print `DYSH_BENCH_SCRIPT_MS=<ms>` (and optionally `DYSH_BENCH_STAGE_MS[<stage>]=<ms>`
lines), which `bench_common.MarkerDTime` does. In `--mode cold` the data path may be a file or a
directory; the script must read it from `$DYSH_BENCH_DATA_PATH` (see `bench_common.resolve_data`).

`--verify` reuses the stdout of the first timed run instead of running the script again.

## Golden captures

Where GBTIDL is unavailable, `--verify` compares dysh's verification values
(`RMS_*`, `TSYS_FDNUM_*` stdout markers) against `scripts/<name>/golden.txt`,
using an absolute tolerance of 1 mK (see PR 0.5 in `perfenhance_plan.md`).
Regenerate goldens with:

```bash
uv run python run_bench.py --capture-golden
```

Capture goldens only from a known-good baseline (e.g., the release branch a
performance PR is based on), and only on the canonical dataset — golden values
are dataset-dependent. A performance PR must not change these values; if one
does, that is a correctness regression, not a benchmark update.
