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

> **Note:** `verify.py` for `hi_survey`/`argus_vanecal` currently compares with `rtol=0.02` (2%
> relative). PR 0.5 in `perfenhance_plan.md` switches this to an absolute 1 mK tolerance before
> PR 1 begins — relative tolerance is unreliable near zero and too loose far from zero for a
> quantity whose meaningful noise floor is fixed in temperature units.

## Running

From this directory:

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
