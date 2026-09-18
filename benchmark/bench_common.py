"""
Helpers shared by the dysh micro-benchmarks and the workflow benchmarks.

The stdout marker protocol is the common contract between the two systems (and with the
GBTIDL scripts in ``workflows/scripts``)::

    DYSH_BENCH_STAGE_MS[<stage>]=<milliseconds>
    DYSH_BENCH_SCRIPT_MS=<milliseconds>

``workflows/run_bench.py`` parses these markers, so any script that uses `MarkerDTime` can be
run under it. This module deliberately imports nothing from a dysh version newer than
`dysh.util.timers`, so the benchmarks run unchanged against both the "before" and "after"
branches of a performance PR.
"""

import os
import time
from statistics import mean, median, stdev

from dysh.util.timers import DTime

DATA_PATH_ENV = "DYSH_BENCH_DATA_PATH"


class MarkerDTime(DTime):
    """A `DTime` that also prints ``DYSH_BENCH_*`` markers when `report` is called.

    The table, ECSV output, and profiler behavior are exactly those of `DTime`; the markers
    are purely additive. Repeated tag names are made unique (``name``, ``name#2``, ...) because
    ``run_bench.py`` keeps one value per stage name and drops stages that are not present in
    every run.

    Use `span` for an aggregate marker that overlaps several tags (e.g. a loop total).
    """

    def span(self, name, since):
        """Record an aggregate stage from the tag ``since`` (exclusive) to the latest tag.

        The marker is printed by `report` after the per-tag markers.
        """
        if not self.active:
            return
        start = next((row[1] for row in reversed(self.stats) if row[0] == since), None)
        if start is None:
            raise ValueError(f"no tag named {since!r}")
        self._spans = [*getattr(self, "_spans", []), (name, (self.stats[-1][1] - start) / 1e6)]

    def report(self, debug=False):
        super().report(debug=debug)
        if not self.active:
            return
        counts = {}
        for prev, cur in zip(self.stats, self.stats[1:], strict=False):
            name = cur[0]
            counts[name] = counts.get(name, 0) + 1
            label = name if counts[name] == 1 else f"{name}#{counts[name]}"
            print(f"DYSH_BENCH_STAGE_MS[{label}]={(cur[1] - prev[1]) / 1e6:.3f}")
        for name, ms in getattr(self, "_spans", []):
            print(f"DYSH_BENCH_STAGE_MS[{name}]={ms:.3f}")
        print(f"DYSH_BENCH_SCRIPT_MS={(time.perf_counter_ns() - self.stats[0][1]) / 1e6:.3f}")


def add_dtime_args(parser):
    """Add the options every `DTime`-based benchmark shares.

    `DTime` reads ``out``, ``append``, ``overwrite``, ``profile``, ``statslines`` and ``sortkey``
    from ``vars(args)`` by bare key, so a driver missing any of them fails with a `KeyError`.
    Adding them through this function guarantees they exist.
    """
    # fmt: off
    parser.add_argument("--out",        "-o", action="store",      help="output filename (astropy Table)", required=False)
    parser.add_argument("--append",     "-a", action="store_true", help="append to previous output file (astropy Table)")
    parser.add_argument("--overwrite",  "-w", action="store_true", help="overwrite a previous output file (astropy Table)")
    parser.add_argument("--profile",    "-p", action="store_true", help="run the profiler")
    parser.add_argument("--statslines", "-e", action="store",      help="number of profiler statistics lines to print", default=25)
    parser.add_argument("--sortkey",    "-x", action="store",      help="How to sort the profiler statistics, 'cumulative' or 'time'", default="cumulative")
    parser.add_argument("--memory",     "-m", action="store_true", help="track memory usage")
    # fmt: on
    return parser


def resolve_data(canonical=None, **dysh_data_kwargs):
    """Return the benchmark data path.

    Resolution order: ``$DYSH_BENCH_DATA_PATH`` (how ``run_bench.py`` points a script at a
    cold-cache copy), then ``canonical`` if it exists on this host, then
    `dysh.util.files.dysh_data` called with the remaining arguments, e.g.
    ``resolve_data(example="getps")``.
    """
    env_path = os.environ.get(DATA_PATH_ENV)
    if env_path:
        return env_path
    if canonical is not None and os.path.exists(canonical):
        return canonical
    from dysh.util.files import dysh_data

    return dysh_data(**dysh_data_kwargs)


def summarize(times):
    """Summary statistics of a list of timings, in the units of the input.

    Returns
    -------
    dict
        ``n``, ``mean``, ``std``, ``median``, ``min``, ``max`` and the raw ``times``.
    """
    return {
        "n": len(times),
        "mean": mean(times),
        "std": stdev(times) if len(times) > 1 else 0.0,
        "median": median(times),
        "min": min(times),
        "max": max(times),
        "times": list(times),
    }
