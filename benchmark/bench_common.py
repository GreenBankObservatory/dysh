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

from bench_stats import summarize
from dysh.util.timers import DTime

__all__ = ["DATA_PATH_ENV", "MarkerDTime", "add_dtime_args", "resolve_data", "summarize"]

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
        """Record an aggregate stage that overlaps several tags, e.g. a loop total.

        The stage runs from the most recent tag named `since` (exclusive) to the most recent tag
        made so far, so call it after the last tag that belongs to the span. The marker is printed by
        `report`, after the per-tag markers, as ``DYSH_BENCH_STAGE_MS[<name>]=<ms>``. Does nothing if
        the timer is inactive.

        Parameters
        ----------
        name : str
            Stage name for the marker.
        since : str
            Name of the tag that marks the start of the span. The span begins when that tag was made.

        Raises
        ------
        ValueError
            If no tag named `since` has been made.
        """
        if not self.active:
            return
        start = next((row[1] for row in reversed(self.stats) if row[0] == since), None)
        if start is None:
            raise ValueError(f"no tag named {since!r}")
        self._spans = [*getattr(self, "_spans", []), (name, (self.stats[-1][1] - start) / 1e6)]

    def report(self, debug=False):
        """Print the timing table (or write it to file), then the benchmark markers.

        The table, ECSV output, and profiler report are those of `DTime.report`. Afterwards one
        ``DYSH_BENCH_STAGE_MS[<tag>]=<ms>`` line is printed per tag, giving the time since the previous
        tag (the first tag is timed from when the timer was created), followed by any `span` stages and
        ``DYSH_BENCH_SCRIPT_MS=<ms>``, the time since the timer was created. That excludes interpreter
        startup and imports that happen before the timer is created. Repeated tag names are printed as
        ``name``, ``name#2``, ``name#3``, ...

        Parameters
        ----------
        debug : bool, optional
            Passed to `DTime.report`, which prints every raw entry.

        Raises
        ------
        Exception
            From `DTime.report`, if the ``out`` file exists and neither ``append`` nor ``overwrite``
            was requested.
        """
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


def add_dtime_args(parser, memory=True):
    """Add the options every `DTime`-based benchmark shares.

    `DTime` reads ``out``, ``append``, ``overwrite``, ``profile``, ``statslines`` and ``sortkey``
    from ``vars(args)`` by bare key, so a driver missing any of them fails with a `KeyError`.
    Adding them through this function guarantees they exist. The options are ``--out/-o``,
    ``--append/-a``, ``--overwrite/-w``, ``--profile/-p``, ``--statslines/-e`` (default 25),
    ``--sortkey/-x`` (default ``"cumulative"``) and ``--memory/-m``. The caller must not reuse
    those short flags for other options.

    Parameters
    ----------
    parser : `argparse.ArgumentParser`
        Parser to add the options to.
    memory : bool, optional
        Add ``--memory/-m``. Pass `False` for a driver that uses ``-m`` for something else.
        ``--memory`` is accepted but not currently read by any driver. Default is `True`.

    Returns
    -------
    `argparse.ArgumentParser`
        The same `parser`, for chaining.
    """
    # fmt: off
    parser.add_argument("--out",        "-o", action="store",      help="output filename (astropy Table)", required=False)
    parser.add_argument("--append",     "-a", action="store_true", help="append to previous output file (astropy Table)")
    parser.add_argument("--overwrite",  "-w", action="store_true", help="overwrite a previous output file (astropy Table)")
    parser.add_argument("--profile",    "-p", action="store_true", help="run the profiler")
    parser.add_argument("--statslines", "-e", action="store",      help="number of profiler statistics lines to print", default=25)
    parser.add_argument("--sortkey",    "-x", action="store",      help="How to sort the profiler statistics, 'cumulative' or 'time'", default="cumulative")
    if memory:
        parser.add_argument("--memory", "-m", action="store_true", help="track memory usage")
    # fmt: on
    return parser


def resolve_data(canonical=None, **dysh_data_kwargs):
    """Return the benchmark data path.

    Resolution order: ``$DYSH_BENCH_DATA_PATH`` (how ``run_bench.py`` points a script at a
    cold-cache copy), then `canonical` if it exists on this host, then `dysh.util.files.dysh_data`
    called with the remaining keyword arguments, e.g. ``resolve_data(example="getps")``.

    Parameters
    ----------
    canonical : str, optional
        Preferred path on hosts where it exists, e.g. the GBO location of the dataset.
    **dysh_data_kwargs
        Keyword arguments for `dysh.util.files.dysh_data`, used when neither of the above applies.

    Returns
    -------
    str or `~pathlib.Path` or None
        A `str` if taken from the environment or `canonical`; otherwise whatever `dysh_data`
        returns (a `~pathlib.Path`, or `None` if the data could not be found).
    """
    env_path = os.environ.get(DATA_PATH_ENV)
    if env_path:
        return env_path
    if canonical is not None and os.path.exists(canonical):
        return canonical
    from dysh.util.files import dysh_data

    return dysh_data(**dysh_data_kwargs)
