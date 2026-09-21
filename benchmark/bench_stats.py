"""
Summary statistics for benchmark timings.

Kept free of dysh imports so that ``workflows/run_bench.py`` (which does not import dysh) and the
micro-benchmarks can share one implementation.
"""

from statistics import mean, median, stdev

__all__ = ["summarize"]


def summarize(times):
    """Summarize a list of timings.

    Parameters
    ----------
    times : list of float
        Timings, in any single unit. Must not be empty.

    Returns
    -------
    dict
        ``n``, ``mean``, ``std`` (sample standard deviation, 0 for a single value), ``median``,
        ``min``, ``max``, and a copy of the raw ``times``, all in the units of the input.

    Raises
    ------
    statistics.StatisticsError
        If `times` is empty.
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
