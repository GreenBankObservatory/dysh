"""
Shared utilities for dysh benchmarks.
"""

import gc
import json
import logging
import os
import subprocess
import sys
import time
import tracemalloc
import warnings
from datetime import datetime
from statistics import mean, stdev

import psutil

# Configure logging
logger = logging.getLogger("dysh.benchmark")


def setup_logging(verbose=False):
    """Set up benchmark logging. Default INFO, -v for DEBUG."""
    level = logging.DEBUG if verbose else logging.INFO
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter("%(message)s"))
    logger.addHandler(handler)
    logger.setLevel(level)

    # Also configure dysh's logger to same level
    from dysh.log import init_logging

    verbosity = 3 if verbose else 2  # 3=DEBUG, 2=INFO
    init_logging(verbosity=verbosity)

    return logger


def get_process_memory_mb():
    """Get current process memory usage in MB (RSS)."""
    process = psutil.Process(os.getpid())
    return process.memory_info().rss / (1024 * 1024)


def get_git_info():
    """Get current git commit hash and branch."""
    try:
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"], stderr=subprocess.DEVNULL).decode().strip()
        branch = (
            subprocess.check_output(["git", "rev-parse", "--abbrev-ref", "HEAD"], stderr=subprocess.DEVNULL)
            .decode()
            .strip()
        )
        return {"commit": commit[:12], "branch": branch}
    except Exception:
        return {"commit": "unknown", "branch": "unknown"}


def time_operation(func, n_iterations=5, warmup=1, silent_errors=False, track_memory=True, name=None):
    """
    Time an operation multiple times and return statistics.

    Parameters
    ----------
    func : callable
        Function to time
    n_iterations : int
        Number of timed iterations
    warmup : int
        Number of warmup iterations (not timed)
    silent_errors : bool
        If True, don't raise errors during warmup
    track_memory : bool
        If True, track memory usage (default True)
    name : str, optional
        Name of operation for logging

    Returns
    -------
    dict
        Statistics including mean_ms, std_ms, min_ms, max_ms, n_iterations,
        and memory stats if track_memory=True
    """
    op_name = name or func.__name__

    # Warmup runs
    logger.debug(f"    [{op_name}] warming up ({warmup} iter)...")
    for i in range(warmup):
        gc.collect()
        try:
            func()
            logger.debug(f"    [{op_name}] warmup {i + 1}/{warmup} complete")
        except Exception as e:
            logger.debug(f"    [{op_name}] warmup {i + 1}/{warmup} error: {e}")
            if not silent_errors:
                raise

    # Timed runs
    times = []
    memory_deltas = []
    peak_memory_allocs = []

    logger.debug(f"    [{op_name}] timing ({n_iterations} iter)...")
    for i in range(n_iterations):
        gc.collect()

        if track_memory:
            mem_before = get_process_memory_mb()
            tracemalloc.start()

        start = time.perf_counter()
        result = func()
        elapsed = time.perf_counter() - start

        if track_memory:
            current, peak = tracemalloc.get_traced_memory()
            tracemalloc.stop()
            mem_after = get_process_memory_mb()
            memory_deltas.append(mem_after - mem_before)
            peak_memory_allocs.append(peak / (1024 * 1024))  # Convert to MB

        elapsed_ms = elapsed * 1000
        times.append(elapsed_ms)
        logger.debug(
            f"    [{op_name}] iter {i + 1}/{n_iterations}: {elapsed_ms:.1f} ms, peak={peak_memory_allocs[-1] if track_memory else 0:.1f} MB"
        )
        del result

    stats = {
        "mean_ms": round(mean(times), 2),
        "std_ms": round(stdev(times), 2) if len(times) > 1 else 0,
        "min_ms": round(min(times), 2),
        "max_ms": round(max(times), 2),
        "n_iterations": n_iterations,
    }

    if track_memory and memory_deltas:
        stats["memory"] = {
            "mean_delta_mb": round(mean(memory_deltas), 2),
            "peak_alloc_mb": round(max(peak_memory_allocs), 2),
        }

    mem_str = f", peak={stats['memory']['peak_alloc_mb']:.1f} MB" if track_memory and memory_deltas else ""
    logger.debug(f"    [{op_name}] result: {stats['mean_ms']:.1f} ms (±{stats['std_ms']:.1f}){mem_str}")
    return stats


def create_results_dict(quick_mode=False, extra_metadata=None):
    """
    Create a results dictionary with standard metadata.

    Parameters
    ----------
    quick_mode : bool
        Whether running in quick mode
    extra_metadata : dict, optional
        Additional metadata to include

    Returns
    -------
    dict
        Results dictionary with metadata section
    """
    results = {
        "metadata": {
            "timestamp": datetime.now().isoformat(),
            "git": get_git_info(),
            "python_version": sys.version.split()[0],
            "quick_mode": quick_mode,
        },
        "benchmarks": {},
    }
    if extra_metadata:
        results["metadata"].update(extra_metadata)
    return results


def save_results(results, output_path):
    """Save results to JSON file."""
    with open(output_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to: {output_path}")


def compare_results(file1, file2):
    """
    Compare two benchmark result files and print a summary.

    Parameters
    ----------
    file1 : str
        Path to baseline results JSON
    file2 : str
        Path to candidate results JSON
    """
    with open(file1) as f:
        results1 = json.load(f)
    with open(file2) as f:
        results2 = json.load(f)

    print("\n" + "=" * 70)
    print("BENCHMARK COMPARISON")
    print("=" * 70)

    print(f"\nBaseline: {file1}")
    print(f"  Git: {results1['metadata']['git']['branch']} @ {results1['metadata']['git']['commit']}")
    print(f"  Timestamp: {results1['metadata']['timestamp']}")

    print(f"\nCandidate: {file2}")
    print(f"  Git: {results2['metadata']['git']['branch']} @ {results2['metadata']['git']['commit']}")
    print(f"  Timestamp: {results2['metadata']['timestamp']}")

    print("\n" + "-" * 70)
    print(f"{'Benchmark':<35} {'Baseline':>12} {'Candidate':>12} {'Speedup':>12}")
    print("-" * 70)

    common = set(results1["benchmarks"].keys()) & set(results2["benchmarks"].keys())

    speedups = []
    for name in sorted(common):
        b1 = results1["benchmarks"][name]
        b2 = results2["benchmarks"][name]

        if "error" in b1 or "error" in b2:
            continue
        if "skipped" in b1 or "skipped" in b2:
            continue

        t1 = b1.get("mean_ms", 0)
        t2 = b2.get("mean_ms", 0)

        if t1 > 0 and t2 > 0:
            speedup = t1 / t2
            speedups.append(speedup)

            if speedup > 1.1:
                speedup_str = f"{speedup:.2f}x FASTER"
            elif speedup < 0.9:
                speedup_str = f"{speedup:.2f}x slower"
            else:
                speedup_str = f"{speedup:.2f}x"

            print(f"{name:<35} {t1:>10.1f}ms {t2:>10.1f}ms {speedup_str:>12}")

    print("-" * 70)
    if speedups:
        avg_speedup = mean(speedups)
        print(f"{'Average speedup:':<35} {' ':>12} {' ':>12} {avg_speedup:.2f}x")
    print("=" * 70)


def print_summary(results):
    """Print a summary of benchmark results."""
    print("\n" + "=" * 70)
    print("RESULTS SUMMARY")
    print("=" * 70)
    print(f"{'Benchmark':<30} {'Time':>12} {'Peak Mem':>12} {'Mem Delta':>12}")
    print("-" * 70)

    for name, data in results["benchmarks"].items():
        if "error" in data:
            print(f"{name:<30} ERROR - {data['error']}")
        elif "skipped" in data:
            print(f"{name:<30} SKIPPED - {data.get('reason', 'not supported')}")
        elif "mean_ms" in data:
            time_str = f"{data['mean_ms']:.1f} ms"
            mem = data.get("memory", {})
            peak_str = f"{mem.get('peak_alloc_mb', 0):.1f} MB" if mem else "n/a"
            delta_str = f"{mem.get('mean_delta_mb', 0):+.1f} MB" if mem else "n/a"
            print(f"{name:<30} {time_str:>12} {peak_str:>12} {delta_str:>12}")


def suppress_warnings():
    """Suppress common warnings during benchmarks."""
    warnings.filterwarnings("ignore", category=UserWarning)
