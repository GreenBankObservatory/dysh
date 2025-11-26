#!/usr/bin/env python
"""
Benchmark script for MapReduce-style operations in dysh.

This benchmarks the typical workflow for mapping/gridding data:
- GBTOffline initialization
- gettp for off scans
- timeaverage operations
- getsigref for signal scans
- ScanBlock operations

Usage:
    # Run benchmark with default project
    python benchmark/bench_mapreduce.py -o results.json

    # Specify a different project
    python benchmark/bench_mapreduce.py --project TGBT25B_608_07 -o results.json

    # Quick benchmark (fewer iterations)
    python benchmark/bench_mapreduce.py --quick -o results.json

    # Compare two result files
    python benchmark/bench_mapreduce.py --compare results_main.json results_branch.json
"""

import argparse
import sys
from pathlib import Path

import numpy as np

# Add benchmark directory to path for imports
sys.path.insert(0, str(Path(__file__).parent))

from bench_utils import (
    compare_results,
    create_results_dict,
    get_git_info,
    logger,
    print_summary,
    save_results,
    setup_logging,
    suppress_warnings,
    time_operation,
)


def benchmark_gbtoffline_init(project, n_iterations=3, use_index_file=True):
    """Benchmark GBTOffline initialization."""
    from dysh.fits.gbtfitsload import GBTOffline

    kwargs = {}
    if not use_index_file:
        kwargs["index_file_threshold"] = float('inf')

    def init_offline():
        sdf = GBTOffline(project, **kwargs)
        return sdf

    return time_operation(init_offline, n_iterations=n_iterations, warmup=1)


def benchmark_gettp(sdf, scan, ifnum, plnum, fdnum, n_iterations=5):
    """Benchmark gettp operation."""
    def do_gettp():
        return sdf.gettp(scan=scan, ifnum=ifnum, plnum=plnum, fdnum=fdnum)

    return time_operation(do_gettp, n_iterations=n_iterations, warmup=1)


def benchmark_timeaverage(scanblock, n_iterations=5):
    """Benchmark timeaverage operation on a ScanBlock."""
    def do_timeaverage():
        return scanblock.timeaverage()

    return time_operation(do_timeaverage, n_iterations=n_iterations, warmup=1)


def benchmark_spectrum_average(spec1, spec2, n_iterations=10):
    """Benchmark averaging two spectra."""
    def do_average():
        return spec1.average(spec2)

    return time_operation(do_average, n_iterations=n_iterations, warmup=1)


def benchmark_getsigref(sdf, scans, ref, ifnum, plnum, fdnum, n_iterations=3):
    """Benchmark getsigref operation."""
    def do_getsigref():
        return sdf.getsigref(scan=scans, ref=ref, ifnum=ifnum, plnum=plnum, fdnum=fdnum)

    return time_operation(do_getsigref, n_iterations=n_iterations, warmup=1)


def benchmark_scanblock_write(scanblock, output_path, n_iterations=3):
    """Benchmark ScanBlock write operation."""
    def do_write():
        scanblock.write(output_path)
        return True

    result = time_operation(do_write, n_iterations=n_iterations, warmup=1)

    # Clean up test files
    import os
    if os.path.exists(output_path):
        os.remove(output_path)

    return result


def benchmark_full_mapreduce_pipeline(sdf, off_scans, signal_scans, ifnum, plnum, fdnum, n_iterations=1):
    """
    Benchmark the full MapReduce pipeline for one beam/polarization.

    This is the end-to-end time for:
    1. gettp on off scans
    2. timeaverage each
    3. average the off spectra
    4. getsigref with signal scans
    """
    def do_pipeline():
        # Get off scans and average them
        off_specs = []
        for scan in off_scans:
            off_sb = sdf.gettp(scan=scan, ifnum=ifnum, plnum=plnum, fdnum=fdnum)
            off_specs.append(off_sb.timeaverage())

        # Average off spectra together
        ref = off_specs[0]
        for spec in off_specs[1:]:
            ref = ref.average(spec)

        # Get signal with reference
        result = sdf.getsigref(scan=signal_scans, ref=ref, ifnum=ifnum, plnum=plnum, fdnum=fdnum)
        return result

    return time_operation(do_pipeline, n_iterations=n_iterations, warmup=1)


def run_benchmarks(project, quick=False, use_index_file=True):
    """Run MapReduce benchmarks."""
    results = create_results_dict(
        quick_mode=quick,
        extra_metadata={"project": project, "use_index_file": use_index_file}
    )

    n_iter = 2 if quick else 5
    n_iter_slow = 1 if quick else 3

    suppress_warnings()

    # Import here to catch import errors gracefully
    try:
        from dysh.fits.gbtfitsload import GBTOffline
    except ImportError as e:
        print(f"ERROR: Cannot import dysh: {e}")
        sys.exit(1)

    kwargs = {}
    if not use_index_file:
        kwargs["index_file_threshold"] = float('inf')

    logger.info(f"\n=== GBTOffline Initialization ===")
    logger.info(f"  - Loading project {project}...")

    try:
        results["benchmarks"]["init"] = benchmark_gbtoffline_init(project, n_iterations=n_iter_slow, use_index_file=use_index_file)
    except Exception as e:
        logger.error(f"  ERROR: {e}")
        results["benchmarks"]["init"] = {"error": str(e)}
        return results

    # Load once for subsequent benchmarks
    sdf = GBTOffline(project, **kwargs)
    results["metadata"]["n_rows"] = len(sdf._index)
    results["metadata"]["n_files"] = len(sdf._sdf) if hasattr(sdf, '_sdf') else 1

    logger.info(f"  Loaded {results['metadata']['n_rows']} rows from {results['metadata']['n_files']} file(s)")

    # Get available scans to determine what we can benchmark
    available_scans = sorted(sdf._index["SCAN"].unique())
    logger.info(f"  Available scans: {available_scans[:10]}{'...' if len(available_scans) > 10 else ''}")

    # Use scans from the MapReduce.py example if available, otherwise adapt
    # Off scans: 19, 58
    # Signal scans: 20-53, 59-77
    off_scans = [19, 58]
    signal_scans = np.r_[20:54, 59:78]

    # Check if expected scans exist
    off_scans = [s for s in off_scans if s in available_scans]
    signal_scans = [s for s in signal_scans if s in available_scans]

    if not off_scans:
        # Fall back to first two scans
        off_scans = available_scans[:2] if len(available_scans) >= 2 else available_scans[:1]
    if not signal_scans:
        # Fall back to remaining scans
        signal_scans = [s for s in available_scans if s not in off_scans][:10]

    logger.info(f"  Using off scans: {off_scans}")
    logger.info(f"  Using signal scans: {signal_scans[:5]}{'...' if len(signal_scans) > 5 else ''} ({len(signal_scans)} total)")

    # Get available ifnum, plnum, fdnum
    ifnums = sorted(sdf._index["IFNUM"].unique())
    plnums = sorted(sdf._index["PLNUM"].unique())
    fdnums = sorted(sdf._index["FDNUM"].unique())

    ifnum = ifnums[0] if ifnums else 0
    plnum = plnums[0] if plnums else 0
    fdnum = fdnums[0] if fdnums else 0

    logger.info(f"  Using ifnum={ifnum}, plnum={plnum}, fdnum={fdnum}")

    # Benchmark individual operations
    logger.info(f"\n=== Individual Operation Benchmarks ===")

    if off_scans:
        logger.info(f"  - gettp (scan {off_scans[0]})...")
        try:
            results["benchmarks"]["gettp"] = benchmark_gettp(
                sdf, off_scans[0], ifnum, plnum, fdnum, n_iterations=n_iter
            )
        except Exception as e:
            results["benchmarks"]["gettp"] = {"error": str(e)}

        # Get a ScanBlock for timeaverage benchmark
        logger.info(f"  - timeaverage...")
        try:
            test_sb = sdf.gettp(scan=off_scans[0], ifnum=ifnum, plnum=plnum, fdnum=fdnum)
            results["benchmarks"]["timeaverage"] = benchmark_timeaverage(test_sb, n_iterations=n_iter)
        except Exception as e:
            results["benchmarks"]["timeaverage"] = {"error": str(e)}

        # Benchmark spectrum averaging
        if len(off_scans) >= 2:
            logger.info(f"  - spectrum average...")
            try:
                spec1 = sdf.gettp(scan=off_scans[0], ifnum=ifnum, plnum=plnum, fdnum=fdnum).timeaverage()
                spec2 = sdf.gettp(scan=off_scans[1], ifnum=ifnum, plnum=plnum, fdnum=fdnum).timeaverage()
                results["benchmarks"]["spectrum_average"] = benchmark_spectrum_average(
                    spec1, spec2, n_iterations=n_iter
                )
            except Exception as e:
                results["benchmarks"]["spectrum_average"] = {"error": str(e)}

    # Benchmark getsigref
    if off_scans and signal_scans:
        logger.info(f"  - getsigref ({len(signal_scans)} signal scans)...")
        try:
            # Create reference spectrum
            ref = sdf.gettp(scan=off_scans[0], ifnum=ifnum, plnum=plnum, fdnum=fdnum).timeaverage()
            if len(off_scans) > 1:
                ref2 = sdf.gettp(scan=off_scans[1], ifnum=ifnum, plnum=plnum, fdnum=fdnum).timeaverage()
                ref = ref.average(ref2)

            results["benchmarks"]["getsigref"] = benchmark_getsigref(
                sdf, signal_scans, ref, ifnum, plnum, fdnum, n_iterations=n_iter_slow
            )
            results["benchmarks"]["getsigref"]["n_signal_scans"] = len(signal_scans)
        except Exception as e:
            results["benchmarks"]["getsigref"] = {"error": str(e)}

    # Benchmark full pipeline
    if off_scans and signal_scans:
        logger.info(f"\n=== Full Pipeline Benchmark ===")
        logger.info(f"  - Full MapReduce pipeline (1 beam, 1 pol)...")
        try:
            results["benchmarks"]["full_pipeline"] = benchmark_full_mapreduce_pipeline(
                sdf, off_scans, signal_scans, ifnum, plnum, fdnum, n_iterations=n_iter_slow
            )
            results["benchmarks"]["full_pipeline"]["n_off_scans"] = len(off_scans)
            results["benchmarks"]["full_pipeline"]["n_signal_scans"] = len(signal_scans)
        except Exception as e:
            results["benchmarks"]["full_pipeline"] = {"error": str(e)}

    # Benchmark multiple beams if available (like the 7-beam KFPA example)
    if len(fdnums) > 1 and off_scans and signal_scans:
        logger.info(f"\n=== Multi-beam Pipeline Benchmark ({len(fdnums)} beams) ===")
        logger.info(f"  - Processing all {len(fdnums)} beams...")
        try:
            def do_all_beams():
                from dysh.spectra import ScanBlock
                all_results = ScanBlock()
                for fd in fdnums:
                    # Get off reference
                    offs = []
                    for scan in off_scans:
                        offs.append(sdf.gettp(scan=scan, ifnum=ifnum, plnum=plnum, fdnum=fd).timeaverage())
                    ref = offs[0]
                    for o in offs[1:]:
                        ref = ref.average(o)
                    # Get signal
                    sb = sdf.getsigref(scan=signal_scans, ref=ref, ifnum=ifnum, plnum=plnum, fdnum=fd)
                    all_results.append(sb)
                return all_results

            results["benchmarks"]["all_beams_pipeline"] = time_operation(
                do_all_beams, n_iterations=n_iter_slow, warmup=1
            )
            results["benchmarks"]["all_beams_pipeline"]["n_beams"] = len(fdnums)
        except Exception as e:
            results["benchmarks"]["all_beams_pipeline"] = {"error": str(e)}

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Benchmark MapReduce-style operations in dysh",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument(
        "-o", "--output",
        help="Output JSON file for results"
    )
    parser.add_argument(
        "--project",
        default="TGBT25B_608_07",
        help="GBT project name (default: TGBT25B_608_07)"
    )
    parser.add_argument(
        "--quick",
        action="store_true",
        help="Run quick benchmark with fewer iterations"
    )
    parser.add_argument(
        "--compare",
        nargs=2,
        metavar=("BASELINE", "CANDIDATE"),
        help="Compare two result files instead of running benchmarks"
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose logging"
    )
    parser.add_argument(
        "--no-index-file",
        action="store_true",
        help="Disable .index file usage (force reading from FITS)"
    )

    args = parser.parse_args()

    setup_logging(verbose=args.verbose)

    if args.compare:
        compare_results(args.compare[0], args.compare[1])
        return

    use_index_file = not args.no_index_file

    print("dysh MapReduce Benchmark")
    print("=" * 50)
    print(f"Git: {get_git_info()['branch']} @ {get_git_info()['commit']}")
    print(f"Project: {args.project}")
    print(f"Quick mode: {args.quick}")
    print(f"Use .index file: {use_index_file}")

    results = run_benchmarks(args.project, quick=args.quick, use_index_file=use_index_file)

    print_summary(results)

    if args.output:
        save_results(results, args.output)
    else:
        print("\nTip: Use -o results.json to save results for comparison")


if __name__ == "__main__":
    main()
