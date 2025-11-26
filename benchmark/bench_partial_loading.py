#!/usr/bin/env python
"""
Benchmark script for dysh operations.

This script measures the performance of key operations like getps, gettp, and
getfs. It works on both main and feature branches - run it in different checkouts
and compare the JSON outputs to measure speedups.

Usage:
    # Run benchmark and save results
    python benchmark/bench_partial_loading.py -o results.json

    # Compare two result files
    python benchmark/bench_partial_loading.py --compare results_main.json results_branch.json

    # Quick benchmark (fewer iterations)
    python benchmark/bench_partial_loading.py --quick -o results.json
"""

import argparse
import inspect
import sys
from pathlib import Path

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


def find_testdata_root():
    """Find the testdata directory."""
    # Try relative to script location
    script_dir = Path(__file__).parent.parent
    testdata = script_dir / "testdata"
    if testdata.exists():
        return testdata

    # Try current working directory
    testdata = Path.cwd() / "testdata"
    if testdata.exists():
        return testdata

    return None


def get_test_datasets():
    """Get available test datasets with proper operation mappings."""
    testdata = find_testdata_root()
    if testdata is None:
        return {}

    datasets = {}

    # Position Switch dataset (from positionswitch.ipynb)
    # AGBT05B_047_01 - scans 51, 53, 55, 57
    ps_file = testdata / "AGBT05B_047_01/AGBT05B_047_01.raw.acs/AGBT05B_047_01.raw.acs.fits"
    if ps_file.exists():
        datasets["getps"] = {
            "path": str(ps_file),
            "size_mb": round(ps_file.stat().st_size / (1024 * 1024), 1),
            "scans": [51, 53, 55, 57],
            "ifnum": 0,
            "plnum": 0,
            "fdnum": 0,
        }

    # Large file for rawspectra benchmarks
    large_file = testdata / "TGBT25B_603_12/TGBT25B_603_12.raw.vegas/TGBT25B_603_12.raw.vegas.A.fits"
    if large_file.exists():
        datasets["large"] = {
            "path": str(large_file),
            "size_mb": round(large_file.stat().st_size / (1024 * 1024), 1),
        }

    # Medium file
    medium_file = testdata / "AGBT05B_047_01/AGBT05B_047_01.raw.acs/AGBT05B_047_01.raw.acs.fits"
    if medium_file.exists():
        datasets["medium"] = {
            "path": str(medium_file),
            "size_mb": round(medium_file.stat().st_size / (1024 * 1024), 1),
        }

    # Nodding dataset (from nodding.ipynb) - check for testtrim version
    nod_file = testdata / "TGBT22A_503_02/TGBT22A_503_02.raw.vegas/TGBT22A_503_02.raw.vegas.testtrim.fits"
    if nod_file.exists():
        datasets["getnod"] = {
            "path": str(nod_file),
            "size_mb": round(nod_file.stat().st_size / (1024 * 1024), 1),
            "scan": 62,
            "ifnum": 0,
            "plnum": 0,
        }

    # Frequency switch - TGBT21A_504_01 has frequency switched data
    fs_file = testdata / "TGBT21A_504_01/TGBT21A_504_01.raw.vegas/TGBT21A_504_01.raw.vegas.A.fits"
    if fs_file.exists():
        datasets["getfs"] = {
            "path": str(fs_file),
            "size_mb": round(fs_file.stat().st_size / (1024 * 1024), 1),
            "scan": 20,
            "ifnum": 0,
            "plnum": 0,
            "fdnum": 0,
        }

    return datasets


def benchmark_gbtfitsload_init(filepath, n_iterations=3, use_index_file=True):
    """Benchmark GBTFITSLoad initialization."""
    from dysh.fits.gbtfitsload import GBTFITSLoad

    kwargs = {}
    if not use_index_file:
        kwargs["index_file_threshold"] = float('inf')

    def load_file():
        sdf = GBTFITSLoad(filepath, **kwargs)
        return sdf

    return time_operation(load_file, n_iterations=n_iterations, warmup=1)


def benchmark_getps(dataset, n_iterations=3, use_index_file=True):
    """Benchmark getps operation (position-switched calibration)."""
    from dysh.fits.gbtfitsload import GBTFITSLoad

    filepath = dataset["path"]
    scans = dataset.get("scans", [51])
    ifnum = dataset.get("ifnum", 0)
    plnum = dataset.get("plnum", 0)
    fdnum = dataset.get("fdnum", 0)

    kwargs = {}
    if not use_index_file:
        kwargs["index_file_threshold"] = float('inf')

    sdf = GBTFITSLoad(filepath, **kwargs)

    def do_getps():
        return sdf.getps(scan=scans, fdnum=fdnum, ifnum=ifnum, plnum=plnum, calibrate=True)

    result = time_operation(do_getps, n_iterations=n_iterations, warmup=1)
    result["scans"] = scans
    result["file_size_mb"] = dataset["size_mb"]

    del sdf
    return result


def benchmark_gettp(dataset, n_iterations=5, use_index_file=True):
    """Benchmark gettp operation (total power)."""
    from dysh.fits.gbtfitsload import GBTFITSLoad

    filepath = dataset["path"]
    scans = dataset.get("scans", [51])
    scan = scans[0] if isinstance(scans, list) else scans
    ifnum = dataset.get("ifnum", 0)
    plnum = dataset.get("plnum", 0)
    fdnum = dataset.get("fdnum", 0)

    kwargs = {}
    if not use_index_file:
        kwargs["index_file_threshold"] = float('inf')

    sdf = GBTFITSLoad(filepath, **kwargs)

    def do_gettp():
        return sdf.gettp(scan=scan, fdnum=fdnum, ifnum=ifnum, plnum=plnum, calibrate=True)

    result = time_operation(do_gettp, n_iterations=n_iterations, warmup=1)
    result["scan"] = scan

    del sdf
    return result


def benchmark_getfs(dataset, n_iterations=3, use_index_file=True):
    """Benchmark getfs operation (frequency-switched calibration)."""
    from dysh.fits.gbtfitsload import GBTFITSLoad

    filepath = dataset["path"]
    scan = dataset.get("scan", 20)
    ifnum = dataset.get("ifnum", 0)
    plnum = dataset.get("plnum", 0)
    fdnum = dataset.get("fdnum", 0)

    kwargs = {}
    if not use_index_file:
        kwargs["index_file_threshold"] = float('inf')

    sdf = GBTFITSLoad(filepath, **kwargs)

    def do_getfs():
        return sdf.getfs(scan=scan, fdnum=fdnum, ifnum=ifnum, plnum=plnum, calibrate=True)

    try:
        result = time_operation(do_getfs, n_iterations=n_iterations, warmup=1)
        result["scan"] = scan
        result["file_size_mb"] = dataset["size_mb"]
    except Exception as e:
        result = {"error": str(e)}

    del sdf
    return result


def benchmark_getnod(dataset, n_iterations=3, use_index_file=True):
    """Benchmark getnod operation (nodding calibration)."""
    from dysh.fits.gbtfitsload import GBTFITSLoad

    filepath = dataset["path"]
    scan = dataset.get("scan", 62)
    ifnum = dataset.get("ifnum", 0)
    plnum = dataset.get("plnum", 0)

    kwargs = {}
    if not use_index_file:
        kwargs["index_file_threshold"] = float('inf')

    sdf = GBTFITSLoad(filepath, **kwargs)

    def do_getnod():
        return sdf.getnod(scan=scan, ifnum=ifnum, plnum=plnum)

    try:
        result = time_operation(do_getnod, n_iterations=n_iterations, warmup=1)
        result["scan"] = scan
        result["file_size_mb"] = dataset["size_mb"]
    except Exception as e:
        result = {"error": str(e)}

    del sdf
    return result


def benchmark_rawspectra_full(filepath, n_iterations=3, use_index_file=True):
    """Benchmark rawspectra loading all rows."""
    from dysh.fits.gbtfitsload import GBTFITSLoad

    kwargs = {}
    if not use_index_file:
        kwargs["index_file_threshold"] = float('inf')

    sdf = GBTFITSLoad(filepath, **kwargs)
    total_rows = len(sdf._index)

    # Check API signature for compatibility
    sig = inspect.signature(sdf.rawspectra)
    has_fitsindex = "fitsindex" in sig.parameters

    def full_load():
        if has_fitsindex:
            return sdf.rawspectra(bintable=0, fitsindex=0)
        else:
            return sdf.rawspectra(bintable=0)

    result = time_operation(full_load, n_iterations=n_iterations, warmup=1)
    result["total_rows"] = total_rows

    del sdf
    return result


def run_benchmarks(datasets, quick=False, use_index_file=True):
    """Run all benchmarks and return results."""
    results = create_results_dict(quick_mode=quick, extra_metadata={"use_index_file": use_index_file})
    results["datasets"] = {k: {"path": v["path"], "size_mb": v["size_mb"]} for k, v in datasets.items()}

    n_iter = 2 if quick else 5
    n_iter_slow = 1 if quick else 3

    suppress_warnings()

    # Benchmark getps (position switch) - the most common operation
    if "getps" in datasets:
        ds = datasets["getps"]
        logger.info(f"\n=== Position Switch Benchmarks ({ds['size_mb']} MB) ===")

        logger.info("  init...")
        results["benchmarks"]["init"] = benchmark_gbtfitsload_init(
            ds["path"], n_iterations=n_iter_slow, use_index_file=use_index_file
        )

        logger.info("  getps...")
        results["benchmarks"]["getps"] = benchmark_getps(ds, n_iterations=n_iter_slow, use_index_file=use_index_file)

        logger.info("  gettp...")
        results["benchmarks"]["gettp"] = benchmark_gettp(ds, n_iterations=n_iter, use_index_file=use_index_file)

        logger.info("  rawspectra_full...")
        results["benchmarks"]["rawspectra_full"] = benchmark_rawspectra_full(
            ds["path"], n_iterations=n_iter_slow, use_index_file=use_index_file
        )

    # Benchmark getfs (frequency switch)
    if "getfs" in datasets:
        ds = datasets["getfs"]
        logger.info(f"\n=== Frequency Switch Benchmarks ({ds['size_mb']} MB) ===")

        logger.info("  getfs...")
        results["benchmarks"]["getfs"] = benchmark_getfs(ds, n_iterations=n_iter_slow, use_index_file=use_index_file)

    # Benchmark getnod (nodding)
    if "getnod" in datasets:
        ds = datasets["getnod"]
        logger.info(f"\n=== Nodding Benchmarks ({ds['size_mb']} MB) ===")

        logger.info("  getnod...")
        results["benchmarks"]["getnod"] = benchmark_getnod(ds, n_iterations=n_iter_slow, use_index_file=use_index_file)

    # Large file benchmarks (if available and not quick mode)
    if "large" in datasets and not quick:
        ds = datasets["large"]
        logger.info(f"\n=== Large File Benchmarks ({ds['size_mb']} MB) ===")

        logger.info("  init_large...")
        results["benchmarks"]["init_large"] = benchmark_gbtfitsload_init(
            ds["path"], n_iterations=n_iter_slow, use_index_file=use_index_file
        )

        logger.info("  rawspectra_full_large...")
        results["benchmarks"]["rawspectra_full_large"] = benchmark_rawspectra_full(
            ds["path"], n_iterations=n_iter_slow, use_index_file=use_index_file
        )

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Benchmark dysh operations",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument(
        "-o", "--output",
        help="Output JSON file for results"
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
        "--list-datasets",
        action="store_true",
        help="List available test datasets and exit"
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

    datasets = get_test_datasets()

    if args.list_datasets:
        print("Available test datasets:")
        for name, info in datasets.items():
            print(f"  {name}: {info['path']} ({info['size_mb']} MB)")
        return

    if not datasets:
        print("ERROR: No test datasets found. Please ensure testdata directory exists.")
        sys.exit(1)

    use_index_file = not args.no_index_file

    print("dysh Benchmark")
    print("=" * 50)
    print(f"Git: {get_git_info()['branch']} @ {get_git_info()['commit']}")
    print(f"Quick mode: {args.quick}")
    print(f"Use .index file: {use_index_file}")
    print(f"Datasets found: {list(datasets.keys())}")

    results = run_benchmarks(datasets, quick=args.quick, use_index_file=use_index_file)

    print_summary(results)

    if args.output:
        save_results(results, args.output)
    else:
        print("\nTip: Use -o results.json to save results for comparison")


if __name__ == "__main__":
    main()
