#!/usr/bin/env python
"""
Benchmark script for dysh operations.

This script measures the performance of flag operations on various size SDFITS files.
 It works on both main and feature branches - run it in different checkouts
and compare the JSON outputs to measure speedups.

Usage:
    # Run benchmark and save results
    python benchmark/bench_flag.py -o results.json

    # Compare two result files
    python benchmark/bench_flag.py --compare results_main.json results_branch.json

    # Quick benchmark (fewer iterations)
    python benchmark/bench_flag.py --quick -o results.json
"""

import argparse
import inspect
import sys
import numpy as np
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

global_sdf = None

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

    from dysh.util.files import dysh_data

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
            "multi": False,
        }

    # Large file for rawspectra benchmarks
    large_file = testdata / "TGBT25B_603_12/TGBT25B_603_12.raw.vegas/TGBT25B_603_12.raw.vegas.A.fits"
    if large_file.exists():
        datasets["large"] = {
            "path": str(large_file),
            "size_mb": round(large_file.stat().st_size / (1024 * 1024), 1),
            "scans": [152],
            "multi": False,
        }

    # Huge file for all benchmarks except nod
    huge_file = dysh_data(example="mixed-fs-ps/data/AGBT16B_225_05/AGBT16B_225_05.raw.vegas")
    if huge_file.exists():
        datasets["huge"] = {
            "path": str(huge_file),
            "size_mb": round(huge_file.stat().st_size / (1024 * 1024), 1),
            "huge": True,
            "multi": True,
            "scans":np.arange(6,42).tolist(),
        }

    # Medium file
    medium_file = testdata / "AGBT05B_047_01/AGBT05B_047_01.raw.acs/AGBT05B_047_01.raw.acs.fits"
    if medium_file.exists():
        datasets["medium"] = {
            "path": str(medium_file),
            "size_mb": round(medium_file.stat().st_size / (1024 * 1024), 1),
            "scans": [51,52,53,54,55,56,57,58],
            "multi": False,
        }

    # Nodding dataset (from nodding.ipynb) - check for testtrim version
    nod_file = testdata / "TGBT22A_503_02/TGBT22A_503_02.raw.vegas/TGBT22A_503_02.raw.vegas.testtrim.fits"
    if nod_file.exists():
        datasets["getnod"] = {
            "path": str(nod_file),
            "size_mb": round(nod_file.stat().st_size / (1024 * 1024), 1),
            "scans": [62,63],
            "ifnum": 0,
            "plnum": 0,
            "multi": False,
        }

    # Frequency switch - TGBT21A_504_01 has frequency switched data
    fs_file = testdata / "TGBT21A_504_01/TGBT21A_504_01.raw.vegas/TGBT21A_504_01.raw.vegas.A.fits"
    if fs_file.exists():
        datasets["getfs"] = {
            "path": str(fs_file),
            "size_mb": round(fs_file.stat().st_size / (1024 * 1024), 1),
            "scans": [20],
            "ifnum": 0,
            "plnum": 0,
            "fdnum": 0,
            "multi": False,
        }

    return datasets



def benchmark_gbtfitsload_init(filepath, n_iterations=3, use_index_file=False, backend=None):
    """Benchmark GBTFITSLoad initialization."""
    from dysh.fits.gbtfitsload import GBTFITSLoad

    kwargs = {}
    if not use_index_file:
        kwargs["index_file_threshold"] = float("inf")

    def load_file():
        global global_sdf
        global_sdf = GBTFITSLoad(filepath,skipflags=True, flag_vegas=False, fitsbackend=backend, **kwargs)
        return global_sdf

    return time_operation(load_file, n_iterations=n_iterations, warmup=1)

def benchmark_write_flag(dataset, sdf, out, n_iterations=3, multi=False):
    """Benchmark flag operation ."""
    def do_write():
        sdf.write(out, overwrite=True, flags=True, multifile=multi)
    result = time_operation(do_write, n_iterations=n_iterations, warmup=1)
    result["scans"] = dataset.get("scans")
    result["file_size_mb"] = dataset["size_mb"]
    return result

def benchmark_flag(dataset, n_iterations=3, use_index_file=False, backend=None):
    """Benchmark flag operation ."""

    scans = dataset.get("scans") # should error out if scans not defined

    kwargs = {}
    if not use_index_file:
        kwargs["index_file_threshold"] = float("inf")

    def do_flag():
        global global_sdf
        chans = global_sdf.nchan(0)
        global_sdf.flag_channel(channel=[[0,100],[chans-100,chans-1]])
        global_sdf.flag(scan=scans[0],channel=np.arange(0,chans,11).tolist())
        global_sdf.apply_flags()

    result = time_operation(do_flag, n_iterations=n_iterations, warmup=1)
    result["scans"] = scans
    result["file_size_mb"] = dataset["size_mb"]
    return result


def run_benchmarks(datasets, fits, quick=False, use_index_file=False, backend=None, huge=False, 
                   multifile=False):
    """Run all benchmarks and return results."""
    results = create_results_dict(
        quick_mode=quick, extra_metadata={"use_index_file": use_index_file, "FITSBackend": backend}
    )
    results["datasets"] = {k: {"path": v["path"], "size_mb": v["size_mb"]} for k, v in datasets.items()}

    n_iter = 2 if quick else 5
    n_iter_slow = 1 if quick else 3

    suppress_warnings()

    # Benchmark all files
    for k in datasets:
        if k == "huge" and not huge:
            continue
        ds = datasets[k]
        logger.info(f"\n=== {k} Flag Benchmarks ({ds['size_mb']} MB) ===")

        logger.info("  init...")
        results["benchmarks"][k+" init"] = benchmark_gbtfitsload_init(
            ds["path"], n_iterations=n_iter_slow, use_index_file=use_index_file, backend=backend
        )

        logger.info("  flag...")
        results["benchmarks"][k+" flag"] = benchmark_flag(
            ds, n_iterations=n_iter_slow, use_index_file=use_index_file, backend=backend
        )
        logger.info("  write flag...")
        results["benchmarks"][k+" write"] = benchmark_write_flag(ds,
            global_sdf,fits, n_iterations=n_iter_slow, multi=multifile, 
        )


    return results


def main():
    parser = argparse.ArgumentParser(
        description="Benchmark dysh operations", formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__
    )
    parser.add_argument("-o", "--output", help="Output JSON file for results")
    parser.add_argument("--quick", action="store_true", help="Run quick benchmark with fewer iterations")
    parser.add_argument(
        "--compare",
        nargs=2,
        metavar=("BASELINE", "CANDIDATE"),
        help="Compare two result files instead of running benchmarks",
    )
    parser.add_argument("--list-datasets", action="store_true", help="List available test datasets and exit")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging")
    parser.add_argument(
        "--no-index-file", action="store_true", help="Disable .index file usage (force reading from FITS)"
    )
    parser.add_argument(
        "--backend",
        "-b",
        action="store",
        help="FITSBackend to use for getting raw spectra. either 'fitsio', 'astropy', or None",
        default=None,
    )
    parser.add_argument(
        "--huge",
        action="store_true",
        help="Include huge file (77GB) in testing.  Cannot be used with --quick",
        default=False
    )  
    parser.add_argument(
        "--multifile",
        "-m",
        action="store_true",
        help="When writing an SDF that has multiple input files, write multiple output files",
        default=False
    ) 
    parser.add_argument(
        "--fits",
        "-f",
        action="store",
        help="Output SDFITS file(s). If 'multifile=True', outputs will be e.g., file1.fits, file2.fits, etc.",
        required=True
    ) 
    args = parser.parse_args()
    if args.quick and args.huge:
        print("ERROR: Options --quick and --huge are mutually exclusive")
        sys.exit(1)

    fbe = str(args.backend)

    setup_logging(verbose=args.verbose)

    if args.compare:
        compare_results(args.compare[0], args.compare[1])
        return

    datasets = get_test_datasets()

    if args.list_datasets:
        print("Available test datasets:")
        for name, info in datasets.items():
            if name == "huge" and not args.huge:
                continue
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
    print(f"Huge mode: {args.huge}")
    print(f"Multifile: {args.multifile}")
    print(f"FITS Backend: {fbe}")
    print(f"Datasets found: {list(datasets.keys())}")

    results = run_benchmarks(
        datasets, args.fits, quick=args.quick, use_index_file=use_index_file, 
        backend=args.backend, huge=args.huge, multifile=args.multifile
    )

    if args.verbose:
        print(f"{results=}")
    print_summary(results)

    if args.output:
        save_results(results, args.output)
    else:
        print("\nTip: Use -o results.json to save results for comparison")


if __name__ == "__main__":
    main()
