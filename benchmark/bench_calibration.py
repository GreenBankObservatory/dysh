#!/usr/bin/env python
#
# Test per-spectrum overhead of the calibration pipeline: first vs warm getps
# (lazy-load penalty), getspec/make_spectrum with and without WCS, timeaverage,
# and common Spectrum operations. Companion to bench_getps.py; used to produce
# before/after numbers for the performance PR series (see workflows/README.md).

import argparse
import sys

from dysh.fits.gbtfitsload import GBTFITSLoad
from dysh.util.files import dysh_data
from dysh.util.timers import DTime

benchname = "calibration"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog=sys.argv[0])
    # fmt: off
    parser.add_argument("--key",         "-k", action="store",       help="input dysh_data key", default="getps")
    parser.add_argument("--scan",        "-c", action="store",       help="scan number for getps", default=51, type=int)
    parser.add_argument("--loop",        "-l", action="store",       help="number of warm getps+timeaverage loops", default=4, type=int)
    parser.add_argument("--skipflags",   "-s", action="store_true",  help="skip reading flags")
    parser.add_argument("--out",         "-o", action="store",       help="output filename (astropy Table)", required=False)
    parser.add_argument("--append",      "-a", action="store_true",  help="append to previous output file (astropy Table)", required=False)
    parser.add_argument("--overwrite",   "-w", action="store_true",  help="overwrite a previous output file (astropy Table)", required=False)
    parser.add_argument("--profile",     "-p", action="store_true",  help="run the profiler")
    parser.add_argument("--statslines",  "-e", action="store",       help="number of profiler statistics lines to print", default=25)
    parser.add_argument("--sortkey",     "-x", action="store",       help="How to sort the profiler statistics, 'cumulative' or 'time'", default="cumulative")
    parser.add_argument("--memory",      "-m", action="store_true",  help="track memory usage")
    # fmt: on
    args = parser.parse_args()
    print(f"using {args}")

    data_cols = ["nspec"]
    data_units = [""]
    data_types = [int]
    dt = DTime(benchname=benchname, data_cols=data_cols, data_units=data_units, data_types=data_types, args=vars(args))

    nspec = 0
    dt.tag("init", [nspec])

    f1 = dysh_data(example=args.key)
    print("Loading ", f1)
    sdf = GBTFITSLoad(f1, skipflags=args.skipflags)
    dt.tag("load", [nspec])

    # First getps pays any lazy-load penalty; subsequent ones are warm.
    sb = sdf.getps(scan=args.scan, fdnum=0, ifnum=0, plnum=0)
    nspec = sb[0].nint
    dt.tag("getps_first", [nspec])

    for i in range(1, args.loop + 1):
        sb = sdf.getps(scan=args.scan, fdnum=0, ifnum=0, plnum=0)
        dt.tag(f"getps_warm{i}", [nspec])

    # Per-spectrum construction cost: WCS + SkyCoord + meta deepcopy per call.
    scan0 = sb[0]
    for i in range(nspec):
        s = scan0.getspec(i)
    dt.tag("getspec_wcs", [nspec])

    for i in range(nspec):
        s = scan0.getspec(i, use_wcs=False)
    dt.tag("getspec_nowcs", [nspec])

    for i in range(1, args.loop + 1):
        ta = sb.timeaverage()
        dt.tag(f"timeaverage{i}", [nspec])

    # Common Spectrum operations on the time-averaged result.
    ta2 = sb.timeaverage()
    avg = ta.average(ta2)
    dt.tag("average_spectra", [nspec])

    st = ta.stats()
    dt.tag("stats", [nspec])

    sm = ta.smooth(kernel="gaussian", width=16)
    dt.tag("smooth", [nspec])

    ta.baseline(degree=2, remove=False)
    dt.tag("baseline", [nspec])

    dt.tag("report", [nspec])
    dt.close()
    dt.report()

    print("final", dt.total() / 1000, " sec")
