#!/usr/bin/env python
# Test performance of data writing using the EDGE data.  This is not
# a standard benchmark test.
import argparse
import cProfile
import os
import pstats
import sys
import time
from pstats import SortKey

from astropy.table import Table, vstack

from dysh.fits.gbtfitsload import GBTFITSLoad


def mkdir(name, clean=True):
    """simpler frontend for making a directory that might also already exist
    clean=True:    also remove files inside
    """
    os.makedirs(name, exist_ok=True)
    if clean:
        fns = os.listdir(name)
        for fn in fns:
            print(f"Removing {fn} from {name}")
            os.unlink(os.path.join(name, fn))


progname = "bench_edge"
benchname = "edge L-band"
data_dir = "/lma1/teuben/GBT/dysh_data/sdfits/"  # should migrate to use dysh_data or $DYSH_DATA

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog=progname)
    parser.add_argument("--file", "-f", action="store", help="input filename", required=True)
    parser.add_argument("--out", "-o", action="store", help="output filename (astropy Table)", required=False)
    parser.add_argument(
        "--append", "-a", action="store_true", help="append to previous output file (astropy Table)", required=False
    )
    parser.add_argument(
        "--overwrite",
        "-w",
        action="store_true",
        help="overwrite a previous output file (astropy Table)",
        required=False,
    )
    parser.add_argument("--profile", "-p", action="store_true", help="run the profiler")
    parser.add_argument("--dosomething", "-d", action="store_true", help="do an optional action")
    parser.add_argument("--skipflags", "-s", action="store_true", help="skip reading flags")
    parser.add_argument("--quit", "-q", action="store_true", help="quit early")
    # parser.add_argument("--index", "-i", action="store_true", help="create dysh index table (pandas)")
    args = parser.parse_args()
    print(f"using {args}")

    if args.quit:
        sys.exit(0)

    timestr = ""
    i = 0

    # output table colnames, units, and dtypes
    table_cols = ["name", "time"]
    table_units = []
    table_dtypes = [str, int]
    table = Table(names=table_cols, meta={"name": f"Dysh Benchmark {benchname}"}, units=table_units, dtype=table_dtypes)
    if args.profile:
        pr = cProfile.Profile()
        pr.enable()

    time_stats = []
    time_data = []
    time_data.append(time.perf_counter_ns())
    time_stats.append(["start", time_data[-1]])

    skipflags = args.skipflags
    print("Loading ", data_dir + args.file)
    print("Skipflags ", skipflags)

    sdf1 = GBTFITSLoad(data_dir + args.file, skipflags=skipflags)
    time_data.append(time.perf_counter_ns())
    time_stats.append(["load", time_data[-1]])
    if args.dosomething:
        # do someting
        mkdir("edge1")
        scans1 = [56, 57, 58]
        sdf1.write("edge1/file.fits", scan=scans1, ifnum=1, plnum=0, intnum=0, overwrite=True)
        #  6 rows
        time_data.append(time.perf_counter_ns())
        time_stats.append(["write1", time_data[-1]])

        mkdir("edge2")
        scans2 = list(range(56, 65))
        sdf1.write("edge2/file.fits", scan=scans2, ifnum=1, plnum=0, overwrite=True)
        #  3090 rows
        time_data.append(time.perf_counter_ns())
        time_stats.append(["write2", time_data[-1]])

    time_data.append(time.perf_counter_ns())
    time_stats.append(["end", time_data[-1]])
    print((time_stats[-1][1] - time_stats[0][1]) / 1e9)

    if args.profile:
        pr.disable()
        ps = pstats.Stats(pr).sort_stats(SortKey.CUMULATIVE)
        ps.print_stats(25)

    for i in range(len(time_stats)):
        table.add_row(time_stats[i])

    for c in table.columns:
        if table[c].info.dtype == float:
            table[c].info.format = "0.1f"

    if args.out is not None:
        if os.path.exists(args.out):
            if args.append:
                oldtable = Table.read(args.out, format="ipac")
                table2 = vstack([oldtable, table])
            elif args.overwrite:
                table2 = table
            else:
                raise Exception(f"{args.out} exists. Use -w to overwrite.")
        else:
            table2 = table
        table2.write(args.out, format="ipac", overwrite=True)
    else:
        table.pprint_all()

    # one final time?
    time_data.append(time.perf_counter_ns())
    time_stats.append(["done", time_data[-1]])
    print((time_stats[-1][1] - time_stats[0][1]) / 1e9)
