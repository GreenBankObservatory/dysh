#!/usr/bin/env python
import argparse
import cProfile
import os
import pstats
import sys
import time
from pstats import SortKey

from astropy.table import Table,vstack
import numpy as np
from dysh.fits.gbtfitsload import GBTFITSLoad

progname="perf_skel"
data_dir = "/data/gbt/examples/"
benchname = "FOOBAR BENCH"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog=progname)
    parser.add_argument("--file", "-f", action="store", help="input filename", required=True)
    parser.add_argument("--out", "-o", action="store", help="output filename (astropy Table)", required=False)
    parser.add_argument("--append", "-a", action="store", help="append to previous output file (astropy Table)", required=False)
    parser.add_argument("--overwrite", "-w", action="store", help="overwrite a previous output file (astropy Table)", required=False)
    parser.add_argument("--profile", "-p", action="store_true", help="run the profiler")
    parser.add_argument("--dosomething", "-d", action="store_true", help="do an optional action")
    #parser.add_argument("--index", "-i", action="store_true", help="create dysh index table (pandas)")
    args = parser.parse_args()
    print(f"using {args}")

    timestr = ""
    i = 0
    # output table colnames, units, and dtypes
    table_cols = [ ]
    table_units = []
    table_dtypes = []
    table = Table(names=table_cols, meta={"name": f"Dysh Benchmark {benchname}"}, units=table_units, dtype=table_dtypes)
    pr = cProfile.Profile()
    time_stats = {}
    time_data = []
    time_data.append(time.perf_counter_ns())
    time_stats["start"] = time_data[-1]
    pr.enable()
    sdf = GBTFITSLoad(data_dir+args.file)
    time_data.append(time.perf_counter_ns())
    time_stats["load"] =  time_data[-1]
    if args.dosomething:
        # do someting
        time_data.append(time.perf_counter_ns())
        time_stats["something"] =  time_data[-1]
    pr.disable()
    # table.add_row([...])

    ps = pstats.Stats(pr).sort_stats(SortKey.CUMULATIVE)
    ps.print_stats(25)

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
