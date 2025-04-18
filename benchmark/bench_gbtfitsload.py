#!/usr/bin/env python
#
import argparse
import cProfile
import os
import pstats
import sys
import time
from pstats import SortKey

import numpy as np
from astropy.table import Table,vstack
from dysh.fits.gbtfitsload import GBTFITSLoad
from dysh.util.files import dysh_data

def mkdir(name, clean=True):
    """ simpler frontend for making a directory that might also already exist
        clean=True:    also remove files inside
    """
    os.makedirs(name, exist_ok = True)
    if clean:
        fns = os.listdir(name)
        for fn in fns:
            print(f"Removing {fn} from {name}")
            os.unlink(os.path.join(name,fn))



progname  = "bench_gbtfitsload"
benchname = "GBTFITSLoad"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog=progname)
    parser.add_argument("--dobench",     "-d", action="store_true",  help="do the benchmark test")
    parser.add_argument("--key",         "-k", action="store",       help="input dysh_data key. ", default="multismall")
    parser.add_argument("--numfiles",    "-n", action="store",       help="number of SDFITS files to load for multifile data", default=1)
    parser.add_argument("--loop",        "-l", action="store",       help="number of times to loop", default=4)
    parser.add_argument("--skipflags",   "-s", action="store_true",  help="skip reading flags")
    parser.add_argument("--out",         "-o", action="store",       help="output filename (astropy Table)", required=False)
    parser.add_argument("--append",      "-a", action="store_true",  help="append to previous output file (astropy Table)", required=False)
    parser.add_argument("--overwrite",   "-w", action="store_true",  help="overwrite a previous output file (astropy Table)", required=False)
    parser.add_argument("--profile",     "-p", action="store_true",  help="run the profiler")
    parser.add_argument("--statslines",  "-e", action="store",       khelp="number of profiler statistics lines to print", default=25)
    parser.add_argument("--quit",        "-q", action="store_true",  help="quit early")    
    #parser.add_argument("--noindex",     "-n", action="store_true",  help="do not create dysh index table (pandas)")
    args = parser.parse_args()
    print(f"using {args}")

    if args.quit:
        sys.exit(0)

    if args.nocalibrate and args.timeaverage:
        raise Exception("You must calibrate if you want to time average")

    timestr = ""
    i = 0
    
    # output table colnames, units, and dtypes
    table_cols = ["name", "time"]
    table_units = []
    table_dtypes = [str, int]
    table = Table(names=table_cols, meta={"name": f"Dysh Benchmark {benchname} {args.key}"}, units=table_units, dtype=table_dtypes)
    if args.profile:
        pr = cProfile.Profile()
        pr.enable()
        
    time_stats = []
    time_data = []
    time_data.append(time.perf_counter_ns())
    time_stats.append(["start", time_data[-1]])

    f1 = dysh_data(acceptance=args.key)
    # use secret GBTFITSLOAD nfiles kwarg to limit number of files loaded.
    nfiles = int(args.numfiles)
    print(f"Loading {nfiles} from {f1}")
    if args.dotest:
        for i in range(1,int(args.loop)+1):
            sdf1 = GBTFITSLoad(f1, skipflags=args.skipflags, nfiles=nfiles)
            time_data.append(time.perf_counter_ns())
            time_stats.append([f"load{i}-{nfiles}", time_data[-1]])

    time_data.append(time.perf_counter_ns())
    time_stats.append(["end", time_data[-1]])
    print((time_stats[-1][1]-time_stats[0][1])/1e9)
        
    if args.profile:
        pr.disable()
        ps = pstats.Stats(pr).sort_stats(SortKey.CUMULATIVE)
        #ps = pstats.Stats(pr).sort_stats(SortKey.CUMULATIVE, SortKey.TIME)
        ps.print_stats(int(args.statslines))

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
    print('final',(time_stats[-1][1]-time_stats[0][1])/1e9,'sec')

    for i in range(1,len(time_stats)):
        dt = (time_stats[i][1]-time_stats[i-1][1])/1e6
        label = time_stats[i][0]
        print(label, dt)
