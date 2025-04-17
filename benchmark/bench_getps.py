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



progname  = "bench_getps"
benchname = "positionswitch"
data_dir  = "/lma1/teuben/GBT/dysh_data/sdfits/"      # should migrate to use dysh_data or $DYSH_DATA

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog=progname)
    # parser.add_argument("--file",        "-f", action="store",       help="input filename", required=True)
    parser.add_argument("--out",         "-o", action="store",       help="output filename (astropy Table)", required=False)
    parser.add_argument("--append",      "-a", action="store_true",  help="append to previous output file (astropy Table)", required=False)
    parser.add_argument("--overwrite",   "-w", action="store_true",  help="overwrite a previous output file (astropy Table)", required=False)
    parser.add_argument("--profile",     "-p", action="store_true",  help="run the profiler")
    parser.add_argument("--dosomething", "-d", action="store_true",  help="do an optional action")
    parser.add_argument("--skipflags",   "-s", action="store_true",  help="skip reading flags")
    parser.add_argument("--quit",        "-q", action="store_true",  help="quit early")    
    #parser.add_argument("--index", "-i", action="store_true", help="create dysh index table (pandas)")
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

    f1 = dysh_data(example="test1")     # position switch example from notebooks/examples
    print("Loading ",f1)
    sdf1 = GBTFITSLoad(f1, skipflags=args.skipflags)
    time_data.append(time.perf_counter_ns())
    time_stats.append(["load", time_data[-1]])
    if args.dosomething:
        scans = [51,53,55,57]
        scans = [51]
        #sb = sdf1.getps(scan=scans, fdnum=0, ifnum=0, plnum=0)
        #time_data.append(time.perf_counter_ns())
        #time_stats.append(["getps1s", time_data[-1]])
        ps1 = sdf1.getps(scan=scans, fdnum=0, ifnum=0, plnum=0).timeaverage()
        #ps = sb.timeaverage()
        time_data.append(time.perf_counter_ns())
        time_stats.append(["getps1t", time_data[-1]])
        #sb = sdf1.getps(scan=scans, fdnum=0, ifnum=0, plnum=0)
        #time_data.append(time.perf_counter_ns())
        #time_stats.append(["getps2s", time_data[-1]])
        ps2 = sdf1.getps(scan=scans, fdnum=0, ifnum=0, plnum=0).timeaverage()
        #ps = sb.timeaverage()
        time_data.append(time.perf_counter_ns())
        time_stats.append(["getps2t", time_data[-1]])
        #sb = sdf1.getps(scan=scans, fdnum=0, ifnum=0, plnum=0)
        #time_data.append(time.perf_counter_ns())
        #time_stats.append(["getps3s", time_data[-1]])
        ps3 = sdf1.getps(scan=scans, fdnum=0, ifnum=0, plnum=0).timeaverage()
        #ps = sb.timeaverage()
        time_data.append(time.perf_counter_ns())
        time_stats.append(["getps3t", time_data[-1]])
        #sb = sdf1.getps(scan=scans, fdnum=0, ifnum=0, plnum=0)
        #time_data.append(time.perf_counter_ns())
        #time_stats.append(["getps4s", time_data[-1]])
        ps4 = sdf1.getps(scan=scans, fdnum=0, ifnum=0, plnum=0).timeaverage()
        #ps = sb.timeaverage()
        time_data.append(time.perf_counter_ns())
        time_stats.append(["getps4t", time_data[-1]])
      

    time_data.append(time.perf_counter_ns())
    time_stats.append(["end", time_data[-1]])
    print((time_stats[-1][1]-time_stats[0][1])/1e9)
        
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
    print('final',(time_stats[-1][1]-time_stats[0][1])/1e9,'sec')

    for i in range(1,len(time_stats)):
        dt = (time_stats[i][1]-time_stats[i-1][1])/1e6
        label = time_stats[i][0]
        print(label, dt)
