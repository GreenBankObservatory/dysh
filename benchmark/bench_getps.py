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
    parser.add_argument("--dobench",     "-d", action="store_true",  help="do the benchmark test")
    parser.add_argument("--key",         "-k", action="store",       help="input dysh_data key", default="test1")
    parser.add_argument("--timeaverage", "-t", action="store_true",  help="time average the Scanblocks to make a Spectrum")
    parser.add_argument("--nocalibrate", "-n", action="store_true",  help="DON'T calibrate the data", default=False)
    parser.add_argument("--loop",        "-l", action="store",       help="number of times to loop", default=4)
    parser.add_argument("--skipflags",   "-s", action="store_true",  help="skip reading flags")
    parser.add_argument("--out",         "-o", action="store",       help="output filename (astropy Table)", required=False)
    parser.add_argument("--append",      "-a", action="store_true",  help="append to previous output file (astropy Table)", required=False)
    parser.add_argument("--overwrite",   "-w", action="store_true",  help="overwrite a previous output file (astropy Table)", required=False)
    parser.add_argument("--profile",     "-p", action="store_true",  help="run the profiler")
    parser.add_argument("--statslines",  "-e", action="store",      help="number of profiler statistics lines to print", default=25)
    parser.add_argument("--quit",        "-q", action="store_true",  help="quit early")    
    parser.add_argument("--justtable",   "-j", action="store_true",  help="just print the existin table and exit")    
    #parser.add_argument("--index", "-i", action="store_true", help="create dysh index table (pandas)")
    args = parser.parse_args()
    print(f"using {args}")

    if args.quit:
        sys.exit(0)

    if args.justtable:
        if args.out is None:
            raise Exception("You must provide the table filename with -o FILENAME")
        table = Table.read(args.out,format="ascii.ecsv")
        table.pprint_all()
        sys.exit(0)

    if args.nocalibrate and args.timeaverage:
        raise Exception("You must calibrate if you want to time average")

    timestr = ""
    i = 0
    
    # output table colnames, units, and dtypes
    table_cols = ["name", "time", "skipflags"]
    table_units = ["","ms",""]
    table_dtypes = [str, float, str]
    table = Table(names=table_cols, meta={"name": f"Dysh Benchmark {benchname}"}, units=table_units, dtype=table_dtypes)
    if args.profile:
        pr = cProfile.Profile()
        pr.enable()
        
    time_stats = []
    time_data = []
    time_data.append(time.perf_counter_ns())
    #time_stats.append(["start", time_data[-1]])
    time_stats.append(["start",  0, ""])

    f1 = dysh_data(example=args.key)     # 'test1' = position switch example from notebooks/examples
    print("Loading ",f1)
    sdf1 = GBTFITSLoad(f1, skipflags=args.skipflags)
    time_data.append(time.perf_counter_ns())
    sk=str(args.skipflags)
    time_stats.append(["load", (time_data[-1]-time_data[-2])/1E6, sk])
    calibrate = not args.nocalibrate
    if args.dobench:
        #scans = [51,53,55,57]
        scans = [51]
        for i in range(1,int(args.loop)+1):
            sb = sdf1.getps(scan=scans, fdnum=0, ifnum=0, plnum=0, calibrate=calibrate)
            time_data.append(time.perf_counter_ns())
            #time_stats.append([f"getps{i}s", time_data[-1]])
            time_stats.append([f"getps{i}s", (time_data[-1]-time_data[-2])/1E6, sk])
            #ps1 = sdf1.getps(scan=scans, fdnum=0, ifnum=0, plnum=0).timeaverage()
            if args.timeaverage:
                ps = sb.timeaverage()
                time_data.append(time.perf_counter_ns())
                #time_stats.append([f"getps{i}t", time_data[-1]])
                time_stats.append([f"getps{i}t", (time_data[-1]-time_data[-2])/1E6, sk])
      

    time_data.append(time.perf_counter_ns())
    time_stats.append(["end-start", (time_data[-1]-time_data[0])/1E6,""])
    #print((time_stats[-1][1]-time_stats[0][1])/1e9)
        
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
    print('final',(time_data[-1]-time_data[0])/1e9,'sec')

    if False:
        for i in range(1,len(time_stats)):
            dt = (time_stats[i][1]-time_stats[i-1][1])/1e6
            label = time_stats[i][0]
            print(label, dt)
