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
from dysh.util.timers import DTime

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
    parser.add_argument("--justtable",   "-j", action="store_true",  help="just print the existing table and exit")    
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

    if args.profile:
        pr = cProfile.Profile()
        pr.enable()

        
    data_cols  = ["skipflags"]
    data_units = [""]
    data_types = [str]
    dt = DTime(out=args.out,
               data_cols=data_cols, data_units=data_units, data_types=data_types)  # no data=[] supported yet

    sk=str(args.skipflags)

    f1 = dysh_data(example=args.key)     # 'test1' = position switch example from notebooks/examples
    print("Loading ",f1)
    sdf1 = GBTFITSLoad(f1, skipflags=args.skipflags)
    dt.tag("load", [sk])
    calibrate = not args.nocalibrate
    if args.dobench:
        #scans = [51,53,55,57]
        scans = [51]
        for i in range(1,int(args.loop)+1):
            sb = sdf1.getps(scan=scans, fdnum=0, ifnum=0, plnum=0, calibrate=calibrate)
            dt.tag(f"getps{i}s", [sk])
            #ps1 = sdf1.getps(scan=scans, fdnum=0, ifnum=0, plnum=0).timeaverage()
            if args.timeaverage:
                ps = sb.timeaverage()
                dt.tag(f"getps{i}t",[sk])
    dt.tag("done",[sk])
        
    if args.profile:
        pr.disable()
        ps = pstats.Stats(pr).sort_stats(SortKey.CUMULATIVE)
        #ps = pstats.Stats(pr).sort_stats(SortKey.CUMULATIVE, SortKey.TIME)
        ps.print_stats(int(args.statslines))

    dt.tag('report',[sk])
    dt.close()
    dt.report()

    # report total CPU time in sec.
    # This does not include startup time before the DTime() [~2 sec]
    total = dt.total()
    print('final',total/1000,' sec')
    
