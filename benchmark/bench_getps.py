#!/usr/bin/env python
#

import argparse
import os
import sys
import numpy as np
from astropy.table import Table
from dysh.fits.gbtfitsload import GBTFITSLoad
from dysh.util.files import dysh_data
from dysh.util.timers import DTime

benchname = "getps"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog=sys.argv[0])
    # fmt: off
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
    parser.add_argument("--statslines",  "-e", action="store",       help="number of profiler statistics lines to print", default=25)
    parser.add_argument("--sortkey", "-x", action="store", help="How to sort the profiler statistics, 'cumulative' or 'time'", default="cumulative")
    parser.add_argument("--quit",        "-q", action="store_true",  help="quit early")    
    #parser.add_argument("--index", "-i", action="store_true", help="create dysh index table (pandas)")
    # fmt: on    
    args = parser.parse_args()
    print(f"using {args}")

    if args.quit:
        sys.exit(0)

    if args.nocalibrate and args.timeaverage:
        raise Exception("You must calibrate if you want to time average")

    data_cols  = ["skipflags"]
    data_units = [""]
    data_types = [str]
    dt = DTime(benchname=benchname,
               data_cols=data_cols, data_units=data_units, data_types=data_types, # no data=[] supported yet
               args=vars(args))

    sk=str(args.skipflags)

    # reading dataset-1

    f1 = dysh_data(example=args.key)     # 'test1' = position switch example from notebooks/examples
    print("Loading ",f1)
    sdf1 = GBTFITSLoad(f1, skipflags=args.skipflags)
    print('STATS:',sdf1.stats())    
    dt.tag("load", [sk])
    calibrate = not args.nocalibrate
    if args.dobench:
        scans = [51,53,55,57]
        #scans = [51]
        for i in range(1,int(args.loop)+1):
            if not args.timeaverage:
                sb = sdf1.getps(scan=scans, fdnum=0, ifnum=0, plnum=0, calibrate=calibrate)
                dt.tag(f"getps{i}s", [sk])
                #ps1 = sdf1.getps(scan=scans, fdnum=0, ifnum=0, plnum=0).timeaverage()
            else:
                ps = sdf1.getps(scan=scans, fdnum=0, ifnum=0, plnum=0, calibrate=calibrate).timeaverage()
                dt.tag(f"getps{i}t",[sk])

    # close data, do some other silly work
    if False:
        del sdf1
        sdf1 = np.arange(1e5)
        dt.tag("arange 1e5", [sk])
        sdf1 = np.arange(1e6)
        dt.tag("arange 1e6", [sk])
        sdf1 = np.arange(1e7)
        dt.tag("arange 1e7", [sk])
        sdf1 = np.arange(1e8)
        dt.tag("arange 1e8", [sk])
        sdf1 = np.arange(1e5)
        dt.tag("arange 1e5", [sk])

        # read it one more time

        f1 = dysh_data(example=args.key)     # 'test1' = position switch example from notebooks/examples
        print("Loading ",f1)
        sdf1 = GBTFITSLoad(f1, skipflags=args.skipflags)
        dt.tag("load", [sk])
        calibrate = not args.nocalibrate
        if args.dobench:
            #scans = [51,53,55,57]
            scans = [51]
            #sb = list(range(int(args.loop)+1))
            for i in range(1,int(args.loop)+1):
                if not args.timeaverage:
                    sb  = sdf1.getps(scan=scans, fdnum=0, ifnum=0, plnum=0, calibrate=calibrate)
                    dt.tag(f"getps{i}s", [sk])
                    #ps1 = sdf1.getps(scan=scans, fdnum=0, ifnum=0, plnum=0).timeaverage()
                else:
                    ps = sdf1.getps(scan=scans, fdnum=0, ifnum=0, plnum=0, calibrate=calibrate).timeaverage()
                    dt.tag(f"getps{i}t",[sk])

        
    dt.tag('report',[sk])
    dt.close()
    dt.report()

    print('final',dt.total()/1000,' sec')
    
