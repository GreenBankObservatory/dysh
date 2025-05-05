#!/usr/bin/env python
#
#       

import argparse
import os
import sys
import numpy as np
from astropy.table import Table
from dysh.util.timers import DTime

benchname = "sdmath"

def initps(n, mode=1):
    if mode == 1:
        r1 = np.random.normal(0,1,n)
        r2 = np.random.normal(0,1,n)
        r3 = np.random.normal(0,1,n)
        r4 = np.random.normal(0,1,n)
    else:
        r1 = np.arange(0*n,1*n)
        r2 = np.arange(1*n,2*n)
        r3 = np.arange(2*n,3*n)
        r4 = np.arange(3*n,4*n)
    return r1,r2,r3,r4

def getps(r1,r2,r3,r4):
    tc = 10.0
    tsys = tc * r1.mean()/(r1.mean()-r2.mean())
    ta = tsys * ( (r1+r2)/(r3+r4) - 1)
    return ta

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog=sys.argv[0])
    # fmt: off
    parser.add_argument("--dobench",     "-d", action="store_true",  help="do the benchmark test")

    parser.add_argument("--loop",        "-l", action="store",       help="number of times to loop", default=4)
    parser.add_argument("--nchan",       "-n", action="store",       help="number of channels", default=10000)
    parser.add_argument("--mode",        "-m", action="store",       help="initialization mode", default=1)

    
    parser.add_argument("--out",         "-o", action="store",       help="output filename (astropy Table)", required=False)
    parser.add_argument("--append",      "-a", action="store_true",  help="append to previous output file (astropy Table)", required=False)
    parser.add_argument("--overwrite",   "-w", action="store_true",  help="overwrite a previous output file (astropy Table)", required=False)
    parser.add_argument("--profile",     "-p", action="store_true",  help="run the profiler")
    parser.add_argument("--statslines",  "-e", action="store",       help="number of profiler statistics lines to print", default=25)
    parser.add_argument("--quit",        "-q", action="store_true",  help="quit early")    
    # fmt: on    
    args = parser.parse_args()
    print(f"using {args}")

    if args.quit:
        sys.exit(0)

    data_cols  = ["nchan", "mode"]
    data_units = ["", ""]
    data_types = [int, int]
   
    dt = DTime(benchname=benchname,
               data_cols=data_cols, data_units=data_units, data_types=data_types,
               args=vars(args))

    n = int(args.nchan)
    mode = int(args.mode)
    data = [n, mode]
    r1,r2,r3,r4 = initps(n, mode)
    dt.tag("init", data)

    if args.dobench:
        for i in range(1,int(args.loop)+1):
            sp = getps(r1,r2,r3,r4)
            dt.tag(f"math{i}", data)
        print("sum:",sp.sum())
        
    dt.tag('report', data)
    dt.close()
    dt.report()

    print('final',dt.total()/1000,' sec')
    
