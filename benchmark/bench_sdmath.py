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

def initps(n, mode=1, dtype=np.float64):
    # initialize a spectrum
    if mode == 1:
        # random but expensive for large n
        r1 = np.random.normal(0,1,n).astype(dtype)
        r2 = np.random.normal(0,1,n).astype(dtype)
        r3 = np.random.normal(0,1,n).astype(dtype)
        r4 = np.random.normal(0,1,n).astype(dtype)
    else:
        # silly sequence, reproducable between C and Python
        r1 = np.arange(0.0*n,1.0*n,dtype=dtype)
        r2 = np.arange(1.0*n,2.0*n,dtype=dtype)
        r3 = np.arange(2.0*n,3.0*n,dtype=dtype)
        r4 = np.arange(3.0*n,4.0*n,dtype=dtype)
    return r1,r2,r3,r4

def initps2(nscan,nchan, mode=1, dtype=np.float64):
    # initialize a spectrum
    # random but expensive for large n
    if mode == -1:
        data = np.random.normal(0,1,(nscan,4,nchan)).astype(dtype)
    else:
        data = np.arange(0.0,nscan*4*nchan, dtype=dtype).reshape(nscan,4,nchan)
    return data

def getps(r1,r2,r3,r4):
    # basic match to get a 4-phase TP into a PS spectrum
    # r1,r2,r3,r4 =  on_cold, on_hot, off_cold, off_hot
    tc = 10.0
    tsys = tc * r1.sum()/(r2-r1).sum()
    ta = tsys * ( (r1+r2)/(r3+r4) - 1)
    return ta

def getps2(data):
    nscan,ndim,nchan = data.shape
    #print(nscan,ndim,nchan)
    tc = 10.0
    tsys = tc * data[:,0,:].sum(axis=1) / (data[:,1,:] - data[:,0,:]).sum(axis=1)
    ta = tsys[:,np.newaxis] * ( (data[:,0,:] + data[:,1,:]) /  (data[:,2,:] + data[:,3,:]) - 1)
    return ta
    
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog=sys.argv[0])
    # fmt: off
    parser.add_argument("--dobench",     "-d", action="store_true",  help="do the benchmark test")

    parser.add_argument("--loop",        "-l", action="store",       help="number of times to loop", default=4)
    parser.add_argument("--nchan",       "-n", action="store",       help="number of channels", default=100000)
    parser.add_argument("--nscan",       "-s", action="store",       help="number of scans", default=1000)
    parser.add_argument("--mode",        "-m", action="store",       help="initialization mode", default=1)
    parser.add_argument("--timeaverage", "-t", action="store_true",  help="time average as well")

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

    data_cols  = ["nscan", "nchan", "mode"]
    data_units = ["", "", ""]
    data_types = [int, int, int]
   
    dt = DTime(benchname=benchname,
               data_cols=data_cols, data_units=data_units, data_types=data_types,
               args=vars(args))

    nchan = int(args.nchan)
    nscan = int(args.nscan)
    mode = int(args.mode)
    data = [nscan, nchan, mode]

    dtype = np.float32

    if mode > 0:
        # dangerous, might run quicker since there's no block of data, only 4 rows
        r1,r2,r3,r4 = initps(nchan, mode, dtype)   # 224 -> 146
        dt.tag("init", data)

        if args.dobench:
            for i in range(1,int(args.loop)+1):
                for j in range(int(args.nscan)):
                    sp = getps(r1,r2,r3,r4)
                dt.tag(f"math_{i}", data)
            print("mean:",sp.sum())
            print('data:',r1[0],r1[1],r4[-1])
            print('type:',type(r1[0]),type(sp.sum()))
    else:
        # trying a full block of nscan x nchan data
        raw = initps2(nscan, nchan, mode, dtype)    # 680 -> 355
        dt.tag("init",data)

        if args.dobench:
            for i in range(1,int(args.loop)+1):
                if not args.timeaverage:
                    ta =  getps2(raw)
                    dt.tag(f"math2_{i}", data)
                else:
                    taver = getps2(raw).mean(axis=0)
                    dt.tag(f"aver2_{i}", data)
            if not args.timeaverage:
                print('mean:',ta.sum())
            else:
                print('mean:',taver.sum())
            print('data:',raw[0][0][0], raw[0][0][1], raw[-1][-1][-1])
            print('type:',type(raw[0][0][0]))
            
    dt.tag('report', data)
    dt.close()
    dt.report()

    print('final',dt.total()/1000,' sec')
    
