#!/usr/bin/env python
#
import argparse
import cProfile
import os
import pstats
import sys
import time
from pathlib import Path
from pstats import SortKey

import numpy as np
from astropy.table import Table,vstack
import astropy.units as u
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

def filestats(path):
    """Compute the average FITS file size in MB and average flag file length in lines for
       all FITS files and flag files in `path`

       Returns
       -------
            tuple of (# FITS files, FITS size in MB, number of lines)
    """
    if path.is_file():
        pass
    if path.is_dir():
        # get the FITS size 
        nsize = []
        for f in path.glob("*.fits"):
            fstats = os.stat(f)
            nsize.append(fstats.st_size)

        # get the lines in file
        nlines = []
        for f in path.glob("*.flag"):
            with open(f,'rb') as fp:
                nlines.append(sum(1 for _ in fp))

        meandata = np.mean(nsize)*u.byte
        meanlines = int(np.mean(nlines))
        return (len(nsize), meandata.to(u.megabyte).value, meanlines)


progname  = "bench_gbtfitsload"
benchname = "GBTFITSLoad"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog=progname)
    parser.add_argument("--dobench",     "-d", action="store_true",  help="do the benchmark test")
    parser.add_argument("--key",         "-k", action="store",       help="input dysh_data key", default="multismallsmall")
    parser.add_argument("--numfiles",    "-n", action="store",       help="number of SDFITS files to load for multifile data", default=1)
    parser.add_argument("--loop",        "-l", action="store",       help="number of times to loop", default=4)
    parser.add_argument("--skipflags",   "-s", action="store_true",  help="skip reading flags")
    parser.add_argument("--out",         "-o", action="store",       help="output filename (astropy Table)", required=False)
    parser.add_argument("--append",      "-a", action="store_true",  help="append to previous output file (astropy Table)", required=False)
    parser.add_argument("--overwrite",   "-w", action="store_true",  help="overwrite a previous output file (astropy Table)", required=False)
    parser.add_argument("--profile",     "-p", action="store_true",  help="run the profiler")
    parser.add_argument("--statslines",  "-e", action="store",       help="number of profiler statistics lines to print", default=25)
    parser.add_argument("--quit",        "-q", action="store_true",  help="quit early")    
    parser.add_argument("--justtable",   "-j", action="store_true",  help="just print the existin table and exit")    
    #parser.add_argument("--noindex",     "-n", action="store_true",  help="do not create dysh index table (pandas)")
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

    timestr = ""
    i = 0
    
    # output table colnames, units, and dtypes
    table_cols = ["name", "time", "# files", "file size", "flag (lines)", "skipflags"]
    #table_cols = ["name", "time", "skipflags"]
    table_units = ["","ms","","MB","",""]
    table_dtypes = [str, float, int, float, int, str]
    table = Table(names=table_cols, meta={"name": f"Dysh Benchmark {benchname} {args.key}"}, units=table_units, dtype=table_dtypes)
    if args.profile:
        pr = cProfile.Profile()
        pr.enable()
        
    time_stats = []
    time_data = []
    time_data.append(time.perf_counter_ns())
    sk = str(args.skipflags)
    time_stats.append(["start", 0, 0, 0, 0, ""])

    f1 = dysh_data(accept=args.key)
    # use secret GBTFITSLOAD nfiles kwarg to limit number of files loaded.
    nfiles = int(args.numfiles)
    print(f"Loading not more than {nfiles} from {f1}")
    trueNfiles,size_b,nflags= filestats(f1)
    nload = min(nfiles,trueNfiles)
    size_mb = np.round(size_b,2)
    print(f"Will load {nload} of {trueNfiles} files. FITS size per file {size_mb}MB, Flag lines {nflags}")
    if args.dobench:
        for i in range(1,int(args.loop)+1):
            sdf1 = GBTFITSLoad(f1, skipflags=args.skipflags, nfiles=nfiles)
            num_loaded = len(sdf1.files)
            time_data.append(time.perf_counter_ns())
            k = f"load{i}"
            time_stats.append([k, (time_data[-1]-time_data[-2])/1e6, num_loaded, size_mb, nflags, sk])

    time_data.append(time.perf_counter_ns())
    time_stats.append(["end", (time_data[-1]-time_data[-2])/1E6, 0, 0, 0, ""])
    #print("Total (s): ",(time_data[-1]-time_data[0])/1E9)
        
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
                oldtable = Table.read(args.out, format="ascii.ecsv")
                table2 = vstack([oldtable, table])
            elif args.overwrite:
                table2 = table
            else:
                raise Exception(f"{args.out} exists. Use -w to overwrite.")
        else:
            table2 = table
        table2.write(args.out, format="ascii.ecsv", overwrite=True)
    else:
        table.pprint_all()


    # one final time?
    time_data.append(time.perf_counter_ns())
    time_stats.append(["done", time_data[-1]])
    print('final',np.round((time_data[-1]-time_data[0])/1e9,3),'sec')
