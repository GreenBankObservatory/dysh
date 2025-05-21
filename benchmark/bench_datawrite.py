#!/usr/bin/env python
#
import argparse
import os
import sys

import numpy as np
import astropy.units as u
from astropy.table import Table
from dysh.fits.gbtfitsload import GBTFITSLoad
from dysh.util.files import dysh_data
from dysh.util.timers import DTime


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
            with open(f, "rb") as fp:
                nlines.append(sum(1 for _ in fp))

        meandata = np.mean(nsize) * u.byte
        meanlines = int(np.mean(nlines))
        return (len(nsize), meandata.to(u.megabyte).value, meanlines)


progname = "bench_datawrite"
benchname = "DataWrite"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog=progname)
    #parser.add_argument("--key", "-k", action="store", help="input dysh_data key", default="multismallsmall")
    #parser.add_argument(
    #    "--numfiles", "-n", action="store", help="number of SDFITS files to load for multifile data", default=1
    #)
    parser.add_argument("--loop", "-l", action="store", help="number of times to loop", default=4)
    parser.add_argument("--skipflags", "-s", action="store_true", help="skip reading/writing flags")
    parser.add_argument("--writedata",  action="store", help="One of 'sdf' or 'sb' to write SDFITS or ScanBlock data")
    parser.add_argument("--rows", "-r", action="store", help="If `writedata` is sdf, the number of rows to write.")
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
    parser.add_argument(
        "--statslines", "-e", action="store", help="number of profiler statistics lines to print", default=25
    )
    parser.add_argument("--sortkey", "-x", action="store", help="How to sort the profiler statistics, 'cumulative' or 'time'", default="cumulative")
    parser.add_argument("--memory",  "-m", action="store_true",  help="track memory usage")
    parser.add_argument("--quit",    "-q", action="store_true", help="quit early")
    # parser.add_argument("--noindex",     "-n", action="store_true",  help="do not create dysh index table (pandas)")
    args = parser.parse_args()

    valid_write = ['sdf','sb']
    if args.writedata not in valid_write:
        raise ValueError(f"writedata must be one of {valid_write}")

    if args.writedata == 'sdf' and args.rows is None:
        raise ValueError("You must supply number of rows to write (-r)  for SDFITS write test.")

    if args.quit:
        sys.exit(0)

    # output table colnames, units, and dtypes
    # DTime automatically handles name and time, so just the additional columns go here.
    data_cols  = ["#files", "file_size", "totsize", "nchan", "nrow", "nIF", "nFd", "nPol", "#flags", "skipflags", "nwrite"]
    data_units = ["",         "MB",        "MB",     "",      "",    "",    "",     "",     "",       "", ""]
    data_types = [int,         float,      float, int,     int,   int,   int,    int,    int,      bool, int]
    dt = DTime(benchname=benchname, data_cols=data_cols, data_units=data_units, data_types=data_types, args=vars(args))
    if args.profile:
        dt.disable() # don't include startup in profiling

    f1 = dysh_data(example="getpslarge")
    trueNfiles, size_b, nflags = filestats(f1)
    nload =  trueNfiles
    size_mb = np.round(size_b, 2)
    if args.skipflags:
        nf = 0
    else:
        nf = nflags
    sdf = GBTFITSLoad(f1, skipflags=args.skipflags)#, nfiles=nfiles)
    s = sdf.stats()
    dt.tag("load", [s['nfiles'], size_mb, size_mb*s['nfiles'], s['nchan'], s['nrows'], s['ifnum'], s['fdnum'], s['plnum'], nf, args.skipflags, -1])
    print("Doing getps...")
    sb = sdf.getps(scan=np.arange(171,196),ifnum=0,plnum=0,fdnum=0) # All Sco-X OnOffs
    dt.tag("getps", [s['nfiles'], size_mb, size_mb*s['nfiles'], s['nchan'], s['nrows'], s['ifnum'], s['fdnum'], s['plnum'], nf, args.skipflags, len(sb)])
    print(f"Got Scanblock of length {len(sb)}")
    # scanblock write vs SDF write
    mkdir(f"{progname}_tmp",clean=True)
    if args.profile:
        dt.enable()
    print(f"Writing {args.writedata}...")
    for i in range(1, int(args.loop) + 1):
        if args.writedata == 'sdf':
            sdf.write(f"{progname}_tmp/sdfbench.fits",flags=args.skipflags)
            dt.tag(f"{args.writedata}write{i}", [s['nfiles'], size_mb, size_mb*s['nfiles'], s['nchan'], s['nrows'], s['ifnum'], s['fdnum'], s['plnum'], nf, args.skipflags, args.rows])
        else:
            for j in range(len(sb)):
                sb[0:j+1].write(f"{progname}_tmp/sbbench.fits",overwrite=True)
                dt.tag(f"{args.writedata}write{i}", [s['nfiles'], size_mb, size_mb*s['nfiles'], s['nchan'], s['nrows'], s['ifnum'], s['fdnum'], s['plnum'], nf, args.skipflags,j+1])

    dt.close()
    dt.report()
