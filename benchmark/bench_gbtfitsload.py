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


progname = "bench_gbtfitsload"
benchname = "GBTFITSLoad"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog=progname)
    parser.add_argument("--dobench", "-d", action="store_true", help="do the benchmark test")
    parser.add_argument("--key", "-k", action="store", help="input dysh_data key", default="multismallsmall")
    parser.add_argument(
        "--numfiles", "-n", action="store", help="number of SDFITS files to load for multifile data", default=1
    )
    parser.add_argument("--loop", "-l", action="store", help="number of times to loop", default=4)
    parser.add_argument("--skipflags", "-s", action="store_true", help="skip reading flags")
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
    parser.add_argument("--memory",      "-m", action="store_true",  help="track memory usage")
    parser.add_argument("--quit", "-q", action="store_true", help="quit early")
    parser.add_argument("--justtable", "-j", action="store_true", help="just print the existin table and exit")
    # parser.add_argument("--noindex",     "-n", action="store_true",  help="do not create dysh index table (pandas)")
    args = parser.parse_args()

    if args.quit:
        sys.exit(0)

    if args.justtable:
        if args.out is None:
            raise Exception("You must provide the table filename with -o FILENAME")
        table = Table.read(args.out, format="ascii.ecsv")
        table.pprint_all()
        sys.exit(0)

    # output table colnames, units, and dtypes
    # DTime automatically handles name and time, so just the additional columns go here.
    data_cols  = ["#files", "file_size", "totsize", "nchan", "nrow", "nIF", "nFd", "nPol", "#flags", "skipflags"]
    data_units = ["",         "MB",        "MB",     "",      "",    "",    "",     "",     "",       ""]
    data_types = [int,         float,      float, int,     int,   int,   int,    int,    int,      bool]
    dt = DTime(benchname=benchname, data_cols=data_cols, data_units=data_units, data_types=data_types, args=vars(args))

    f1 = dysh_data(accept=args.key)
    # use secret GBTFITSLOAD nfiles kwarg to limit number of files loaded.
    nfiles = int(args.numfiles)
    print(f"Loading not more than {nfiles} from {f1}")
    trueNfiles, size_b, nflags = filestats(f1)
    nload = min(nfiles, trueNfiles)
    size_mb = np.round(size_b, 2)
    print(f"Will load {nload} of {trueNfiles} files. FITS size per file {size_mb}MB, Flag lines {nflags}")
    if args.dobench:
        for i in range(1, int(args.loop) + 1):
            sdf = GBTFITSLoad(f1, skipflags=args.skipflags, nfiles=nfiles)
            s = sdf.stats()
            if args.skipflags:
                nf = 0
            else:
                nf = nflags
            dt.tag(f"load{i}", [s['nfiles'], size_mb, size_mb*s['nfiles'], s['nchan'], s['nrows'], s['ifnum'], s['fdnum'], s['plnum'], nf, args.skipflags])

    dt.close()
    dt.report()
