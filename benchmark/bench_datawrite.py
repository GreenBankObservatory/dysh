#!/usr/bin/env python
#
# Test writing of data write performance
import argparse
import os
import sys

import astropy.units as u
import numpy as np

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
        nsize = [os.stat(path).st_size]
        meandata = np.mean(nsize) * u.byte
        meanlines = 0 
        return (len(nsize), meandata.to(u.megabyte).value, meanlines)        
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
        if len(nlines) > 0:
            meanlines = int(np.mean(nlines))
        else:
            meanlines = 0            
        return (len(nsize), meandata.to(u.megabyte).value, meanlines)


progname = "bench_datawrite"
benchname = "DataWrite"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog=progname)
    # parser.add_argument("--key", "-k", action="store", help="input dysh_data key", default="getpslarge")
    # parser.add_argument(
    #    "--numfiles", "-n", action="store", help="number of SDFITS files to load for multifile data", default=1
    # )
    parser.add_argument("--loop", "-l",      action="store",      help="number of times to loop", default=4)
    parser.add_argument("--skipflags", "-s", action="store_true", help="skip reading/writing flags")
    parser.add_argument("--source", "-S",    action="store",      help="pick a source  1=NGC2415 2=NGC2782 3=Sco-X", default=3)
    parser.add_argument("--writedata",       action="store",      help="One of 'sdf' or 'sb' to write SDFITS or ScanBlock data", default="sdf")
    parser.add_argument("--rows", "-r",      action="store",      help="If `writedata` is sdf, the number of rows to write.")
    parser.add_argument("--out", "-o",       action="store",      help="output filename (astropy Table)", required=False)
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
    parser.add_argument(
        "--sortkey",
        "-x",
        action="store",
        help="How to sort the profiler statistics, 'cumulative' or 'time'",
        default="cumulative",
    )
    parser.add_argument("--memory", "-m", action="store_true", help="track memory usage")
    parser.add_argument("--quit", "-q", action="store_true", help="quit early")
    # parser.add_argument("--noindex",     "-n", action="store_true",  help="do not create dysh index table (pandas)")
    args = parser.parse_args()

    valid_write = ["sdf", "sb", "sp"]
    if args.writedata not in valid_write:
        raise ValueError(f"writedata must be one of {valid_write}")

    if args.writedata == "sdf" and args.rows is None:
        raise ValueError("You must supply number of rows to write (-r)  for SDFITS write test.")

    if args.quit:
        sys.exit(0)

    # output table colnames, units, and dtypes
    # DTime automatically handles name and time, so just the additional columns go here.
    data_cols = [
        "#files",
        "file_size",
        "totsize",
        "nchan",
        "nrow",
        "nIF",
        "nFd",
        "nPol",
        "#flags",
        "skipflags",
        "nwrite",
    ]
    data_units = ["", "MB", "MB", "", "", "", "", "", "", "", ""]
    data_types = [int, float, float, int, int, int, int, int, int, bool, int]
    dt = DTime(benchname=benchname, data_cols=data_cols, data_units=data_units, data_types=data_types, args=vars(args))
    if args.profile:
        dt.disable()  # don't include startup in profiling

    #f1 = dysh_data(example="getps2")   # for new dysh_data()
    f1 = dysh_data(example="getpslarge") # use 'getps2' for just source=1
    print("FILE:",f1)
    trueNfiles, size_b, nflags = filestats(f1)
    nload = trueNfiles
    size_mb = np.round(size_b, 2)
    if args.skipflags:
        nf = 0
    else:
        nf = nflags
    sdf = GBTFITSLoad(f1, skipflags=args.skipflags)  # , nfiles=nfiles)
    s = sdf.stats()
    dt.tag(
        "load",
        [
            s["nfiles"],
            size_mb,
            size_mb * s["nfiles"],
            s["nchan"],
            s["nrows"],
            s["ifnum"],
            s["fdnum"],
            s["plnum"],
            nf,
            args.skipflags,
            -1,
        ],
    )
    if   args.source == "1":
        scans = np.arange(152,154)   # NGC2415
    elif args.source == "2":
        scans = np.arange(156,160)   # NGC2782
    elif args.source == "3":
        scans = np.arange(171,197)   # Sco-X
    else:
        scans = np.arange(171,197)   # Sco-X
        
    print(f"Doing getps for source {args.source} scans={scans}")
    sb = sdf.getps(scan=scans, ifnum=0, plnum=0, fdnum=0)  # selected OnOffs
    dt.tag(
        "getps",
        [
            s["nfiles"],
            size_mb,
            size_mb * s["nfiles"],
            s["nchan"],
            s["nrows"],
            s["ifnum"],
            s["fdnum"],
            s["plnum"],
            nf,
            args.skipflags,
            len(sb),
        ],
    )
    print(f"Got Scanblock of length {len(sb)}")
    # scanblock write vs SDF write
    mkdir(f"{progname}_tmp", clean=True)
    if args.profile:
        dt.enable()
    print(f"Writing {args.writedata}...")
    for i in range(1, int(args.loop) + 1):
        if args.writedata == "sdf":
            sdf.write(f"{progname}_tmp/sdfbench.fits",  flags=args.skipflags, scan=scans, overwrite=True)
            dt.tag(
                f"{args.writedata}write{i}",
                [
                    s["nfiles"],
                    size_mb,
                    size_mb * s["nfiles"],
                    s["nchan"],
                    s["nrows"],
                    s["ifnum"],
                    s["fdnum"],
                    s["plnum"],
                    nf,
                    args.skipflags,
                    args.rows,
                ],
            )
        elif args.writedata == "sb":
            for j in range(len(sb)):
                sb[0 : j + 1].write(f"{progname}_tmp/sbbench.fits", overwrite=True)
                dt.tag(
                    f"{args.writedata}write{i}",
                    [
                        s["nfiles"],
                        size_mb,
                        size_mb * s["nfiles"],
                        s["nchan"],
                        s["nrows"],
                        s["ifnum"],
                        s["fdnum"],
                        s["plnum"],
                        nf,
                        args.skipflags,
                        j + 1,
                    ],
                )
        elif args.writedata == "sp":
            sp = sb.timeaverage()
            sp.write(f"{progname}_tmp/spbench.fits", overwrite=True)
            dt.tag(
                f"{args.writedata}write",
                [
                    s["nfiles"],
                    size_mb,
                    size_mb * s["nfiles"],
                    s["nchan"],
                    s["nrows"],
                    s["ifnum"],
                    s["fdnum"],
                    s["plnum"],
                    nf,
                    args.skipflags,
                    1,
                ],
            )

    dt.close()
    dt.report()
