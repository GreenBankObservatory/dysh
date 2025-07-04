#!/usr/bin/env python
#
# Test performance of OTF calibration

import argparse
import os
import sys

from dysh.fits.gbtfitsload import GBTFITSLoad
from dysh.spectra import ScanBlock
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


benchname = "otf"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog=sys.argv[0])
    # fmt: off
    # parser.add_argument("--file",        "-f", action="store",       help="input filename", required=True)
    parser.add_argument("--dobench",     "-d", action="store_true",  help="do the benchmark test")
    parser.add_argument("--key",         "-k", action="store",       help="input dysh_data key", default="test1")
    parser.add_argument("--timeaverage", "-t", action="store_true",  help="time average the Scanblocks to make a Spectrum")
    parser.add_argument("--nocalibrate", "-n", action="store_true",  help="DON'T calibrate the data", default=False)
    parser.add_argument("--loop",        "-l", action="store",       help="number of times to loop", default=4)
    parser.add_argument("--feeds",       "-f", action="store",       help="number of feeds to use (1..16)", default=16)
    parser.add_argument("--skipflags",   "-s", action="store_true",  help="skip reading flags")
    parser.add_argument("--big",         "-b", action="store_true",  help="use the original big file")

    parser.add_argument("--out",         "-o", action="store",       help="output filename (astropy Table)", required=False)
    parser.add_argument("--append",      "-a", action="store_true",  help="append to previous output file (astropy Table)", required=False)
    parser.add_argument("--overwrite",   "-w", action="store_true",  help="overwrite a previous output file (astropy Table)", required=False)
    parser.add_argument("--profile",     "-p", action="store_true",  help="run the profiler")
    parser.add_argument("--statslines",  "-e", action="store",       help="number of profiler statistics lines to print", default=25)
    parser.add_argument("--sortkey",     "-x", action="store",       help="How to sort the profiler statistics, 'cumulative' or 'time'", default="cumulative")
    parser.add_argument("--memory",      "-m", action="store_true",  help="track memory usage")
    parser.add_argument("--quit",        "-q", action="store_true",  help="quit early")
    #parser.add_argument("--index", "-i", action="store_true", help="create dysh index table (pandas)")
    # fmt: on
    args = parser.parse_args()
    print(f"using {args}")

    if args.quit:
        sys.exit(0)

    if args.nocalibrate and args.timeaverage:
        raise Exception("You must calibrate if you want to time average")

    data_cols = ["skipflags"]
    data_units = [""]
    data_types = [str]
    dt = DTime(benchname=benchname, data_cols=data_cols, data_units=data_units, data_types=data_types, args=vars(args))

    sk = str(args.skipflags)

    do_L = True  # hardcoded for now
    if do_L:
        # NGC6946 in L-band
        f1 = dysh_data(example="mapping-L/data/TGBT17A_506_11.raw.vegas")
        scans = list(range(14, 28))
    else:
        # NGC5954 double galaxy in EDGE survey
        # too slow?
        f1 = dysh_data(accept="AGBT21B_024_20/AGBT21B_024_20.raw.vegas")  # @todo why does just AGBT21B_024_20  not work
        scans = list(range(22, 57))  # just one DecLatMap of 35 scans, no Vane/Sky
    print("Loading ", f1)

    sdf1 = GBTFITSLoad(f1, skipflags=args.skipflags)
    print("STATS:", sdf1.stats())
    dt.tag("load1", [sk])
    if args.big:
        sdf2 = sdf1
        dt.tag("write1", [sk])
        dt.tag("load2", [sk])
    else:  # need to create smaller file first, 13 "on" scans, plus a single "off"
        if do_L:
            mkdir("ngc6946")
            sdf1.write("ngc6946/file.fits", scan=scans, overwrite=True, fdnum=0, ifnum=0, plnum=0)

            dt.tag("write1", [sk])
            sdf2 = GBTFITSLoad("ngc6946", skipflags=args.skipflags)
            print("STATS:", sdf2.stats())
            dt.tag("load2", [sk])
            del sdf1
        else:
            mkdir("ngc5954a")
            sdf1.write("ngc5954a/file.fits", scan=scans, overwrite=True)
            dt.tag("write1", [sk])
            sdf2 = GBTFITSLoad("ngc5954a", skipflags=args.skipflags)
            print("STATS:", sdf2.stats())
            dt.tag("load2", [sk])
            del sdf1
        # sdf1 = np.arange(100)   # dummy space
        # dt.tag("mem",[sk])

    # note this loop is simpler than the more realistic one in test_otf.py
    if args.dobench:
        if do_L:
            sb = ScanBlock()
            for s in scans[:-1]:
                sb1 = sdf2.getsigref(scan=s, ref=scans[-1], fdnum=0, ifnum=0, plnum=0)[0]
                sb.append(sb1)
            dt.tag("getsigref", [sk])
            sb.write("otf.fits", overwrite=True)  #  14388480 bytes,   840 spectra (should be 793 per pol)
            dt.tag("write2", [sk])

        else:
            calibrate = not args.nocalibrate
            intnums = list(range(4, 64))
            intnums = [0]
            scan = [22]
            for i in range(1, int(args.loop) + 1):
                for f in range(int(args.feeds)):  # loop over all feeds
                    sb = sdf2.gettp(scan=scan, fdnum=f, ifnum=0, plnum=0, intnum=intnums, calibrate=True, cal=False)
                dt.tag(f"gettp{i}s", [sk])

    dt.tag("report", [sk])
    dt.close()
    dt.report()

    print("final", dt.total() / 1000, " sec")
