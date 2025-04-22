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


benchname = "otf"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog=sys.argv[0])
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

    data_cols  = ["skipflags"]
    data_units = [""]
    data_types = [str]
    dt = DTime(benchname=benchname,
               out=args.out, append=args.append, overwrite=args.overwrite, profile=args.profile, statslines=args.statslines, # @todo   use args ?
               data_cols=data_cols, data_units=data_units, data_types=data_types)  # no data=[] supported yet

    sk=str(args.skipflags)

    f1 = dysh_data(accept='AGBT21B_024_20/AGBT21B_024_20.raw.vegas')  # @todo why does just AGBT21B_024_20  not work
    print("Loading ",f1)
    sdf1 = GBTFITSLoad(f1, skipflags=args.skipflags)
    dt.tag("load1", [sk])
    if args.big:
        sdf2 = sdf1
        dt.tag("write1",[sk])
        dt.tag("load2",[sk])        
    else:
        mkdir('ngc5954a')
        scans = list(range(22,57))  # just one DecLatMap of 35 scans, no Vane/Sky
        sdf1.write('ngc5954a/file.fits', scan=scans,overwrite=True)
        dt.tag("write1",[sk])
        sdf2 = GBTFITSLoad('ngc5954a', skipflags=args.skipflags)
        dt.tag("load2",[sk])
        
    calibrate = not args.nocalibrate
    if args.dobench:
        scan = [22]
        for i in range(1,int(args.loop)+1):
            for f in range(int(args.feeds)):    # loop over all feeds
                sb = sdf2.gettp(scan=scan, fdnum=f, ifnum=0, plnum=0, intnum=0, calibrate=True, cal=False)
            dt.tag(f"gettp{i}s", [sk])

    dt.tag('report',[sk])
    dt.close()
    dt.report()

    print('final',dt.total()/1000,' sec')
    

__result__ = """

# -s 
 load1 12059.1      True
write1  8107.9      True   varies to 9154
 load2  3003.0      True   varies to 3156
report     0.0      True

# -s -d -f 16 -l 1 -b
gettp1s 65762.6      True

# -s -d -f 1 -l 4 -b
gettp1s  4306.2      True
gettp2s  4220.9      True
gettp3s  4166.5      True
gettp4s  4254.6      True

# -s -d -f 1 -l 4
gettp1s   972.9      True
gettp2s   983.8      True
gettp3s  1156.3      True
gettp4s   985.4      True


-b shows that reading small chunks from the big file is expensive, in this case 4s vs. 1s
which is roughly the ratio of rows



"""
