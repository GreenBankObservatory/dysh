#!/usr/bin/env python 
from sdfitsload import SDFITSLoad, get_size, Obsblock,baseline
import os
import sys
import time
import cProfile
import pstats
from pstats import SortKey
import numpy as np
from astropy.table import Table,vstack
from specutils import Spectrum1D, SpectrumList,SpectralRegion
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="revisedstructure")
    parser.add_argument('--file','-f', action='store', help='input filename', required=True)
    parser.add_argument('--wcs',action='store_true', help='create WCS for each Spectrum1D')
    parser.add_argument('--fix',action='store_true', help='fix warnings about keywords in WCS creation')
    parser.add_argument('--profile', '-p',  action='store_true', help='run the profiler')
    parser.add_argument('--baseline','-b', action='store_true', help='remove baselines')
    parser.add_argument('--maxload','-m', action='store', help='maximum number of spectra to load (create obsblocks)',default=1E16)
#    parser.add_argument('--baseline','-b', action='store', help='remove baselines',default=None)
    args = parser.parse_args()
    print(f"using {args}")

    examples = "/data/gbt/examples/"
    files = [
    #"rxco-W/data/TSCAL_220105_W.raw.vegas/TSCAL_220105_W.raw.vegas.E.fits",
    #"onoff-L/data/TGBT21A_501_11.raw.vegas.fits", 
    #"nod-KFPA/data/TGBT22A_503_02.raw.vegas/TGBT22A_503_02.raw.vegas.F.fits",
    "misc/ngc5291.fits",
    #"misc/W3OH.fits",
    #"misc/IC1481.fits",
    #"mapping-L/data/TGBT17A_506_11.raw.vegas/TGBT17A_506_11.raw.vegas.A.fits", 
    #"mixed-fs-ps/data/AGBT16B_225_05/AGBT16B_225_05.raw.vegas/AGBT16B_225_05.raw.vegas.B.fits",

    ]
    timestr = ""
    files = [args.file]
    size = np.zeros(len(files))
    nhdu = np.zeros(len(files))
    nrows = np.zeros(len(files))
    nchanl = np.zeros(len(files))
    t0 = np.zeros(len(files))
    t0a = np.zeros(len(files))
    t1 = np.zeros(len(files))
    t2 = np.zeros(len(files))
    t3 = np.zeros(len(files))
    tb1 = np.zeros(len(files))
    tb2 = np.zeros(len(files))
    tb3 = np.zeros(len(files))
    t5 = np.zeros(len(files))
    hdu = np.zeros(len(files))
    i=0
    names =['File', 'Size','N_hdu','HDU','N_rows','N_chan', 'Load','Create_Obsblocks','Baseline_1', 'Baseline_2', 'Baseline_3', 'Total']
    units = ['', 'MB',   '',    '',    '',    '' ,      'ms',      'ms',          'ms',            'ms'      ,   'ms',          'ms' ]
    dtypes = [str, float, int,    int, int,   int,      float,     float,            float,          float,         float,        float]
    table = Table( names=names, meta={'name':'SDFITSLoad Timing'},
                   units=units,dtype=dtypes)
    for fn in files:
        basename = os.path.basename(fn)
        pr = cProfile.Profile()
        inf = f'{examples}{fn}'
        size[i] = os.path.getsize(fn)/1048576
        t0[i] = time.perf_counter_ns()
        pr.enable()
        s = SDFITSLoad(fn)
        t0a[i] = time.perf_counter_ns()
        for h in range(1,len(s._hdu)):
            hdu[i] = h
            t1[i] = time.perf_counter_ns()
            s._loadlists(fix=args.fix,wcs=args.wcs,hdu=int(hdu),maxspect=float(args.maxload))
            t2[i] = time.perf_counter_ns()
            nchanl[i] = len(s._obsblock[0][0].spectral_axis)
            if args.baseline:# is not None:
                j=0
                for o in s._obsblock:
                    p = o[0]
                    x = o[0].spectral_axis
                    nchan = len(x)
                    cdelt = x[1]-x[0]
                    center = x[nchan // 2]
                    width = 0.25*nchan*cdelt
                    exclude = SpectralRegion.from_center(center,width)
                    print("EXCLUDE ",exclude,file=sys.stderr)
                    baseline(speclist=o,order=1,exclude=exclude,maxspec=1000)
                    tb1[i] = time.perf_counter_ns()
                    print("Done O1",file=sys.stderr)
                    baseline(speclist=o,order=2,exclude=exclude,maxspec=1000)
                    tb2[i] = time.perf_counter_ns()
                    print("Done O2",file=sys.stderr)
                    baseline(speclist=o,order=3,exclude=exclude,maxspec=1000)
                    tb3[i] = time.perf_counter_ns()
                    print("Done O3",file=sys.stderr)
                    j=j+1
            t3[i] = time.perf_counter_ns()
            pr.disable()
            #s.summary()
            nhdu[i] = len(s._hdu)-1
            nrows[i] = s.nrows(h-1)
            load = (t0a-t0)/1E6
            obs  = (t2-t1)/1E6
            base1  = (tb1-t2)/1E6
            base2  = (tb2-tb1)/1E6
            base3  = (tb3-tb2)/1E6
            total = (t3-t0)/1E6
            table.add_row([basename,size,nhdu,hdu,nrows,nchanl,load,obs,base1,base2,base3,total])
        del s

        #s = io.StringIO()
        #sortby = (SortKey.TIME,SortKey.CUMULATIVE)
        ps = pstats.Stats(pr).sort_stats(SortKey.CUMULATIVE)
        ps.print_stats(25)
        i+=1


    for c in table.columns:
        if table[c].info.dtype == float:
            table[c].info.format = '0.1f'

    if os.path.exists("tab.out"):
        oldtable = Table.read("tab.out",format='ipac')
        table2 = vstack([oldtable,table])
    else:
        table2 = table
    table2.write("tab.out",format='ipac',overwrite=True)
    #t2.pprint_all()

