#!/usr/bin/env python 
from sdfitsload import SDFITSLoad, get_size, Obsblock,baseline
import os
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
    parser.add_argument('--file','-f', action='store', help='input filename')
    parser.add_argument('--wcs',action='store_true', help='create WCS for each Spectrum1D')
    parser.add_argument('--fix',action='store_true', help='fix warnings about keywords in WCS creation')
    parser.add_argument('--profile',action='store_true', help='run the profiler')
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
    t0 = np.zeros(len(files))
    t1 = np.zeros(len(files))
    t2 = np.zeros(len(files))
    t3 = np.zeros(len(files))
    i=0
    for fn in files:
        pr = cProfile.Profile()
        inf = f'{examples}{fn}'
        size[i] = os.path.getsize(fn)/1048576
        t0[i] = time.perf_counter_ns()
        pr.enable()
        s = SDFITSLoad(fn)
        t1[i] = time.perf_counter_ns()
        s._loadlists(fix=args.fix,wcs=args.wcs)
        t2[i] = time.perf_counter_ns()
        for o in s._obsblock:
            p = o[0]
            x = o[0].spectral_axis
            nchan = len(x)
            cdelt = x[1]-x[0]
            center = x[nchan // 2]
            width = 0.25*nchan*cdelt
            exclude = SpectralRegion.from_center(center,width)
            print("EXCLUDE ",exclude)
            baseline(speclist=o,order=1,exclude=exclude)
        t3[i] = time.perf_counter_ns()
        pr.disable()
        #s.summary()
        nhdu[i] = len(s._hdu)
        nrows[i] = np.sum(s._nrows)
        del s

        #s = io.StringIO()
        #sortby = (SortKey.TIME,SortKey.CUMULATIVE)
        ps = pstats.Stats(pr).sort_stats(SortKey.CUMULATIVE)
        ps.print_stats(20)
        i+=1

#        timestr+=f'   {size:.0f}\t{nhdu}\t{nrows}\t{(t1-t0)/1E6:.1f}\t{(t2-t1)/1E6:.1f}\t{(t2-t0)/1E6:.1f}\n'
    
    load = (t1-t0)/1E6
    obs  = (t2-t1)/1E6
    base  = (t3-t2)/1E6
    total = (t3-t0)/1E6
    names =['Size','N_hdu','N_rows','Load','Create_Obsblocks','Baseline', 'Total']
    units = ['MB',   '',      '',    'ms',      'ms',          'ms',       'ms']
    dtypes = [float, int, int, float,float,float,float]
    table = Table(data=[size,nhdu,nrows,load,obs,base,total],
                   names=names, meta={'name':'SDFITSLoad Timing'},
                   units=units,dtype=dtypes)
    for c in table.columns:
        if table[c].info.dtype == float:
            table[c].info.format = '0.1f'

    if os.path.exists("tab.out"):
        oldtable = Table.read("tab.out",format='ipac')
        t2 = vstack([oldtable,table])
    else:
        t2 = table
    t2.write("tab.out",format='ipac',overwrite=True)
    #t2.pprint_all()

