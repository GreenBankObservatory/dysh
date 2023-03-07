from sdfitsload import SDFITSLoad, get_size
import os
import time
import cProfile
import pstats
from pstats import SortKey
import numpy as np

if __name__ == "__main__":
    examples = "/data/gbt/examples/"
    files = [
    "mapping-L/data/TGBT17A_506_11.raw.vegas/TGBT17A_506_11.raw.vegas.A.fits", 
    "onoff-L/data/TGBT21A_501_11.raw.vegas.fits", 
    "rxco-W/data/TSCAL_220105_W.raw.vegas/TSCAL_220105_W.raw.vegas.A.fits", 
    "rxco-W/data/TSCAL_220105_W.raw.vegas/TSCAL_220105_W.raw.vegas.E.fits",
    "misc/ngc5291.fits",
    "misc/W3OH.fits",
    "misc/IC1481.fits",
    ]
    timestr = ""
    for fn in files:#[3:]:
        pr = cProfile.Profile()
        inf = f'{examples}{fn}'
        size = os.path.getsize(inf)/1048576
        #size = get_size(inf)/1048576
        #print(f'\n{inf} size: {size:.1f} MB')
        t0 = time.perf_counter_ns()
        pr.enable()
        s = SDFITSLoad(inf)
        t1 = time.perf_counter_ns()
        s._loadlists(fix=True)
        pr.disable()
        t2 = time.perf_counter_ns()
        #s.summary()
        nhdu = len(s._hdu)
        nrows = np.sum(s._nrows)

        #s = io.StringIO()
        sortby = SortKey.CUMULATIVE
        ps = pstats.Stats(pr).sort_stats(sortby)
        ps.print_stats(10)

        timestr+=f'   {size:.0f}\t{nhdu}\t{nrows}\t{(t1-t0)/1E6:.1f}\t{(t2-t1)/1E6:.1f}\t{(t2-t0)/1E6:.1f}\n'
        #del s
    print('\n# Timing (ms)')
    print('# Size\tNhdu\tNrows\tLoad\tCreate Obsblocks\t\tTotal')
    print(timestr)

    #print(s.getvalue())

