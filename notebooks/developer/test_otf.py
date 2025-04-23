#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 16 11:06:02 2025

@author: teuben
"""


import os
import numpy as np
import matplotlib.pyplot as plt

import astropy.units as u
from astropy.io import fits
from specutils import Spectrum1D
from specutils.manipulation import box_smooth, gaussian_smooth, trapezoid_smooth

import matplotlib.pyplot as plt

from dysh.fits.sdfitsload import SDFITSLoad
from dysh.fits.gbtfitsload import GBTFITSLoad
import dysh.util as util
from dysh.util.selection import Selection
from dysh.util.files import dysh_data

#%%  helper functions

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
            

#%%    preparing ngc0001.fits since the data container is large  (1.3 GB)

f1 = dysh_data('AGBT21B_024_01')   
sdf1 = GBTFITSLoad(f1, skipflags=True)
sdf1 = GBTFITSLoad(f1)    # ,skipflags=True)
#    with skipflags:   CPU times: user 9min 1s, sys: 3.49 s, total: 9min 5s Wall time: 8min 56s
# without skipflags:   CPU times: user 4.07 s, sys: 539 ms, total: 4.61 s



# NGC0001 is scans 

if False:
    scans=list(range(21,60))+list(range(65,104))
    sdf1.write('ngc0001.fits',scan=scans, overwrite=True, multifile=False)
    # 442MB
    sdf2 = GBTFITSLoad("ngc0001.fits")
    # CPU times: user 3.91 s, sys: 128 ms, total: 4.04 s
    mkdir("ngc0001")
    sdf1.write('ngc0001/file.fits',scan=scans, overwrite=True)
    sdf3 = GBTFITSLoad("ngc0001")
    sb1 = sdf3.gettp(scan=23,fdnum=0,ifnum=0,plnum=0,calibrate=True,cal=False)
    # there are 68 integrations in here
    
#%%
   
start = 23
end   = 57
noff = 4  

sdf = sdf3
tsys = 200

fp = open("ngc0001.txt","w")

for s in list(range(start,end+1)):
    sb1 = sdf3.gettp(scan=s,fdnum=0,ifnum=0,plnum=0,calibrate=True,cal=False)[0]
    nrows = sb1.nrows
    int_0 = list(range(0,noff)) + list(range(nrows-noff,nrows))
    off = sdf.gettp(scan=s,fdnum=0,ifnum=0,plnum=0,intnum=int_0,calibrate=True,cal=False)[0].timeaverage()
    int_1 = list(range(noff,nrows-noff))
    for i in int_1:
        on = sdf.gettp(scan=s,fdnum=0,ifnum=0,plnum=0,intnum=i,calibrate=True,cal=False)[0].timeaverage()
        ta = tsys * (on-off)/off
        print(s,i)
        fp.write("%d %d %s\n" % (s,i,ta.stats(qac=True)))
fp.close()
        
#  CPU times: user 1h 23min 47s, sys: 5min 18s, total: 1h 29min 6s    Wall time: 1h 29min 4s  




#%%    preparing ngc5954.fits since the data container is large (923 MB)

f1 = dysh_data('AGBT21B_024_20')  # AGBT21B_024_20.raw.vegas 
print(f1)
sdf1 = GBTFITSLoad(f1, skipflags=True)
# CPU times: user 8.14 s, sys: 778 ms, total: 8.92 s
sdf1 = GBTFITSLoad(f1)    # ,skipflags=True)
#    with skipflags:   CPU times: user 9min 1s, sys: 3.49 s, total: 9min 5s Wall time: 8min 56s
# without skipflags:   CPU times: user 4.07 s, sys: 539 ms, total: 4.61 s



# NGC0001 is scans 

if False:
    scans=list(range(20,59))+list(range(64,103))
    mkdir("ngc5954")
    sdf1.write('ngc5954/file.fits',scan=scans, overwrite=True)
    sdf3 = GBTFITSLoad("ngc5954")
    sb1 = sdf3.gettp(scan=23,fdnum=0,ifnum=0,plnum=0,calibrate=True,cal=False)
    # there are 68 integrations in here
    
#%%    DecLatMap

#   DecLatMap:    sky/vane  20,21  and 57,58    data:  22,56
#   one row for 60 points and 16 beams took almost 6'
   
start = 22
end   = 56
#end   = 22
noff = 4  

sdf = sdf3
tsys = 200

fp = open("ngc5954.txt","w")

for s in list(range(start,end+1)):
    for f in range(16):
        sb1 = sdf3.gettp(scan=s,fdnum=f,ifnum=0,plnum=0,calibrate=True,cal=False)[0]
        nrows = sb1.nrows
        nrows = nrows - 1
        int_0 = list(range(0,noff)) + list(range(nrows-noff,nrows))
        int_1 = list(range(noff,nrows-noff))
        off = sdf.gettp(scan=s,fdnum=f,ifnum=0,plnum=0,intnum=int_0,calibrate=True,cal=False)[0].timeaverage()
        sb1 = sdf.gettp(scan=s,fdnum=f,ifnum=0,plnum=0,intnum=int_1,calibrate=True,cal=False)[0]
        for i in range(len(int_1)):
            on = sb1.calibrated(i)
            ta = tsys * (on-off)/off
            c2 = on.meta["CRVAL2"]
            c3 = on.meta["CRVAL3"]
            print(f,s,i)
            fp.write("%d %d %d %g %g %s\n" % (f,s,i,c2,c3,ta.stats(qac=True)))
fp.close()
        
#  CPU times: user 1h 23min 47s, sys: 5min 18s, total: 1h 29min 6s    Wall time: 1h 29min 4s  


#%%    DecLatMap sky coverage

# need:   SCAN,FDNUM,CRVAL2,CRVAL3,INTNUM
kw=['SCAN', 'FDNUM', 'INTNUM', 'CRVAL2', 'CRVAL3']
    
#   DecLatMap:    sky/vane  20,21  and 57,58    data:  22,56
#   one row for 60 points and 16 beams took almost 6'
   
start = 22
end   = 56
end   = 22

sdf = sdf3
tsys = 200

fp = open("coverage.txt","w")

for s in list(range(start,end+1)):
    for f in range(16):
        sb1 = sdf3.gettp(scan=s,fdnum=f,ifnum=0,plnum=0,calibrate=True,cal=False)[0]
        nrows = sb1.nrows
        print(nrows)
        for i in range(nrows):
            sp = sb1.calibrated(i)
            c2 = sp.meta["CRVAL2"]
            c3 = sp.meta["CRVAL3"]
            print(f,s,i, c2,c3)
            fp.write("%d %d %d %g %g %s\n" % (f,s,i,c2,c3,ta.stats(qac=True)))
fp.close()
        
#%%  DecLatMap select an RA,DEC to average the spectrum

#%%   counting sizes


def sizes(sdf):
    s = {}
    s['nfiles'] = len(sdf.filenames())
    s['nrows'] = len(sdf._index)
    s['nchan'] = len(sdf.getspec(0).flux)
    s['scan'] = len(sdf._index['SCAN'].unique())
    s['fdnum'] = len(sdf._index['FDNUM'].unique())
    s['ifnum'] = len(sdf._index['IFNUM'].unique())
    s['plnum'] = len(sdf._index['PLNUM'].unique())
    s['intnum'] = len(sdf._index['INTNUM'].unique())
    s['sig'] = len(sdf._index['SIG'].unique())
    s['cal'] = len(sdf._index['CAL'].unique())
    return s
    