#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 16 11:06:02 2025

@author: teuben

there are 3 examples of OTF:
    
    otf1    L-band example (best developed, last example in this notebook)
    otf2    EDGE NGC0001 (first example)
    otf3    EDGE NGC5954 (second example)
"""


import os
import numpy as np
import matplotlib.pyplot as plt

import astropy.units as u
from astropy.io import fits
from specutils import Spectrum1D
from specutils.manipulation import box_smooth, gaussian_smooth, trapezoid_smooth


from dysh.fits.sdfitsload import SDFITSLoad
from dysh.fits.gbtfitsload import GBTFITSLoad
import dysh.util as util
from dysh.util.selection import Selection
from dysh.util.files import dysh_data

from dysh.spectra import ScanBlock

from astropy.io import fits
#%%
ks=['DATE-OBS','SCAN', 'IFNUM', 'PLNUM', 'FDNUM', 'INTNUM', 'CAL', 'PROCSEQN']


#  some more liberal panda dataframe display options
import pandas as pd
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 1000)
#%%  debugging

import dysh
dysh.log.init_logging(3)   # 0=ERROR 1=WARNING 2=INFO 3=DEBUG

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
            
            
def skyplot(sdf, scan=None, fdnum=8, size=None, title=None):
    """ for a selection of scans plot the skymap coverage for a given feed
    If size is given, plot will be rescaled in arcsec from -size to size
    in both RA and DEC around the mean position in X and Y (RA,DEC)
    
    @todo  change to using fig.,ax. instead of plt.
    @todo  allow sdf to be a file, so it can bypass the sdf.write() and
           grab sdf1 from the file

    @todo  add source and <ra,dec> to title
    @todo  add color bar with scan number as the scale
    """
    # make a temporary fits file with the selected scans
    sdf.write('skyplot.fits', scan=scan, overwrite=True, multifile=False)
    #
    sdf1 = GBTFITSLoad('skyplot.fits')
    f=sdf1._index['FDNUM']
    mask = f==fdnum
    if mask.sum() == 0:
        print(f"no data for scan={scan}, fdnum={fdnum}")
        return

    x=sdf1._index['CRVAL2'][mask]
    y=sdf1._index['CRVAL3'][mask]
    s=sdf1._index['SCAN'][mask]

    mz = (x!=0) | (y!=0)
    x = x[mz]
    y = y[mz]
    s = s[mz]
    
    if size is None:
        #  this will make RA run the wrong way, see heatmap
        plt.scatter(x, y, c=s)
        #return (x,y,s)
    else:
        xc = x.mean()
        yc = y.mean()
        dx = (x-xc)*np.cos(np.radians(yc)) * 3600
        dy = (y-yc) * 3600
        sc=plt.scatter(dx,dy,c=s)
        cbar = plt.colorbar(sc)
        cbar.set_label("scan")
        plt.xlim(size,-size)
        plt.ylim(-size,size)
        plt.xlabel('dRA (arcsec)')
        plt.ylabel('dDEC (arcsec)')
        plt.gca().set_aspect('equal')
        #return (xc,yc,s)                
             
    if title is not None:
        plt.title(title)             
    
def heatmap(sdf, cell=6, scans=None):
    kw=['SCAN', 'FDNUM', 'INTNUM', 'CRVAL2', 'CRVAL3']
    df=sdf._index[kw]
    if scans is not None:
        print(f"Selection on scans={scans}")
        df = df[df['SCAN'].isin(scans)]
    x=df['CRVAL2']
    y=df['CRVAL3']
    
    # it seems when EXPOSURE is small (0.018 in my example) the x and y are 0
    # note these are not VANE/SKY.
    
    mask = np.where(~np.isclose(x,0) & ~np.isclose(y,0))[0]
    nxy = len(x)
    n0 = nxy - len(mask)
    print(f'{n0}/{nxy} positions at (0,0)')
    x=df.iloc[mask]['CRVAL2']
    y=df.iloc[mask]['CRVAL3']
    #
    xc = x.mean()
    yc = y.mean()
    y = (y - yc)*3600.0
    x = (x - xc)*3600.0 * np.cos(np.pi*yc/180)
    xmin = x.min()
    xmax = x.max()
    ymin = y.min()
    ymax = y.max()
    bmax = max(-xmin, -ymin, xmax, ymax)
    print(f"bmax={bmax}")
    # @todo  fix edges to be symmetric around (0,0)
    xedges = np.arange(-bmax-cell, bmax+cell, cell)
    yedges = np.arange(-bmax-cell, bmax+cell, cell)
    print(xedges)
    
    heatmap, xedges, yedges = np.histogram2d(x, y, bins=(xedges, yedges))
    dx = xedges[1]-xedges[0]
    dy = yedges[1]-yedges[0]
    print("dx,dy",dx,dy,"arcsec")

    plt.figure()
    plt.imshow(heatmap.T, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], origin='lower', aspect='equal', cmap='viridis')
    plt.colorbar(label=f'Visits per  cell')

    plt.xlabel(f'dX, cell= {dx:.2f} arcsec')
    plt.ylabel(f'dY  cell= {dx:.2f} arcsec')
    plt.title(f'{sdf.filename}   {xc:.6f}  {yc:.6f}')

#%%    OTF2::    preparing ngc0001.fits since the data container is large  (1.3 GB)

f1 = dysh_data('AGBT21B_024_01')   
sdf1 = GBTFITSLoad(f1, skipflags=True)
# sdf1 = GBTFITSLoad(f1)    # ,skipflags=True)
# skipflags=False:   CPU times: user 9min 1s, sys: 3.49 s, total: 9min 5s Wall time: 8min 56s
# skipflags=True:    CPU times: user 4.07 s, sys: 539 ms, total: 4.61 s
sdf1.summary()


if False:
    scans=list(range(21,60))+list(range(65,104))
    sdf1.write('ngc0001.fits',scan=scans, overwrite=True, multifile=False)
    # 442MB
    sdf2 = GBTFITSLoad("ngc0001.fits")
    sdf2.summary()
    # CPU times: user 3.91 s, sys: 128 ms, total: 4.04 s
    mkdir("ngc0001")
    sdf1.write('ngc0001/file.fits',scan=scans, overwrite=True)
    sdf3 = GBTFITSLoad("ngc0001")
    sb1 = sdf3.gettp(scan=23,fdnum=0,ifnum=0,plnum=0,calibrate=True,cal=False)
    # there are 68 integrations in each scan
    scans=list(range(23,58)) + list(range(67,102))
    scans=list(range(23,58)) 
    mkdir("ngc0001sky")
    sdf1.write('ngc0001sky/file.fits',scan=scans, overwrite=True)
    sdf4 =  GBTFITSLoad("ngc0001sky")
    sdf4.summary()
    
    #
    sdf1.write("ngc0001_vane", scan=[21,22], overwrite=True, multifile=False)
    sdf5 = GBTFITSLoad("ngc0001_vane")
    sdf5.summary()
    
    
#%% skyplot NGC0001

skyplot(sdf2, scan=list(range(23,58)),  size=60, title="NGC0001 DecLatMap")
skyplot(sdf2, scan=list(range(67,76)),  size=60, title="NGC0001 RALongMap")

skyplot(sdf4,scan=list(range(23,58)),fdnum=0,title="NGC0001sky", size=60)

# UGC01659
skyplot(sdf1, scan=list(range(117,152)),  size=60, title="UGC01659 DecLatMap")
skyplot(sdf1, scan=list(range(161,196)),  size=60, title="UGC01659 RALongMap")
    
#%% work on ngc0001/file*fits

# this one can have the INTNUM bug   
sdf2 = GBTFITSLoad('ngc0001.fits', skipflags=True)    
sdf2.summary()


#%%  speed up for NGC0001 by usinng the smaller sdf4


f1 = dysh_data(accept='AGBT21B_024_01/AGBT21B_024_01.raw.vegas')
sdf1 = GBTFITSLoad(f1, skipflags=True)    # ~14sec

sdf1 = sdf4    # just 23..57


#%%  EDGE + waterfall

sdf = sdf4    # sdf1 is slow here, pick sdf4)

start = 23
end   = 57
#end  = 24
scans = list(range(start,end+1))
nchan = 1024
fdnum = sdf.udata("FDNUM")
#fdnum = [6,8,10,12]
nedge = 8
ncube = 0    # set to -1 if no waterfall needed

sb = ScanBlock()
sp1 = 0
waterfall = 0

for j,fd in enumerate(fdnum):    # 7 min in short map without waterfall;  but was 45mins
    print("Loop/Feed",j,fd)
    tpsb = sdf.gettp(scan=scans,fdnum=fd,ifnum=0,plnum=0,calibrate=True,cal=False)  # 6.5 sec
    nrows = tpsb[0].nrows
    nrows2 = nrows - 2  * nedge
    off = np.empty((nedge*2,nchan), dtype=float)
    exp_off = np.empty((nedge*2,nchan), dtype=float)
    on = np.empty((nrows-nedge*2,nchan), dtype=float)
    for i,(t,s) in enumerate(zip(tpsb,scans)):
        print(" loop/scan",j,i,s)
        off[:nedge] = t._calibrated[:nedge]
        off[nedge:] = t._calibrated[-nedge:]
        exp_off[:nedge] = t.exposure[:nedge,np.newaxis]
        exp_off[nedge:] = t.exposure[-nedge:,np.newaxis]
        on = t._calibrated[nedge:-nedge]
        off_avg = np.average(off, axis=0, weights=exp_off**2.)
        ta = (on - off_avg[np.newaxis,:])/(off_avg[np.newaxis,:])
        sb1 = sdf.gettp(scan=s,fdnum=fd,ifnum=0,plnum=0,calibrate=True,cal=False,intnum=list(range(nedge,nrows-nedge)))[0]  # 2.2 sec
        sb1._calibrated = ta
        
        if ncube >= 0:
            if ncube == 0:
                # set up waterfall cube
                nz = len(fdnum)               # feeds:    normally 16
                nz = 16
                ny = ta.shape[0] * len(scans) # time:     nscans*nint, should be nrows-2*nedge
                nx = ta.shape[1]              # channel:  nchan
                waterfall = np.zeros( (nz,ny,nx), dtype=float)
                ncube = 1
                sp1 = sb1.timeaverage()   # placeholder for baseline subtraction
            for k in range(nrows2):
                l = i*nrows2 + k
                sp1._data =  ta[k,:]
                #sp1.baseline(7, remove=True)
                #sp1.baseline(7,exclude=[(350,650)], remove=True)
                #sp1.baseline(7,include=[(10,350),(650,1000)], remove=True)
                #sp1._data = sp1._data - sp1._data.mean()
                sb1._calibrated[j] = sp1._data
                #sp1.undo_baseline()                # not needed
                waterfall[fd,l,:] = sp1._data
                
        sb.append(sb1)            
            
#%% write waterfall

hdu = fits.PrimaryHDU(waterfall*100)    # 100K, pending proper calibration using vane/sky
hdu.writeto('otf2-waterfall8.fits', overwrite=True)
        
#%% write calibrated spectra
# Save to SDFITS and then run the gbtgridder on it - this takes too long.....1h52m
sb.write("otf-test8.fits", overwrite=True)   # should this not INFO how much will be written?

# gbtgridder --size 32 32  --channels 10:1010 -o test3 --clobber --auto otf-test.fits

#%% bench : one feed  ~ 98sec
sb = ScanBlock()

fdnum = [7]

for fd in fdnum:
    print("Feed",fd)
    tpsb = sdf1.gettp(scan=scans,fdnum=fd,ifnum=0,plnum=0,calibrate=True,cal=False)  # 6.5 sec
    off = np.empty((4*2,nchan), dtype=float)
    exp_off = np.empty((4*2,nchan), dtype=float)
    on = np.empty((68-4*2,nchan), dtype=float)
    for i,(t,s) in enumerate(zip(tpsb,scans)):
        print("scan",s)
        off[:4] = t._calibrated[:4]
        off[4:] = t._calibrated[-4:]
        exp_off[:4] = t.exposure[:4,np.newaxis]
        exp_off[4:] = t.exposure[-4:,np.newaxis]
        on = t._calibrated[4:-4]
        off_avg = np.average(off, axis=0, weights=exp_off**2.)     # < 1ms
        ta = (on - off_avg[np.newaxis,:])/(off_avg[np.newaxis,:])  # < 1ms
        itpsb = sdf1.gettp(scan=s,fdnum=fd,ifnum=0,plnum=0,calibrate=True,cal=False,intnum=list(range(4,64)))[0]  # 2.2 sec
        itpsb._calibrated = ta
        sb.append(itpsb)


#   2.2sec for 60 integrations:    -> 36ms/intnum
#   calling gettp() doesn't make a difference if 

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






#%%    OTF3:: preparing ngc5954.fits since the data container is large (923 MB)

f1 = dysh_data('AGBT21B_024_20')  # AGBT21B_024_20.raw.vegas 
print(f1)
sdf1 = GBTFITSLoad(f1, skipflags=True)
# CPU times: user 8.14 s, sys: 778 ms, total: 8.92 s
sdf1 = GBTFITSLoad(f1)    # ,skipflags=True)
#    with skipflags:   CPU times: user 9min 1s, sys: 3.49 s, total: 9min 5s Wall time: 8min 56s
# without skipflags:   CPU times: user 4.07 s, sys: 539 ms, total: 4.61 s

sdf1.summary()
# 22..56 is the decmap
# There are 156 rows, of which 35 (22%) are the decmap

# NGC0001 is scans 

if False:
    scans=list(range(20,59))+list(range(64,103))    # this includes the matching sky/vane
    mkdir("ngc5954")
    sdf1.write('ngc5954/file.fits',scan=scans, overwrite=True)
    sdf3 = GBTFITSLoad("ngc5954")
    sb1 = sdf3.gettp(scan=23,fdnum=0,ifnum=0,plnum=0,calibrate=True,cal=False)
    # there are 68 integrations in here
    

    
#%%    DecLatMap:   bad slow way]

#   DecLatMap:    sky/vane  20,21  and 57,58    data:  22,56
#   one row for 60 points and 16 beams took almost 6' this way
   
start = 22
end   = 56
end   = 22
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

#%%


#   skyplot(sdf3, scan=list(range(22,57)),  size=60, title="NGC5954 DecLatMap")
#   skyplot(sdf3, scan=list(range(66,101)), size=60, title="NGC5954 RALongMap")


#%%
sdf3.write('junk1.fits', scan=list(range(22,57)), overwrite=True, multifile=False)
sdf3a = GBTFITSLoad('junk1.fits')
sdf3a.summary()

#%% collect

sdf = sdf3a

f=sdf._index['FDNUM']
mask = f==8

x=sdf._index['CRVAL2'][mask]
y=sdf._index['CRVAL3'][mask]
s=sdf._index['SCAN'][mask]

mz = (x!=0) | (y!=0)

#%%    DecLatMap sky coverage

# need:   SCAN,FDNUM,CRVAL2,CRVAL3,INTNUM
kw=['SCAN', 'FDNUM', 'INTNUM', 'CRVAL2', 'CRVAL3']
    
#   DecLatMap:    sky/vane  20,21  and 57,58    data:  22,56
#   one row for 60 points and 16 beams took almost 6'
   
start = 22
end   = 56
end   = 22

sdf = sdf3a
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
    

#%%


#%%

heatmap(sdf4)

#%%


f1 = dysh_data('AGBT21B_024_41')  # AGBT21B_024_20.raw.vegas 
print(f1)
sdf1 = GBTFITSLoad(f1, skipflags=True)  #

# 	UGC03960   31..65  83..117
scans = list(range(31,66)) + list(range(83,118))
print(scans)
heatmap(sdf1, cell=6, scans=scans)

#%%  testps with getsigref

f1 = dysh_data(test='getps')
sdf1 = GBTFITSLoad(f1)

sp1 = sdf1.getps(fdnum=0, ifnum=0, plnum=0).timeaverage()

sp2 = sdf1.getsigref(152,153,fdnum=0, ifnum=0, plnum=0).timeaverage()


sp1.stats()         # rms ~ 0.67
(sp1-sp2).stats()   # rms ~ 2e-8


#%%  OTF1::  L-band

__help__ = """
    SCAN         OBJECT VELOCITY       PROC  ...  # INT # FEED     AZIMUTH   ELEVATIO
0      6          3C286      0.0      OnOff  ...     13      1  248.365704   72.55312
1      7          3C286      0.0      OnOff  ...     13      1  250.038525  73.677706
2      8         SgrB2M     57.0      Track  ...     61      1  142.384473  12.254537
3      9           W33B     64.0      Track  ...     61      1  132.205692   17.97864
4     10  G30.589-0.044     38.0      Track  ...     61      1  115.481079  25.552326
5     11  G31.412+0.307     97.0      Track  ...     61      1  115.803415   27.09864
6     12  G35.577-0.029     49.0      Track  ...     61      1   112.31392  29.033662
7     13  G40.622-0.137     32.0      Track  ...     61      1   107.77357  31.315933
8     14        NGC6946     45.0  DecLatMap  ...     61      1   38.962631  41.364443
9     15        NGC6946     45.0  DecLatMap  ...     61      1   38.991018   41.44875
10    16        NGC6946     45.0  DecLatMap  ...     61      1   38.991814  41.535999
11    17        NGC6946     45.0  DecLatMap  ...     61      1   39.019325   41.62031
12    18        NGC6946     45.0  DecLatMap  ...     61      1   39.019752  41.707519
13    19        NGC6946     45.0  DecLatMap  ...     61      1   39.046779  41.791907
14    20        NGC6946     45.0  DecLatMap  ...     61      1   39.046452  41.879147
15    21        NGC6946     45.0  DecLatMap  ...     61      1   39.072974  41.963652
16    22        NGC6946     45.0  DecLatMap  ...     61      1   39.072156  42.050849
17    23        NGC6946     45.0  DecLatMap  ...     61      1   39.098252  42.135459
18    24        NGC6946     45.0  DecLatMap  ...     61      1   39.096777  42.222654
19    25        NGC6946     45.0  DecLatMap  ...     61      1   39.122446  42.307343
20    26        NGC6946     45.0  DecLatMap  ...     61      1    39.12023  42.394543
21    27        NGC6946     45.0      Track  ...     61      1   39.088117  42.149313

[22 rows x 13 columns]

"""

f1 = dysh_data(example='mapping-L/data/TGBT17A_506_11.raw.vegas')    # 2.3GB
print(f1)

sdf1 = GBTFITSLoad(f1)   #   ~  5 secs
sdf1.summary()


#%%  testing one scan
# note INTNUM has odd values for scan=6

if False:
    mkdir("otf1")
    scans=[20,27]
    sdf1.write('otf1/file.fits',scan=scans, overwrite=True, ifnum=0, plnum=0)

#  scans 14..26 is the DecLatMap  (procseq=1..13)
#  scan 27 is a Track on the reference position

# ifnum: 0..4    
# fdnum: 0
# plnum: 0..1

sdf1.getsigref(scan=26, ref=27, fdnum=0, ifnum=0, plnum=0)
# ValueError: cannot reindex on an axis with duplicate labels
# error gone when fix 425 was merged


#%%



#%%   prepare a smaller sdfits just for NGC6946

ifnum=0        # needs to be 0
plnum=0        # pick 0 and 1 to average
fdnum=0        # fixed
nint = 61      # nunber of integrations per scan  (a.k.a. nrows/2)
intnum = 30    # pick something in the middle for a test
sig = 20
ref = 27

# make smaller sdfits
scans = list(range(14,28))
mkdir("otf1")
# scans = [sig, ref]
sdf1.write('otf1/file.fits',scan=scans, overwrite=True, fdnum=fdnum, ifnum=ifnum, plnum=plnum)

                     
sdf2 = GBTFITSLoad('otf1')
sdf2.summary()

# skyplot(sdf2, scan=list(range(14,27)),  fdnum=0, size=1500, title="NGC6946 DecLatMap")

#  quick test
sp2 = sdf2.getsigref(scan=sig, ref=ref, fdnum=fdnum, ifnum=ifnum, plnum=plnum, intnum=intnum).timeaverage()
sp2.plot(xaxis_unit="km/s")

# region between -2000 and -500   and 500..2000 can be used for baseline

#%% pick sdf1 or sdf2 for speed testing the calibration portion

# on the big one it takes ~ 1'12"
# small one it takes ~      6"

scans = list(range(14,27))
nscan = len(scans)
fdnum = 0
ifnum = 0
plnum = 0

# test to get nint, nchan
sb = sdf1.getsigref(scans[0], ref, fdnum=fdnum, ifnum=ifnum, plnum=plnum)[0]
nint = sb.nint - 1     # last integration always problematic
nchan = sb.nchan
print(f"Using nchan={nchan} and nint={nint}")

# sdf = sdf1    # big one (2300MB)
sdf = sdf2    # just the relevant scans 14..27  (5MB)

ny = nscan*nint
nx = nchan
waterfall = np.zeros( (ny,nx), dtype=float)
sp1 = 0

#  mode=0,1   ~5sec    - no real baseline subtraction
#  mode=2     ~60 sec  - baseline subtraction
mode = 2

sb = ScanBlock()

print("Using", sdf.filename, "for", scans)

# bmask for mode=1
bmask = np.arange(nchan)
bmask = (bmask>500) & (bmask<1500) | (bmask>2500) & (bmask<3500)
for i,s in enumerate(scans):
    print(i,s)
    sb1 =  sdf.getsigref(scan=s, ref=ref, fdnum=fdnum, ifnum=ifnum, plnum=plnum)[0]
    if sp1 == 0:
        sp1 = sb1.timeaverage()   # now spectrum is in sp1.flux.data
    for j in range(nint):
        print("Loop",i,j,s)
        k = i*nint + j
        waterfall[k] = sb1._calibrated[j]
        if mode == 1:
            # poor man's baseline reduction
            smean = sb1._calibrated[j][bmask].mean()
            waterfall[k] = sb1._calibrated[j] - smean
        elif mode == 2:
            # sp1 = sb1.timeaverage()
            #print("DATA0", sb1._calibrated[j][100])
            sp1._data = sb1._calibrated[j]     # _data    value?
            sp1.baseline(1, include=[(500, 1500), (2500, 3550)], remove=True, model='polynomial')
            #print(sp1.baseline_model)
            sb1._calibrated[j] = sp1._data
            sp1.undo_baseline()
            #sp1.hack()
            waterfall[k] = sb1._calibrated[j]
            #print("DATA1", sb1._calibrated[j][100])
    sb.append(sb1)    
    

plt.figure(1)
plt.clf()
plt.imshow(waterfall, vmin=-0.25, vmax=0.25, origin="lower")
plt.xlabel("Channel")
plt.ylabel("Scan/Time")
plt.colorbar(label="Ta*")


#%%

hdu = fits.PrimaryHDU(waterfall)
hdu.writeto('otf1-waterfall.fits', overwrite=True)
    
sb.timeaverage().plot(xaxis_unit='chan',title='one')
    
#%%   subtract a good baseline model from all scanblocks
ta = sb.timeaverage()
ta.plot(xaxis_unit='chan', title='two')
ta.stats(qac=True)

sb.timeaverage().plot(xaxis_unit='chan',title='one again')

if False:
    ta.baseline(1, include=[(500, 1500), (2500, 3550)], remove=True, model='polynomial')
    sb.subtract_baseline(ta.baseline_model, tol=1000, force=True)
else:
    # can also do -500 to 500 km/s
    ta.baseline(1, [(1500,2500)], remove=True, model='polynomial')
    print(np.mean( sb[0]._calibrated[0]))
    sb.subtract_baseline(ta.baseline_model, force=True)
    print(np.mean( sb[0]._calibrated[0]))
    ta2 = sb.timeaverage()
    ta2.plot(xaxis_unit='chan', title='three')
    ta2.stats(qac=True)

#%% write calibrated spectra
sb.write("otf-test2.fits", overwrite=True)    #  300 ms




#%%   gbtgridder example - if you have gbtgridder

#
# gbtgridder --size 32 32  --channels 500:3500 -o test2 --clobber --auto otf-test2.fits
# pixels:  2.9'  100x100 ->  
# confirmed the two cubes from sdf1 and sdf2 are identical 

