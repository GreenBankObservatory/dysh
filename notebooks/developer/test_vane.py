#!/usr/bin/env python3

"""
This script was developed in spyder during the vane calibration work.

From nodding we have two cases:
nod1  fdnum=[10,1]      unbalanced calrows 0 != 6        tp nocal;  TBD
nod3  fdnum=[0,1]       unbalanced calrows 0 != 6        tp nocal;  TBD
"""

import os
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from scipy.stats import norm
from astropy.io import fits
import astropy.units as u

from dysh.fits.sdfitsload import SDFITSLoad
from dysh.fits.gbtfitsload import GBTFITSLoad
from dysh.util.files import dysh_data
from dysh.util.selection import Selection
from dysh.spectra.core import mean_tsys
    
#  useful keys for a mult-beam observation listing

k=['DATE-OBS','SCAN', 'IFNUM', 'PLNUM', 'FDNUM', 'INTNUM', 'PROCSCAN','FEED', 'SRFEED', 'FEEDXOFF', 'FEEDEOFF', 'SIG', 'CAL', 'PROCSEQN', 'PROCSIZE']
ks=['DATE-OBS','SCAN', 'SUBOBSMODE', 'IFNUM', 'PLNUM', 'FDNUM', 'INTNUM', 'CAL', 'PROCSEQN']


#  some more liberal panda dataframe display options
import pandas as pd
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 1000)

# pd.options.display.max_columns = None

#%%  debugging

import dysh
dysh.log.init_logging(3)   # 0=ERROR 1=WARNING 2=INFO 3=DEBUG

#%%  helper functions

def mkdir(name):
    """ simpler frontend for making a directory that might also aready exist
        should also remove files inside
    """
    os.makedirs(name, exist_ok = True)


def vcal(sdf, vane, sky, debug=False):
    """ find the two vane and sky and find tsys for each beam
    """
    kb=['DATE-OBS','SCAN', 'IFNUM', 'PLNUM', 'FDNUM', 'PROCSCAN', 'FEED', 'SRFEED', 'FEEDXOFF', 'FEEDEOFF']
    a = sdf._index[kb]
    feeds = a['FDNUM'].unique()
    feeds.sort()
    for f in feeds:
        b = a.loc[a['FDNUM']==f]
        c1 = b.loc[b['SCAN']==vane]
        c2 = b.loc[b['SCAN']==sky]
        print("FEED",f,c1,c2)
        
        # take aver of all the c1's (vane) and c2's (sky)
        c1_data = 1.0
        c2_data = 2.0
    return
    
    b=a.loc[a['FEEDXOFF']==0.0]
    c=b.loc[b['FEEDEOFF']==0.0]
    d1=c.loc[c['PROCSCAN']=='BEAM1']
    d2=c.loc[c['PROCSCAN']=='BEAM2']
    #
    if len(d1['FDNUM'].unique()) == 1 and len(d2['FDNUM'].unique()) == 1:
        beam1 = d1['FDNUM'].unique()[0]
        beam2 = d2['FDNUM'].unique()[0]
        fdnum1 = d1['FEED'].unique()[0]
        fdnum2 = d2['FEED'].unique()[0]
        if debug:
            print("beams: ",beam1,beam2,fdnum1,fdnum2)
        return [beam1,beam2]
    else:
        # try one other thing
        if len(c['FEED'].unique()) == 2:
            print("getbeam rescued")
            b = c['FEED'].unique() - 1
            return list(b)
        print("too many in beam1:",d1['FDNUM'].unique())
        print("too many in beam2:",d2['FDNUM'].unique())
        return []

def getbeam(sdf, debug=False):
    """ find the two nodding beams
    """
    kb=['DATE-OBS','SCAN', 'IFNUM', 'PLNUM', 'FDNUM', 'PROCSCAN', 'FEED', 'SRFEED', 'FEEDXOFF', 'FEEDEOFF']
    a = sdf._index[kb]
    b=a.loc[a['FEEDXOFF']==0.0]
    c=b.loc[b['FEEDEOFF']==0.0]
    d1=c.loc[c['PROCSCAN']=='BEAM1']
    d2=c.loc[c['PROCSCAN']=='BEAM2']
    #
    if len(d1['FDNUM'].unique()) == 1 and len(d2['FDNUM'].unique()) == 1:
        beam1 = d1['FDNUM'].unique()[0]
        beam2 = d2['FDNUM'].unique()[0]
        fdnum1 = d1['FEED'].unique()[0]
        fdnum2 = d2['FEED'].unique()[0]
        if debug:
            print("beams: ",beam1,beam2,fdnum1,fdnum2)
        return [beam1,beam2]
    else:
        # try one other thing
        if len(c['FEED'].unique()) == 2:
            print("getbeam rescued")
            b = c['FEED'].unique() - 1
            return list(b)
        print("too many in beam1:",d1['FDNUM'].unique())
        print("too many in beam2:",d2['FDNUM'].unique())
        return []

def vanecal(sdf, vane, sky,  feeds=[], tcal=1.0):
    """ loop over feeds to get tsys factor
        for efficiency of large data, it's better to sdf.write() and use only the
        vane and sky scans'
        Example EDGE:  (1.3GB)
             all 163 scans:     6m33s
             write small sdf:     25s
             just vane/cal:       20s
    """
    tsys = np.zeros(len(feeds), dtype=float)
    if False:
        # needs testing, what to do about ifnum, plnum - force one
        mkdir('vanecal')
        sdf.write('vanecal/file.fits',scan=[vane,sky],overwrite=True)
        sdf1 = GBTFITSLoad('vanecal')
    else:
        sdf1 = sdf
    
    i=0
    for f in feeds:
        v = sdf1.gettp(scan=vane, fdnum=f, calibrate=True, cal=False).timeaverage()
        s = sdf1.gettp(scan=sky,  fdnum=f, calibrate=True, cal=False).timeaverage()
        tsys[i] = mean_tsys(v.flux, s.flux, 1.0) * tcal
        # optional?    += tcal/2.0
        i = i + 1
    return tsys
        
#%%  classic tp/ps

f1 = dysh_data(test="getps")        # OnOff, small one (1 IF, 1 PL, 1 INT)
#f1 = dysh_data(example="getps")    # OnOff, long one (4 IFs, 2 PLs, 151 INTs)
#f1 = dysh_data(test="AGBT18B_354_03/AGBT18B_354_03.raw.vegas") # OffOn

print("Using",f1.parts[-1])     # isn't there a better name for this?


sdf1 = GBTFITSLoad(f1)
sdf1.summary()
sdf1.summary(verbose=True)
sdf1._index[k]

sp1 = sdf1.getps().timeaverage()
sp1.plot()
sp1.stats(qac=True)
# 0.21853874407385052 0.6796743343920357 -3.7057502269744873 4.343878746032715

tp1a = sdf1.gettp(scan=152).timeaverage()
tp1b = sdf1.gettp(scan=153).timeaverage()

tp1a.plot()
tp1b.plot()

tp1a.stats(qac=True)  # '490015338.1141468 152969155.2487644 3311173.5 805602304.0'

#%% NOD EXAMPLE-1   tp_nocal    NOD_BEAMS  10,1   (FEED 11,2)

_help = """
   SCAN    OBJECT VELOCITY   PROC  PROCSEQN    RESTFREQ     DOPFREQ # IF # POL # INT # FEED     AZIMUTH   ELEVATIO
0   281      VANE      0.0  Track         1  111.711282  111.711282    1     1     2     16  310.707976  56.549105
1   282       SKY      0.0  Track         1  111.711282  111.711282    1     1     2     16  310.708028  56.549066
2   290  1-631680      0.0    Nod         2  111.711282  111.711282    1     1     6     16    310.3087  55.619947
3   289  1-631680      0.0    Nod         1  111.711282  111.711282    1     1     6     16  310.344245  55.705345
"""


f1 = dysh_data(accept='AGBT22A_325_15/AGBT22A_325_15.raw.vegas')  # accept='nod1'
sdf=GBTFITSLoad(f1)
# 8 files, 16 beams, each file has 2 beams - 4 scans, VANE/SKY/Nod/Nod
sdf.summary()    # 256 rows
# extract 290 and 289 (note order is odd in sdfits:   290 came before 289
mkdir("vane1")
sdf.write('vane1/file.fits',scan=[281,282], overwrite=True)   # 64 rows

vane1 = GBTFITSLoad('vane1')
vane1.summary()
vane1._index[ks]    # 64 rows

vcal(vane1, 281, 282)

v1 = vane1.getspec(0).flux.value
s1 = vane1.getspec(4).flux.value
t1 = s1/(v1-s1)                          # 0.53 +/0 0.02
np.nanmean(t1[100:900])
np.nanstd(t1[100:900])

v1 = vane1.gettp(scan=281, fdnum=8, calibrate=True, cal=False).timeaverage()
# ok
s1 = vane1.gettp(scan=282, fdnum=8, calibrate=True, cal=False).timeaverage()
# ok

sp = sdf.gettp(scan=281, fdnum=8, calibrate=True, cal=False).timeaverage().plot()
# ok

print(vanecal(vane1, 281, 282, feeds=range(16)))
# [1.21549487 1.07764171 4.13970304 1.15904234 0.9971537  1.122969
# 0.99023236 1.01537602 1.02035436 0.98760956 1.00412189 1.05833454
# 0.99750374 1.02722791 1.00445433 1.00216715]


#%%   EDGE data  (1363MB)

_help = """
     SCAN    OBJECT VELOCITY       PROC  PROCSEQN    RESTFREQ     DOPFREQ # IF # POL # INT # FEED     AZIMUTH   ELEVATIO
0      17      VANE      0.0      Track         1  113.571858  113.571858    1     1    21     16  165.049197  79.052319
1      18       SKY      0.0      Track         1  113.571858  113.571858    1     1    21     16  165.050181  79.052262
2      19   NGC0001      0.0        Nod         1  113.571858  113.571858    1     1   122     16  166.640475   79.14949
3      20   NGC0001      0.0        Nod         2  113.571858  113.571858    1     1   122     16  168.125608  79.199539
4      21      VANE      0.0      Track         1  113.571858  113.571858    1     1    21     16  178.544679  79.395538
5      22       SKY      0.0      Track         1  113.571858  113.571858    1     1    21     16  178.545661  79.395501
6      23   NGC0001      0.0  DecLatMap         1  113.571858  113.571858    1     1    68     16  177.410654  78.233486
"""

sdf2 = GBTFITSLoad(dysh_data('AGBT21B_024_01/AGBT21B_024_01.raw.vegas'))
a = sdf2.summary()  # 208 scans   1363MB
b=a[a["OBJECT"] == "VANE"]
print(b)   # there are 13 vane's in here

# make a smaller dataset for testing
mkdir("edge1")
sdf2.write("edge1/file.fits", fdnum=8, scan=[17,18,19,20], intnum=[0,1], overwrite=True)
sdf2.write("edge1/file.fits", fdnum=8, scan=[17,18,19,20], overwrite=True)
sdf2.write("edge1/file.fits", scan=[17,18,19,20], overwrite=True)                # 4576 rows
#  CPU times: user 23.5 s, sys: 3.02 s, total: 26.5 s   Wall time: 25.2 s

edge1 = GBTFITSLoad("edge1")
edge1.summary()
edge1._index[ks]

tsys = vanecal(edge1, 17, 18, feeds=range(16))
#  sdf2   CPU times: user 6min 29s, sys: 10.5 s, total: 6min 40s   Wall time: 6min 33s
#  edge1  CPU times: user   20.1 s, sys: 271 ms, total:   20.4 s   Wall time:   20.1 s
print(tsys)
# array([1.26915776, 1.14312562, 1.13592276, 1.1619848 , 1.11120068,
#       1.20323269, 1.12788577, 1.13010956, 1.10427554, 1.1221987 ,
#       1.15696975, 1.15716321, 1.0804609 , 1.16266635, 1.0460511 ,
#       1.09217787])
#
 
# plotting a passband
edge1.gettp(scan=17, fdnum=8, calibrate=True, cal=False).timeaverage().plot()
# ok !!

# not working yet, since it wants ON/OFF, needs the GETTP hack
sp2 = sdf2.getnod(scan=19)
# IndexError: index 0 is out of bounds for axis 0 with size 0

#%% NOD EXAMPLE-3 tp_nocal   NOD_BEAMS 0,1


f3 = dysh_data(accept='AGBT15B_244_07/AGBT15B_244_07.raw.vegas')
sdf=GBTFITSLoad(f3)
# 8 fits files,   2 for beams, 4 for IF  - 12 scans (2 CALSEQ)
sdf.summary()
# 11072 rows

mkdir("nod3cal")
sdf.write("nod3cal/file.fits", scan=130, ifnum=0, plnum=0, overwrite=True)
nod3cal = GBTFITSLoad("nod3cal")
nod3cal.summary()
nod3cal._index[ks]   # 82 rows

v1 = nod3cal.gettp(scan=130, fdnum=0, calibrate=True, cal=False).timeaverage()
s1 = nod3cal.gettp(scan=130, fdnum=1, calibrate=True, cal=False).timeaverage()
v1.plot()
s1.plot()
t1 = s1/(v1-s1) 

tsys = vanecal(nod3cal, 130, 130, feeds=[0,1])

mkdir("nod3")
sdf.write('nod3/file.fits', scan=[131,132], ifnum=0, plnum=0, overwrite=True)  #244

 
nod3 = GBTFITSLoad('nod3')
nod3.summary()
nod3._index[k]    # 244 rows
getbeam(nod3)     # [0,1]



nod3.gettp(scan=131,fdnum=0,calibrate=True, cal=False).timeaverage().plot()
nod3.gettp(scan=132,fdnum=0,calibrate=True, cal=False).timeaverage().plot()
nod3.gettp(scan=131,fdnum=1,calibrate=True, cal=False).timeaverage().plot()
# ok

#%% issue 257

_help = """
   SCAN    OBJECT VELOCITY   PROC  PROCSEQN RESTFREQ DOPFREQ # IF # POL # INT # FEED AZIMUTH ELEVATIO
0     1  B0329+54      0.0  Track         1      1.5     1.5    1     2    16      1     0.0      0.0
1     2  B0329+54      0.0  Track         1      1.5     1.5    1     2    16      1     0.0      0.0
2     3  B0329+54      0.0  Track         1      1.5     1.5    1     2    16      1     0.0      0.0
3     4  B0329+54      0.0  Track         1      1.5     1.5    1     2    16      1     0.0      0.0
4     5  B0329+54      0.0  Track         1      1.5     1.5    1     2    16      1     0.0      0.0
"""

f4 = dysh_data("TSCI_RYAN_16/TSCI_RYAN_16.raw.vegas")
sdf4 = GBTFITSLoad(f4)
sdf4.summary()
sdf4._index[ks]   # 192 rows

sdf4.getps()
# Exception: Multiple SUBOBSMODE present, cannot deal with this yet ['TPWCAL', 'TPNOCAL']

# scan 1 is the TPWCAL
sp4a = sdf4.gettp(scan=1, plnum=0)
sp4b = sdf4.gettp(scan=1, plnum=1)               
    
sp4a.timeaverage().plot(ymin=-1e8, ymax=1e9, title="TSCI_RYAN_16 scan=1 plnum=0", xaxis_unit="chan")
sp4b.timeaverage().plot(ymin=-1e8, ymax=1e9, title="TSCI_RYAN_16 scan=1 plnum=1", xaxis_unit="chan")
     
sp = sdf4.gettp(scan=2, plnum=0, calibrate=True, cal=False).timeaverage().plot(ymin=-1e8, ymax=1e9)
# ok

