#!/usr/bin/env python3

"""
Description:
------------
This is a script developed in spyder during the nodding work.

"""

import numpy as np
import numpy.ma as ma
from scipy.stats import norm
import astropy.units as u
from astropy.io import fits

import matplotlib.pyplot as plt
from dysh.fits.sdfitsload import SDFITSLoad
from dysh.fits.gbtfitsload import GBTFITSLoad
import dysh.util as util
from dysh.util.files import dysh_data
from dysh.util.selection import Selection

#  useful keys for a mult-beam observation

k=['DATE-OBS','SCAN', 'IFNUM', 'PLNUM', 'FDNUM', 'INTNUM', 'PROCSCAN','FEED', 'SRFEED', 'FEEDXOFF', 'FEEDEOFF', 'SIG', 'CAL', 'PROCSEQN', 'PROCSIZE']
ks=['DATE-OBS','SCAN', 'IFNUM', 'PLNUM', 'FDNUM', 'INTNUM', 'CAL', 'PROCSEQN']


#  some more liberal panda dataframe display options
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 1000)

# pd.options.display.max_columns = None

#%%  getnod (fake)

def getnod(sdf, scan=None, fdnum=[0,1]):
    """  a fake getnod type function
            scan=  needs list of 2 (for now just a nodding pair)
            fdnum= needs list of 2 (often this is 0,1)
    """
    p1a = sdf.gettp(scan=scan[0],fdnum=fdnum[0])[0].timeaverage()
    p1b = sdf.gettp(scan=scan[1],fdnum=fdnum[0])[0].timeaverage() 
    t1 = (p1a-p1b)/p1b
    p2a = sdf.gettp(scan=scan[1],fdnum=fdnum[1])[0].timeaverage()
    p2b = sdf.gettp(scan=scan[0],fdnum=fdnum[1])[0].timeaverage()
    t2 = (p2a-p2b)/p2b

    t = (t1+t2)/2
    
    return t

def getbeam(sdf):
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
        print("beams: ",beam1,beam2,fdnum1,fdnum2)
        return [beam1,beam2]
    else:
        print("too many in beam1:",d1['FDNUM'].unique())
        print("too many in beam2:",d2['FDNUM'].unique())
        return []

 
#%%  classic tp/ps
   
f1 = dysh_data(test="getps")        # OnOff, small one
#f1 = dysh_data(example="getps1")   # OnOff, long one (4 IFs, ...)
#f1 = dysh_data(test="AGBT18B_354_03/AGBT18B_354_03.raw.vegas") # OffOn

print("Using",f1.parts[-1])     # isn't there a better name for this?


sdf1 = GBTFITSLoad(f1)
sdf1.summary()
sdf1.summary(verbose=True)
sdf1._index[k]

#p1 = sdf1.gettp

#%%

nbox = 127

p1 = sdf1.gettp()
sp0 = p1[0].timeaverage()  # tsys=17.45805259
sp1 = p1[1].timeaverage()  #      17.24000331
p1a = (sp0-sp1)/sp0
sp3 = p1a.smooth('box',nbox)
sp3.plot()

#%%

p2 = sdf1.getps()
sp2 = p2[0].timeaverage()
sp4 = sp2.smooth('box',nbox)
sp4.plot()


#%%   dividing
sp5 = sp4/sp3
sp5.plot()
sp5.stats()    # about 19.7K

#%%  full example of NGC2415

# timing:  loading ~1sec  getps chain: ~8.4s (getps is most of it)
#          getps() on all 5 IFs and 2 PLs takes 73sec, very linear in #IF and #PL
sdf1 = GBTFITSLoad(dysh_data(example="getps1"))
sdf1.getps(plnum=0,ifnum=0).timeaverage().smooth('box',51).plot(xaxis_unit='km/s')


#%% EXAMPLE-0  fs/nod  NOD_BEAMS 2,6

#  example from
#  https://gbtdocs.readthedocs.io/en/latest/how-tos/data_reduction/gbtidl.html
#  TGBT22A_503_02.raw.vegas
#  in here one frequency-switched scan (#64), 
#  and two nodding scans (#62 and #63)
#  7 SDFITS files, one for each beam, are stored in the data directory

# this is a huge file, we preload the minimum number of scans

f0 = dysh_data(example="nod-KFPA/data/TGBT22A_503_02.raw.vegas")     # example='nodfs'
sdf0 = GBTFITSLoad(f0)
sdf0.summary()
getbeam(sdf0)

# Loaded 7 FITS files
# CPU times: user 17.9 s, sys: 5.54 s, total: 23.4 s
# Wall time: 1min 22s
# 18GB VIRT, 5G RES 4G SHR

#     SCAN OBJECT VELOCITY   PROC  PROCSEQN   RESTFREQ    DOPFREQ # IF # POL # INT # FEED     AZIMUTH   ELEVATIO
# 0     60   W3_1    -40.0  OnOff         1  23.959156  23.694495    6     2    31      7  324.227878  38.705977
# 1     61   W3_1    -40.0  OnOff         2  23.959156  23.694495    6     2    31      7  324.152607  39.012111
# 2     62   W3_1    -40.0    Nod         1  23.959156  23.694495    6     2    31      7  324.274257  38.419406
# 3     63   W3_1    -40.0    Nod         2  23.959156  23.694495    6     2    31      7  324.367171  38.285767
# 4     64   W3_1    -40.0  Track         1  23.959156  23.694496    6     2    21      7  324.388427  38.083529
# 5     65   W3_1    -40.0  OnOff         1  23.787811  23.694495    5     2    31      7    324.6698   37.10657
# 6     66   W3_1    -40.0  OnOff         2  23.787811  23.694495    5     2    31      7  324.574402  37.418884
# 7     67   W3_1    -40.0    Nod         1  23.787811  23.694495    5     2    31      7  324.732475  36.827076
# 8     68   W3_1    -40.0    Nod         2  23.787811  23.694495    5     2    31      7  324.831553  36.694947
# 9     69   W3_1    -40.0  Track         1  23.787811  23.694496    5     2    31      7  324.874404  36.463014

!mkdir -p nod0


#  can this use less memory? still seems to sit at 35GB virtual
sdf0.write('nod0/file.fits',scan=[62,63],ifnum=0, plnum=0, overwrite=True)
sdf0.write('nod0/file.fits',scan=[67,68],ifnum=0, plnum=0, overwrite=True)
# 5IF: 191404 -rw-rw-r-- 1 teuben teuben 195989760 Sep 13 14:45 junk620.fits
# 1IF         -rw-rw-r-- 1 teuben teuben  32685120 Sep 13 14:59 junk620.fits

# IF=6 POL=2 INT=31 FEED=2 SCAN=2


#%%

nod0 = GBTFITSLoad('nod0')
nod0.summary()
nod0._index[k]
getbeam(nod0)

sp = getnod(nod0,scan=[62,63])
sp.smooth('box',151).plot()

#%%
#     now with the multifile trick
sdf62 = GBTFITSLoad('junk62')
sdf62.summary()
sdf62._index[k]    # 2975 for 2 beams   1736 for 7 beams   1*2*31*2*2

p1a = sdf62.gettp(scan=62, fdnum=2, plnum=0)[0].timeaverage().smooth('box',5)
p1b = sdf62.gettp(scan=63, fdnum=2, plnum=0)[0].timeaverage().smooth('box',5)
t1 = (p1a-p1b)/p1b
p2a = sdf62.gettp(scan=63, fdnum=6, plnum=0)[0].timeaverage().smooth('box',5)
p2b = sdf62.gettp(scan=62, fdnum=6, plnum=0)[0].timeaverage().smooth('box',5)
t2 = (p2a-p2b)/p2b

t = (t1+t2)/2
ts = t.smooth('box',10)
ts.plot(xaxis_unit="km/s", xmin=440, xmax=530)
    
#%% EXAMPLE-1   tp_nocal    NOD_BEAMS  10,1   (FEED 11,2)

f1 = dysh_data(accept='AGBT22A_325_15/AGBT22A_325_15.raw.vegas')  # accept='nod1'
sdf=GBTFITSLoad(f1)
# 8 files, 16 beams, each file has 2 beams - 4 scans, VANE/SKY/Nod/Nod
sdf.summary()
# extract 290 and 289 (note order is odd in sdfits:   290 came before 289
!mkdir -p nod1
sdf.write('nod1/file.fits',scan=[289,290], overwrite=True)   # 192 rows

nod1 = GBTFITSLoad('nod1')
nod1.summary()
nod1._index[k]    # 24 rows
getbeam(nod1)

sp = getnod(nod1,scan=[289,290])
# Exception: unbalanced calrows 0 != 6
sp.smooth('box',20).plot(xaxis_unit='km/s')

#%% EXAMPLE-2 tp (deprecated)

f2 = dysh_data(accept='TREG_050627/TREG_050627.raw.acs/TREG_050627.raw.acs.fits')
sdf=GBTFITSLoad(f2)
sdf.summary()
# extracts 182 and 183 for testing
sdf.write('junk182.fits', scan=[182,183])

sdf182 = GBTFITSLoad('junk182.fits')
sdf182.summary()
sdf182._index[k]   # 128 rows

# data[beam=2][scan=2][if=2][time=4][pol=2][cal=2]
# deprecarted, don't use

#%% EXAMPLE-3 tp_nocal   NOD_BEAMS 0,1


f3 = dysh_data(accept='AGBT15B_244_07/AGBT15B_244_07.raw.vegas')
sdf=GBTFITSLoad(f3)
# 8 fits files,   2 for beams, 4 for IF  - 12 scans (2 CALSEQ)
sdf.summary()
# 11072 rows
!mkdir -p nod3
sdf.write('nod3/file.fits', scan=[131,132], ifnum=0, plnum=0, overwrite=True)  #244

 
nod3 = GBTFITSLoad('nod3')
nod3.summary()
nod3._index[k]    # 244 rows
getbeam(nod3)

# where did the IF go?   1310 is for one beam (there are 7)
# 
sp = getnod(nod3,scan=[131,132])
# Exception: unbalanced calrows 0 != 6
sp.smooth('box',20).plot(xaxis_unit='km/s')


#%% EXAMPLE-4  tp   NOD_BEAMS 0,1

f4 = dysh_data(accept='TGBT18A_500_06/TGBT18A_500_06.raw.vegas')   # Arp220
sdf = GBTFITSLoad(f4)
# 2 files - 10 scans written in odd order (64,65,66 first)
sdf.summary()

!mkdir -p nod4
sdf.write('nod4/file.fits', scan=[57,58],plnum=1,overwrite=True)   # 240

nod4 = GBTFITSLoad('nod4')
nod4.summary()
nod4._index[k]   # 240 rows

sp = getnod(nod4,scan=[57,58])
sp.smooth('box',20).plot(xaxis_unit='km/s')

# data[int=30,30][pol=2][cal=2] 

#%% EXAMPLE-5  tp (deprecated)

f5 = dysh_data(accept='TSCAL_19Nov2015/TSCAL_19Nov2015.raw.acs/TSCAL_19Nov2015.raw.acs.fits')
sdf=GBTFITSLoad(f5)
# 1 file, 1 scan
sdf.summary()
sdf._index[k]    # 64 rows

# data[scan=1][beam=2][if=2][time=4][pol=2][cal=2]

# NOTE: this data incomplete? only has procseqn=1

# deprecated, don't use it

#%% EXAMPLE-6  tp   NOD_BEAMS 0,1  but both PL's needed

f6 = dysh_data(accept='AGBT17B_319_06/AGBT17B_319_06.raw.vegas')
sdf = GBTFITSLoad(f6)
# 3 files,  one for each IF
sdf.summary()
getbeam(sdf)


!mkdir -p nod6
sdf.write('nod6/file.fits', scan=[9,10], ifnum=0, overwrite=True)   # 248
# scans on AGPeg number 9..54 in nod-pairs

nod6 = GBTFITSLoad('nod6')
nod6.summary()
nod6._index[k]   # 248 rows     2*31*1*2*2

sp = getnod(nod6, scan=[9,10])
sp.smooth('box',151).plot(xaxis_unit='km/s')
 
#  data   [if=3][scan=2][int=31][pol=2][cal=2]
#  scan=9 fdnum=0      scan=10 fdnum=1




#%% EXAMPLE-7  tp   NOD_BEAMS 0,1   (but called FEED=3,7)

f7 = dysh_data(accept='TGBT21A_501_10/TGBT21A_501_10.raw.vegas')   # NGC_2639 
sdf=GBTFITSLoad(f7)
# 2 files
sdf.summary()
getbeam(sdf)

!mkdir -p nod7
sdf.write('nod7/file.fits', scan=[36,37], plnum=1, overwrite=True)  # 128 selected

nod7 = GBTFITSLoad('nod7')
nod7.summary()
nod7._index[k]  

sp = getnod(nod7, scan=[36,37])
sp.smooth('box',151).plot(xaxis_unit='km/s')

#   data   [ scan=2] [int=16] [pol=2] [cal=2]

#%% EXAMPLE-8            NOD_BEAMS 0,1   (FEED 1,4)
 
f8 = dysh_data(accept="AGBT19A_340_07/AGBT19A_340_07.raw.vegas")
sdf = GBTFITSLoad(f8)
# 8 files;   4IF 2PL 6-INT 2FEED   (4IF * 2FD in the 8 files)
sdf.summary()
sdf._index[k]     # 768
getbeam(sdf)


!mkdir -p nod8
sdf.write('nod8/file.fits', scan=[43,44], ifnum=0, plnum=0, overwrite=True)  # 48


nod8= GBTFITSLoad('nod8')
nod8.summary()
nod8._index[k]   # 48 rows    (2*6*2*2)

sp = getnod(nod8, scan=[43,44])
sp.smooth('box',151).plot(xaxis_unit='km/s')


#%% EXAMPLE-9     NOD_BEAMS 0,1   (but unknown PROCSCAN)

f9 = dysh_data(accept="AGBT12A_076_05/AGBT12A_076_05.raw.acs")
sdf = GBTFITSLoad(f9)
# 1 file
sdf.summary()    # 3072 rows       4*2*12*2*2*8
sdf._index[k]   # PROCSCAN="unknown"
getbeam(sdf)     # 
# 4-IF 2-POL 12-INT 2-FEED

!mkdir -p nod9
sdf.write('nod9/file.fits', scan=[12,13], ifnum=3, plnum=1, overwrite=True)   # 96

nod9= GBTFITSLoad('nod9')
nod9.summary()
nod9._index[k]     # 96/    768 rows    [scan=2][if=4] [int=12] [pol=2] [cal=2]
# 2*12*2*2
# very odd order of scans  -  issues/376

sp = getnod(nod9,scan=[12,13])
sp.smooth('box',20).plot(xaxis_unit='km/s')


