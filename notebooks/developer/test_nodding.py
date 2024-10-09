#!/usr/bin/env python3

"""
This script was developed in spyder during the nodding work, 
it also contains a local getnod() , which works but it biased since it used gettp().
10 dataset were tested, see final column how well it works

nod0  fdnum=[2,6]       16GB memory issue?               tp cal; OK
nod1  fdnum=[10,1]      unbalanced calrows 0 != 6        tp nocal;  TBD
nod2                    deprecated                       -
nod3  fdnum=[0,1]       unbalanced calrows 0 != 6        tp noal;  TBD
nod4  fdnum=[0,1]       odd scan order                   tp cal; getnod OK, but weird Arp220 signal
nod5                    deprecated                       -
nod6  fdnum=[0,1]       both PL's needed!                tp cal; getnod failing multiple feeds/scan
nod7  fdnum=[0,1]       last scan short, weak signal     tp cal; OK (weak line?)
nod8  fdnum=[0,1]       no signal?                       tp cal; OK (no signal)
nod9  fdnum=[0,1]       but PROCSCAN unknown,            tp cal; getnod failing multiple feeds/scan   #376,#390

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

#  useful keys for a mult-beam observation listing

k=['DATE-OBS','SCAN', 'IFNUM', 'PLNUM', 'FDNUM', 'INTNUM', 'PROCSCAN','FEED', 'SRFEED', 'FEEDXOFF', 'FEEDEOFF', 'SIG', 'CAL', 'PROCSEQN', 'PROCSIZE']
ks=['DATE-OBS','SCAN', 'IFNUM', 'PLNUM', 'FDNUM', 'INTNUM', 'CAL', 'PROCSEQN']


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
    """
    os.makedirs(name, exist_ok = True)

def getnod(sdf, scan=None, fdnum=[0,1], ps=True):
    """  a fake getnod type function using a series of gettp()
            scan=  needs list of 2 (for now just a nodding pair)
            fdnum= needs list of 2 (often this is 0,1)
         fake because we do the timeaverage in power mode, not spectral mode
         same for averaging their tsys.  short of using getspec() of course.
         By setting ps=False, one can use the "Average BS" mode, instead
         of the default "Average PS" mode.
    """
    s1 = sdf.gettp(scan=scan[0],fdnum=fdnum[0])[0]
    t1s = s1.tsys.mean()   # should use np.nanmean()
    s1 = s1.timeaverage()
    
    r1 = sdf.gettp(scan=scan[1],fdnum=fdnum[0])[0]
    t1r = r1.tsys.mean()
    r1 = r1.timeaverage()
    
    s2 = sdf.gettp(scan=scan[1],fdnum=fdnum[1])[0]
    t2s = s2.tsys.mean()
    s2 = s2.timeaverage()
    
    r2 = sdf.gettp(scan=scan[0],fdnum=fdnum[1])[0]
    t2r = r2.tsys.mean()
    r2 = r2.timeaverage()

    if ps:
        tsys1 = 0.5*(t1s+t1r)
        tsys2 = 0.5*(t2s+t2r)
        #
        t1 = tsys1*(s1-r1)/r1
        t2 = tsys2*(s2-r2)/r2
    else:        
        tsys1 = 0.5*(t1s+t2r)
        tsys2 = 0.5*(t2s+t1r)
        #
        t1 = tsys1*(s1-r2)/r2
        t2 = tsys2*(s2-r1)/r1

    t = (t1+t2)/2
    
    print("fake getnod: tsys=%.2f %.2f (PS mode=%s)" % (tsys1,tsys2,repr(ps)))
    
    return t

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

 
print("Added getnod() and getbeam() for local testing")
#%%  classic tp/ps

f1 = dysh_data(test="getps")        # OnOff, small one (1 IF, 1 PL, 1 INT)
#f1 = dysh_data(example="getps")    # OnOff, long one (4 IFs, 2 PLs, 151 INTs)
#f1 = dysh_data(test="AGBT18B_354_03/AGBT18B_354_03.raw.vegas") # OffOn

print("Using",f1.parts[-1])     # isn't there a better name for this?


sdf1 = GBTFITSLoad(f1)
sdf1.summary()
sdf1.summary(verbose=True)
sdf1._index[k]


#%%

nbox = 127
inc=[[1500,15000],[18000,28000]]

p1 = sdf1.gettp(ifnum=0, plnum=0)
sp0 = p1[0].timeaverage()  # tsys=17.45805259
sp1 = p1[1].timeaverage()  #      17.24000331
tsys1 = p1[0].tsys.mean()
tsys2 = p1[1].tsys.mean()
print('TSYS:',tsys1,tsys2)
tsys = (tsys1+tsys2)/2
p1a = tsys*(sp0-sp1)/sp0
p1a.baseline(1,include=inc,remove=True)
sp3 = p1a.smooth('box',nbox)
sp3._unit = u.K   # hack
sp3.plot(xaxis_unit="km/s", xmin=3800-1000,xmax=3800+1000, ymin=-0.05, ymax=0.15)

#%%

p2 = sdf1.getps(ifnum=0, plnum=0)
sp2 = p2[0].timeaverage()
sp2.baseline(1,include=inc,remove=True)
sp4 = sp2.smooth('box',nbox)
sp4.plot(xaxis_unit="km/s", xmin=3800-1000,xmax=3800+1000, ymin=-0.05, ymax=0.15)


#%%   differencing
sp5 = sp4-sp3
sp5.plot(xaxis_unit="km/s", xmin=3800-1000,xmax=3800+1000)
#  bias in this difference is about  0.002, on a peak of 0.10 K, or 2%

#%%  full example of NGC2415, gives nice spectrum

# timing:  loading ~1sec  getps chain: ~8.4s (getps is most of it)
#          getps() on all 5 IFs and 2 PLs takes 73sec, very linear in #IF and #PL
sdf1 = GBTFITSLoad(dysh_data(example="getps"))
sdf1.getps(plnum=0,ifnum=0).timeaverage().smooth('box',51).plot(xaxis_unit='km/s')

#%% 

GBTFITSLoad('nod0_test60.fits').getspec(0).plot(title='GBTIDL getps',  xaxis_unit="km/s")
GBTFITSLoad('nod0_test62.fits').getspec(0).plot(title='GBTIDL getnod', xaxis_unit="km/s")
GBTFITSLoad('nod0_test64.fits').getspec(0).plot(title='GBTIDL getfs')











#%% EXAMPLE-0  ps/fs/nod  NOD_BEAMS 2,6

#  example taken from TGBT22A_503_02.raw.vegas
#  https://gbtdocs.readthedocs.io/en/latest/how-tos/data_reduction/gbtidl.html
#  TGBT22A_503_02.raw.vegas
#  in here one frequency-switched scan (#64), 
#  and two nodding scans (#62 and #63) as well get position-switches (#60,61)
#  7 SDFITS files, one for each KFPA beam, are stored in the data directory

gbtidl0 = """
filein, "nod-KFPA/data/TGBT22A_503_02.raw.vegas"

sclear
getsigref, 62, 63, fdnum=2, ifnum=0, plnum=0
accum
getsigref, 63, 62, fdnum=6, ifnum=0, plnum=0
accum
ave
fileout,"nod0_test62.fits"
keep

sclear
getfs,64,fdnum=0
fileout,"nod0_test64.fits"
keep

sclear
getps,60,plnum=0,ifnum=0
fileout,"nod0_test60.fits"
keep
"""
#   BUG?   there are no units in these keep files.....

print(gbtidl0)

# this is a huge 19GB file, we preload the minimum number of scans - it also needs 32GB memory to be speedy
f0 = "nod-KFPA/data/TGBT22A_503_02.raw.vegas.trim.fits"  #  4 scans (60..63)
f0 = "nod-KFPA/data/TGBT22A_503_02.raw.vegas"            # 25 scans
f0 = dysh_data(example=f0)                        # AKA example='nod'
sdf0 = GBTFITSLoad(f0)                            # 7 FITS files
sdf0.summary()                                    # 25 scans 60..84 
getbeam(sdf0)                                     # [2,6]

# Sample times:
# CPU times: user 17.9 s, sys: 5.54 s, total: 23.4 s   Wall time: 1min 22s  on 16MB
# CPU times: user 19.3 s, sys: 1.31 s, total: 20.7 s   Wall time: 19.9 s    on 32MB
# CPU times: user 19.4 s, sys: 5.06 s, total: 24.5 s   Wall time: 1min 32s  on 32MB first load
# CPU times: user 18.9 s, sys: 1.3 s, total: 20.2 s    Wall time: 19.3 s    on 32MB 2nd load

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

#%% comparing GETNOD without any smoothing

sp0 = sdf0.getnod(scan=[62,63],ifnum=0,plnum=0).timeaverage()
sp1 = GBTFITSLoad('nod0_test62.fits').getspec(0)   # results from two getsigref's in GBTIDL
d = sp0.flux.value - sp1.flux.value
print("getnod rms diff",np.nanstd(d) )     #  0.0001160868980983191
np.where(np.isnan(d))                      #  [   0, 9216])
sp0.stats(qac=True)    # 0.3748503311988827  0.35504111651677267 -1.2325043436127545 2.1533314871468248
sp1.stats(qac=True)    # 0.37477970123291016 0.35496729612350464 -1.2321275472640991 2.1529178619384766
#%% comparing GETPS

sp0 = sdf0.getps(scan=60,plnum=0,ifnum=0,fdnum=0).timeaverage()
sp1 = GBTFITSLoad('nod0_test60.fits').getspec(0)
d = sp0.flux.value - sp1.flux.value
print("getnod rms diff", np.nanstd(d))    #  1.716955005887903e-05     trim:  6.906450791469074e-08
np.where(np.isnan(d))

# 1.716728447608622e-05


# below older comparisons with smoothing in the middle.

sdf0.write('nod0-scan60.fits',scan=[60,61],ifnum=0,plnum=0,fdnum=0,overwrite=True,multifile=False)

sdf0.write('nod0-scan60b.fits',scan=[62,63],ifnum=0,plnum=0,fdnum=[2,6],overwrite=True,multifile=False)
# IndexError: list index out of range
sdf0.write('nod0-scan60b.fits',scan=[62,63],ifnum=0,plnum=0,fdnum=[2,6],overwrite=True)

#%%  comparing DYSH and GBTIDL (ps and nod)

sdf0.getnod(scan=[62,63],ifnum=0,plnum=0).timeaverage().smooth('box',11).plot(title='dysh getnod',xaxis_unit="km/s")
sdf0.getnod(scan=62,ifnum=0,plnum=1).timeaverage().smooth('box',101).plot(xaxis_unit="km/s")
sdf0.getnod(scan=67,ifnum=0,plnum=0).timeaverage().smooth('box',101).plot(xaxis_unit="km/s")
sdf0.getnod(scan=72,ifnum=0,plnum=0).timeaverage().smooth('box',101).plot(xaxis_unit="km/s")
sdf0.getnod(scan=77,ifnum=0,plnum=0).timeaverage().smooth('box',101).plot(xaxis_unit="km/s")
sdf0.getnod(scan=82,ifnum=0,plnum=0).timeaverage().smooth('box',101).plot(xaxis_unit="km/s")

# [75,76] and [80,81] are very weak

# there is no support for pol a
sp0 = sdf0.getps(scan=60,plnum=0,ifnum=0,fdnum=0).timeaverage().smooth('box',11)
sp1 = sdf0.getps(scan=60,plnum=1,ifnum=0,fdnum=0).timeaverage().smooth('box',11)
w0 = sp0.meta["TSYS"]**-2
w1 = sp1.meta["TSYS"]**-2
print("Tsys for two pols:",sp0.meta["TSYS"],sp1.meta["TSYS"])    # 56.681667064114 57.261110023956014
# pick
sp_p1 = (sp0+sp1)/2                  # straight averaging of polarizations
# or
sp_p1 = (sp0*w0 + sp1*w0)/(w0+w1)    # properly weighted 


sp_p1.plot(title='dysh getps',xaxis_unit="km/s")

# both cases return units "ct", which means we cannot subtract the spectra
sp_p2 = GBTFITSLoad('nod0_test60.fits').getspec(0)
sp_p2 = SDFITSLoad('nod0_test60.fits').getspec(0)

d_p = sp_p1.flux.value - sp_p2.flux.value
print("getps rms diff",d_p[1:-1].std() )     #  0.0014509875806096991

sp_n1 = sdf0.getnod(scan=[62,63],ifnum=0,plnum=0).timeaverage().smooth('box',11)
sp_n2 = GBTFITSLoad('nod0_test62.fits').getspec(0)

d_n = sp_n1.flux.value - sp_n2.flux.value
print("getnod rms diff", d_n[1:-1].std())    #  0.00016841065800563842

#%% old and new gbtidl
old="nod0_low_accuracy/"
new="./"
sp1 = GBTFITSLoad(new+'nod0_test62.fits').getspec(0)
sp2 = GBTFITSLoad(old+'nod0_test62.fits').getspec(0)
(sp1-sp2).stats()    # rms 2.0514587e-08

#%%

#  can this use less memory? still seems to sit at 35GB virtual
mkdir("nod0")
sdf0.write('nod0/file.fits',scan=[62,63],ifnum=0, plnum=0, overwrite=True)
# sdf0.write('nod0/file.fits',scan=[67,68],ifnum=0, plnum=0, overwrite=True)
# 5IF: 191404 -rw-rw-r-- 1 teuben teuben 195989760 Sep 13 14:45 junk620.fits
# 1IF         -rw-rw-r-- 1 teuben teuben  32685120 Sep 13 14:59 junk620.fits

# IF=6 POL=2 INT=31 FEED=2 SCAN=2

mkdir("nod0fs")
sdf0.write('nod0fs/file.fits',scan=[64],ifnum=0, plnum=0, fdnum=0, overwrite=True)


#%% nod0: ps mode

nod0 = GBTFITSLoad('nod0')
nod0.summary()
nod0._index[k]   # 868 rows (124 * 7)
getbeam(nod0)

sp = getnod(nod0,scan=[62,63], fdnum=[2,6])
# tsys=63.63 73.10
sp.smooth('box',51).plot(xaxis_unit="km/s", xmin=-100, xmax=50, title='getnod fake bs=False')

sp = nod0.getnod(scan=62).timeaverage()
sp.stats(qac=True)  # 0.3748503311988827 0.35504111651677267 -1.2325043436127545 2.1533314871468248
sp.stats(qac=True, test='0.3748503311988827 0.35504111651677267 -1.2325043436127545 2.1533314871468248')

sp.smooth('box',51).plot(xaxis_unit="km/s", xmin=-100, xmax=50, title='sdf::getnod')
#%% nod0: bs mode

sp = getnod(nod0,scan=[62,63], fdnum=[2,6], ps=False)
# tsys=68.78 67.94
sp.smooth('box',51).plot(xaxis_unit="km/s", xmin=-100, xmax=50, title='getnod fake bs=True')

#%%
nod0fs = GBTFITSLoad('nod0fs')
nod0fs.summary()
nod0fs._index[k]   # 84 rows   - now 168 ?

p1  = nod0fs.getfs(plnum=0)[0]
s1=p1.timeaverage()
s1.smooth('gauss',5).plot(xaxis_unit="km/s", xmin=-100, xmax=50, title="getfs-1")

p2  = nod0fs.getfs(plnum=1)[0]
s2=p2.timeaverage()
s2.smooth('gauss',5).plot(xaxis_unit="km/s", xmin=-100, xmax=50, title="getfs-2")

s3 = (s1+s2)/2    # @todo note we're doing unweighted average
s3.smooth('gauss',51).plot(xaxis_unit="km/s", xmin=-100, xmax=50, title="getfs-3")
#   @todo   this looks like it's shifted from the figure in gbtidl manual !!!   3.3 MHz ???
s3.smooth('gauss',5).plot(xaxis_unit="GHz", xmin=23.685, xmax=23.705) 
#   why do the negative components looks so strong



#%% issue 342
filename = '/home/sdfits/AGBT24A_329_02/AGBT24A_329_02.raw.vegas'
sdfits = GBTFITSLoad(filename)
sdfits.write("AGBT24A_329_02.raw.sub.vegas.fits", plnum=1, intnum=0, overwrite=False)
sdf = GBTFITSLoad("AGBT24A_329_02.raw.sub.vegas4.fits")
sdf.summary()




#%%  PS mode compared in several ways [issue #342]

mkdir("nod0ps")
sdf0.write('nod0ps/file.fits',scan=[60,61], plnum=0, ifnum=0, fdnum=0, overwrite=True)

mkdir("nod0ps_multi")
sdf0.write('nod0ps_multi/file.fits',scan=[60,61], plnum=0, ifnum=0, fdnum=[0,1], overwrite=True)
GBTFITSLoad('nod0ps_multi')
# ok

mkdir("nod0ps_single")
sdf0.write('nod0ps_multi/file.fits',scan=[60,61], plnum=0, ifnum=0, fdnum=[0,1], overwrite=True, multifile=False)
GBTFITSLoad('nod0ps_single')
# AttributeError: 'Selection' object has no attribute 'TIMESTAMP'


nod0ps = GBTFITSLoad('nod0ps')
nod0ps.summary()
nod0ps._index[k]   # 124 rows   - has a very odd order of things, see issue #342

s1 = nod0ps.getps(plnum=0).timeaverage()
s1.stats(qac=True)    # 2.3194959113603586 0.5859667447709659 -1.017372244576539 4.539838041025146
                      
                      
                      
s2 = sdf0.getps(scan=60, plnum=0, ifnum=0, fdnum=0).timeaverage()
s2.stats(qac=True)    # 2.3194959113603586 0.5859667447709659 -1.017372244576539 4.539838041025146
(s2-s1).stats()['rms'].value     # 0.0
#   this is issue #342 when s2-s1 was not zero

sdf3 = GBTFITSLoad(dysh_data(example="nod-KFPA/data/TGBT22A_503_02.raw.vegas.trim.fits"))
sdf3.summary()  # 992 rows
s3 = sdf3.getps(scan=60, plnum=0, ifnum=0, fdnum=0).timeaverage()
s3.stats(qac=True)    # 2.319488925090192 0.5859659084045957 -1.0173834072641188 4.539820576177497
(s2-s3).stats()['rms'].value # 1.716728447608622e-05



s4 = GBTFITSLoad('nod0_test60.fits').getspec(0)
s4.stats(qac=True)    # 2.319489002227783 0.5859659314155579 -1.0173834562301636 4.539820671081543
d = s3.flux.value - s4.flux.value
np.nanstd(d)   # 6.906450791469074e-08
(s4-s3).stats()['rms'].value 

# these were created on LMA machine 
sdf5 = GBTFITSLoad(dysh_data(example="nod-KFPA/data/TGBT22A_503_02.raw.vegas.trim1.fits" ))
sdf5.summary()    # 992 rows
s5 = sdf5.getps(scan=60, plnum=0, ifnum=0, fdnum=0).timeaverage()
s5.stats(qac=True)    #   2.319488925090192 0.5859659084045957 -1.0173834072641188 4.539820576177497
                      #   2.3194959113603586 0.5859667447709659 -1.017372244576539 4.539838041025146  new trim1
(s5-s1).stats()      # with with the new trim1            
                    
sdf6 = GBTFITSLoad(dysh_data(example="nod-KFPA/data/TGBT22A_503_02.raw.vegas.trim7.fits" ))
sdf6.summary()  # 1736 rows
s6 = sdf6.getps(scan=60, plnum=0, ifnum=0, fdnum=0).timeaverage()
s6.stats(qac=True)    #   2.319488925090192 0.5859659084045957 -1.0173834072641188 4.539820576177497
(s5-s6).stats()       #   1.71672845e-05
  
#   ok, s3, s5, s6 are identical !!!

sdf7 = GBTFITSLoad(dysh_data(example="nod-KFPA/data/TGBT22A_503_02.raw.vegas.trim60.fits" ))
sdf7.summary()  # 10416 rows
s7 = sdf7.getps(scan=60, plnum=0, ifnum=0, fdnum=0).timeaverage()
s7.stats(qac=True)    #   2.323344799372547 0.5913854528123186 -0.9301776419005892 4.659711507349211           
                      #   2.319488925090192 0.5859659084045957 -1.0173834072641188 4.539820576177497     NEW
(s7-s1).stats()['rms'].value   # 1.716728447608622e-05

sdf8 = GBTFITSLoad(dysh_data(example="nod-KFPA/data/TGBT22A_503_02.raw.vegas.trim0.fits" ))
sdf8.summary()  # 124 rows
s8 = sdf8.getps(scan=60, plnum=0, ifnum=0, fdnum=0).timeaverage()
s8.stats(qac=True)    #   2.319488925090192 0.5859659084045957 -1.0173834072641188 4.539820576177497
(s8-s1).stats()['rms'].value  # 1.716728447608622e-05



#%%  now testng getnod() on the ones that are supposed to have 62,63

# original
s2n = sdf0.getnod(scan=[62,63], plnum=0, ifnum=0).timeaverage()
s2n.stats(qac=True)   # 0.3748503311988827 0.35504111651677267 -1.2325043436127545 2.1533314871468248

# trim
s3n = sdf3.getnod(scan=[62,63], plnum=0, ifnum=0).timeaverage()
s3n.stats(qac=True)  # 0.3747961289586152 0.35496024341337706 -1.2321275976563082 2.1529178398243145'
(s2n-s3n).stats(qac=True)    # 5.8894550977673705e-05 0.00011608584854222399 -0.00041872595961756076 0.0005777096994614705

# trim1
s5n = sdf5.getnod(scan=[62,63], plnum=0, ifnum=0).timeaverage()
s5n.stats(qac=True)   # 0.3747961289586152 0.35496024341337706 -1.2321275976563082 2.1529178398243145
                      # 0.3748503311988827 0.35504111651677267 -1.2325043436127545 2.1533314871468248  new trim1
# trim7
s6n = sdf6.getnod(scan=[62,63], plnum=0, ifnum=0).timeaverage()
s6n.stats(qac=True)   # 0.3747961289586152 0.35496024341337706 -1.2321275976563082 2.1529178398243145

(s2n-s5n).stats()     # 0.0




#%%   comparing the different getnod styles

nod0g = GBTFITSLoad('nod0_good')
nod0g.summary()
nod0g._index[k]   # 868 rows
getbeam(nod0g)
sp0g = nod0g.getnod().timeaverage()

nod0b = GBTFITSLoad('nod0_bad')               # old version where rows were not sorted
nod0b.summary()
nod0b._index[k]   # 868 rows
getbeam(nod0b)
sp0b = nod0b.getnod().timeaverage()

sp0a = sdf0.getnod(scan=[62,63],plnum=0,ifnum=0).timeaverage()

sp0i = GBTFITSLoad('nod0_test62.fits').getspec(0)    # GBTIDL ground truth  - units are still wrong a 'ct', not 'K'

nod0c = GBTFITSLoad(dysh_data(example="nod-KFPA/data/TGBT22A_503_02.raw.vegas.trim62.fits" ))
nod0c.summary()   # 10416  rows
sp0c = nod0c.getnod(scan=[62,63],plnum=0,ifnum=0).timeaverage()

nod0d = GBTFITSLoad(dysh_data(example="nod-KFPA/data/TGBT22A_503_02.raw.vegas.trim4.fits" ))
nod0d.summary()   # 20832  rows
sp0d = nod0d.getnod(scan=[62,63],plnum=0,ifnum=0).timeaverage()

sp0b.stats(qac=True)     # 0.3815288525997136  0.3599690570433538  -1.214826609304905  2.240987877616527
sp0g.stats(qac=True)     # 0.3748503311988827  0.35504111651677267 -1.2325043436127545 2.1533314871468248
sp0a.stats(qac=True)     # 0.3748503311988827  0.35504111651677267 -1.2325043436127545 2.1533314871468248
sp0i.stats(qac=True)     # 0.37477970123291016 0.35496729612350464 -1.2321275472640991 2.1529178619384766
sp0c.stats(qac=True)     # 0.3747961289586152  0.35496024341337706 -1.2321275976563082 2.1529178398243145
sp0d.stats(qac=True)     # 0.3747961289586152  0.35496024341337706 -1.2321275976563082 2.1529178398243145

(sp0g-sp0b).stats()['rms'].value   # 0.06589345837914665
(sp0g-sp0a).stats()['rms'].value   # 0.0
(sp0g-sp0i).stats()['rms'].value   # 0.0001160868980983191

#%%
# manual mode
sdf62 = GBTFITSLoad('nod0')
sdf62.summary()
sdf62._index[k]    # 2975 for 2 beams   1736 for 7 beams   1*2*31*2*2

# fake gbtidl method
p1a = sdf62.gettp(scan=62, fdnum=2, plnum=0)[0].timeaverage().smooth('box',5)  # ON
p1b = sdf62.gettp(scan=63, fdnum=2, plnum=0)[0].timeaverage().smooth('box',5)  # OFF

t1 = (p1a-p1b)/p1b
p2a = sdf62.gettp(scan=63, fdnum=6, plnum=0)[0].timeaverage().smooth('box',5)  # ON
p2b = sdf62.gettp(scan=62, fdnum=6, plnum=0)[0].timeaverage().smooth('box',5)  # OFF
t2 = (p2a-p2b)/p2b

tsys=68.365   # from nod0
t = tsys*(t1+t2)/2
ts = t.smooth('box',21)
ts.plot(xaxis_unit="km/s", xmin=-100, xmax=50, title='junk62')

#%%

# would this work too?  
# this is like a classic BS mode, and averaging the two BS scans
t3 = (p1a-p2b)/p2b
t4 = (p2a-p1b)/p1b
tsys=68.365   # from nod0
t5 = tsys*(t3+t4)/2
t5.smooth('box',21).plot(xaxis_unit="km/s", xmin=-100, xmax=50, title='getnod "bs" same time')
#  only get a decent piece of spectrum so baseline subtration works
t6 = t5[2000:4500]
t6.baseline(exclude=[750,1750], degree=3, remove=True)
t6.smooth('box',21).plot(xaxis_unit="km/s", xmin=-100, xmax=50, title='getnod "bs" same time')









#%%

#%% EXAMPLE-1   tp_nocal    NOD_BEAMS  10,1   (FEED 11,2)

f1 = dysh_data(accept='AGBT22A_325_15/AGBT22A_325_15.raw.vegas')  # accept='nod1'
sdf=GBTFITSLoad(f1)
# 8 files, 16 beams, each file has 2 beams - 4 scans, VANE/SKY/Nod/Nod
sdf.summary()
# extract 290 and 289 (note order is odd in sdfits:   290 came before 289
mkdir("nod1")
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
mkdir("nod3")
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

mkdir("nod4")
sdf.write('nod4/file.fits', scan=[57,58],plnum=1,overwrite=True)   # 240

nod4 = GBTFITSLoad('nod4')
nod4.summary()
nod4._index[k]   # 240 rows
getbeam(nod4)

sp = getnod(nod4,scan=[57,58])
# tsys=124.63 110.38
sp.smooth('box',20).plot(xaxis_unit='km/s')

# data[int=30,30][pol=2][cal=2] 

# signal expected at 5500 km/ss, but not in spectral range. Big wide feature at -2500
# restffreq at 48.99 GHz  VLSR at 7 km/s

# OK, order of scans is not relevant
scans=[64,65,66,57,58,59,60,61,62,63]     # signal at 0.1K
scans=[57,58,59,60,61,62,63,64,65,66]
scans=[57,58]                             # signal at 0.3K
scans=[59,60]                             # 0.1
scans=[61,62]                             # nothing
scans=[63,64]                             # nothing
scans=[65,66]                             # 0.2
sdf.getnod(scan=scans,plnum=0).timeaverage().smooth('box',151).plot(xaxis_unit="km/s")
# 0.10K signal at -2500
sdf.getnod(scan=scans,plnum=1).timeaverage().smooth('box',151).plot(xaxis_unit="km/s")
# not clear, perhaps some 0.05K feature at -2500 but bad baselines







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

f6 = dysh_data(accept='AGBT17B_319_06/AGBT17B_319_06.raw.vegas')   # AGPeg
sdf = GBTFITSLoad(f6)
# 3 files,  one for each IF
sdf.summary()
getbeam(sdf)


mkdir("nod6")
sdf.write('nod6/file.fits', scan=[9,10], ifnum=0, overwrite=True)   # 248
# scans on AGPeg number 9..54 in nod-pairs

nod6 = GBTFITSLoad('nod6')
nod6.summary()
nod6._index[k]   # 248 rows     2*31*1*2*2
getbeam(nod6)

sp = getnod(nod6, scan=[9,10])
# tsys=133.72 226.86
sp.smooth('box',151).plot(xaxis_unit='km/s')

nod6.getnod()
# Exception: Only one FDNUM is allowed per Scan, found [0, 1]
 
#  data   [if=3][scan=2][int=31][pol=2][cal=2]
#  scan=9 fdnum=0      scan=10 fdnum=1

scans=[9,10]
sp0 = sdf.getnod(scan=scans,plnum=0)
# Exception: Odd number of scans for getnod









#%% EXAMPLE-7  tp   NOD_BEAMS 0,1   (but called FEED=3,7)

f7 = dysh_data(accept='TGBT21A_501_10/TGBT21A_501_10.raw.vegas')   # NGC_2639 
sdf=GBTFITSLoad(f7)               # 2 files
sdf.summary()                     # scans 36..41
getbeam(sdf)                      # [0,1]

mkdir("nod7")
sdf.write('nod7/file.fits', scan=[36,37], plnum=1, overwrite=True)  # 128 selected

nod7 = GBTFITSLoad('nod7')
nod7.summary()
nod7._index[k]  
getbeam(nod7)

sp = getnod(nod7, scan=[36,37])
# tsys=46.83 51.65
sp.smooth('box',151).plot(xaxis_unit='km/s')

#   data   [ scan=2] [int=16] [pol=2] [cal=2]


# more scans (the last scan is too short, so skip 40,41)
scans=[36,37,38,39]
scans=[36,37]    # weak signal ~ 3345
scans=[38,39]    # no obvious signal
sb0 = sdf.getnod(scan=scans,plnum=0)
sb1 = sdf.getnod(scan=scans,plnum=1)
sp=sb0.timeaverage() +  sb1.timeaverage()
sp.smooth('box',151).plot(xaxis_unit="km/s")









#%% EXAMPLE-8            NOD_BEAMS 0,1   (FEED 1,4)
 
f8 = dysh_data(accept="AGBT19A_340_07/AGBT19A_340_07.raw.vegas")   # NGC1377
sdf = GBTFITSLoad(f8)   # 8 files;   4IF 2PL 6-INT 2FEED   (4IF * 2FD in the 8 files)
sdf.summary()           # scans 43..46 NGC1377
sdf._index[k]           # 768 rows
getbeam(sdf)            # [0,1]


mkdir("nod8")
sdf.write('nod8/file.fits', scan=[43,44], ifnum=0, plnum=0, overwrite=True)  # 48


nod8= GBTFITSLoad('nod8')
nod8.summary()
nod8._index[k]   # 48 rows    (2*6*2*2)
getbeam(nod8)

sp = getnod(nod8, scan=[43,44])
# tsys=43.88 51.29
sp.smooth('box',151).plot(xaxis_unit='km/s')
# no obvious line


sb = sdf.getnod(ifnum=0, plnum=0)
sb.timeaverage().smooth('box',51).plot()

sdf.getnod(ifnum=0, plnum=1).timeaverage().smooth('box',151).plot(xaxis_unit="km/s")
# no obvious signal at expected ~1785






#%% EXAMPLE-9     NOD_BEAMS 0,1   (but unknown PROCSCAN)

f9 = dysh_data(accept="AGBT12A_076_05/AGBT12A_076_05.raw.acs")   # NGC7538IRS1
sdf = GBTFITSLoad(f9)
# 1 file
sdf.summary()    # 3072 rows       4*2*12*2*2*8
sdf._index[k]    # scans 10..17    PROCSCAN="unknown"
getbeam(sdf)     # [0,1]
# 4-IF 2-POL 12-INT 2-FEED

mkdir("nod9")
sdf.write('nod9/file.fits', scan=[12,13], ifnum=3, plnum=1, overwrite=True)   # 96

nod9= GBTFITSLoad('nod9')
nod9.summary()
nod9._index[k]     # 96/    768 rows    [scan=2][if=4] [int=12] [pol=2] [cal=2]
getbeam(nod9)
# 2*12*2*2
# very odd order of scans  -  issues/376

sp = getnod(nod9,scan=[12,13])
# tsys=26.48 24.32
sp.smooth('box',41).plot(xaxis_unit='km/s', xmin=-500,xmax=0, ymin=2.55, ymax=2.70)
# signal at -340, expected at -64


nod9.getnod(fdnum=[0,1]).timeaverage().smooth('box',151).plot(xaxis_unit="km/s")
# Exception: Only one FDNUM is allowed per Scan, found [1, 0]

sdf.getnod(fdnum=[0,1], ifnum=0, plnum=1).timeaverage().smooth('box',151).plot(xaxis_unit="km/s")
# Exception: Only one FDNUM is allowed per Scan, found [1, 0]

sb = sdf.getnod(plnum=1, ifnum=3, fdnum=[0,1])
# Exception: Only one FDNUM is allowed per Scan, found [0, 1]
