#!/usr/bin/env python3
#

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
kw=['DATE-OBS','SCAN', 'SUBOBSMODE', 'IFNUM', 'PLNUM', 'FDNUM', 'INTNUM', 'CAL', 'PROCSEQN', 'CALPOSITION','TCAL','TSYS']

#  some more liberal panda dataframe display options
import pandas as pd
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 1000)

# pd.options.display.max_columns = None

#%%  debugging

# ?? can one only do this once per kernel ??

import dysh
#dysh.log.init_logging(3)   # 0=ERROR 1=WARNING 2=INFO 3=DEBUG
dysh.log.init_logging(1)    # 0 didn't work


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
            
def regress(data_string, accuracy = 1e-9):
    """ what qac_stats does
    """
    s = data_string.split()
    ns = len(s)
    

#%%   full edge

f1=dysh_data('AGBT15B_287/AGBT15B_287_19.raw.vegas')
sdf = GBTFITSLoad(f1)
sdf.summary()


if False:
    # make a small test, HI in the 2nd IF
    sdf.write('edge.fits',scan=[56,57,58], ifnum=1, plnum=0, intnum=0, multifile=False, overwrite=True)

# the full data should have 23x better S/N
sp = []
for s in [56, 59, 62]:         # three triplets
   print("Working on ",s)
   for pl in [0,1]:            # two polarizations
       sp1 = sdf.gettp(scan=s+0,ifnum=1,fdnum=0,plnum=pl).timeaverage()
       sp2 = sdf.gettp(scan=s+1,ifnum=1,fdnum=0,plnum=pl).timeaverage()
       sp3 = sdf.gettp(scan=s+2,ifnum=1,fdnum=0,plnum=pl).timeaverage()
       sp.append((sp1-sp2)/sp2 * sp1.meta["TSYS"])
       sp.append((sp3-sp2)/sp2 * sp3.meta["TSYS"])

final_sp = sp[0].average(sp[1:])
final_sp[20000:30000].stats(qac=True)
#  '0.020347771462288694 0.009417386409339035 -0.25814777126242 0.19478180068687698'     2-pol
final_sp.plot(xaxis_unit="km/s", xmin=-4000, xmax=-3500, ymin=-0.1, ymax=1.6)

#%% small test needs edge.fits

sdf = GBTFITSLoad('edge.fits')
sdf.summary()

sp=[]
s=56
sp1 = sdf.gettp(scan=s+0,ifnum=1,fdnum=0,plnum=0).timeaverage()
sp2 = sdf.gettp(scan=s+1,ifnum=1,fdnum=0,plnum=0).timeaverage()
sp3 = sdf.gettp(scan=s+2,ifnum=1,fdnum=0,plnum=0).timeaverage()
p1 = (sp1-sp2)/sp2 * sp1.meta["TSYS"]
p2 = (sp3-sp2)/sp2 * sp3.meta["TSYS"]
final_sp = p1.average(p2)
final_sp[20000:30000].stats(qac=True)
# '0.09221439817177293 0.15920341703887056 -0.5727973574918991 0.709867703916811'

final_sp.meta["RESTFREQ"] = final_sp.meta["RESTFRQ"] = 1420.405751786

final_sp.plot(xaxis_unit="km/s")
# spectrum is at ~-3720,   vlsr ~ 1720 though,   restfreq wrong?

final_sp.plot(xaxis_unit="km/s", xmin=-4000, xmax=-3500, ymin=-0.3, ymax=1.9)


#%% new getsigref from #546

p1 = sdf.getsigref(scan=56,ref=57,fdnum=0,ifnum=1,plnum=0).timeaverage()
p2 = sdf.getsigref(scan=58,ref=57,fdnum=0,ifnum=1,plnum=0).timeaverage()
final_sp2 = p1.average(p2)
final_sp2[20000:30000].stats(qac=True)
# '0.09244847050922786 0.15949706831348484 -0.5741176754236221 0.7108851075172424'

#%% why can't I subtract these two.

d = final_sp - final_sp2
# UnitConversionError: Can only apply 'subtract' function to quantities with compatible dimensions

# when drilling down to sp.flux.data they also cannot be subtracted

# TypeError: unsupported operand type(s) for -: 'memoryview' and 'memoryview'
