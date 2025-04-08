#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 10:38:12 2024

@author: teuben
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
from dysh.fits.gbtfitsload import GBTOnline
from dysh.fits.gbtfitsload import GBTOffline
from dysh.util.files import dysh_data

from dysh.util import get_project_testdata

#%%  debugging

import dysh
#dysh.log.init_logging(3)   # 0=ERROR 1=WARNING 2=INFO 3=DEBUG
dysh.log.init_logging(0)


#%%  classic getps

f1 = dysh_data(test="getps")        # OnOff, small one (1 IF, 1 PL, 1 INT)
#f1 = dysh_data(example="getps")    # OnOff, long one (4 IFs, 2 PLs, 151 INTs)
#f1 = dysh_data(test="AGBT18B_354_03/AGBT18B_354_03.raw.vegas") # OffOn
fn1 = f1.parts[-1]
print(f"Using {fn1}")


sdf1 = GBTFITSLoad(f1)
sdf1.summary()
sdf1.summary(verbose=True)

#%% trimming for online simulation

sdf=GBTFITSLoad(dysh_data(example='onoff-L/data/TGBT21A_501_11.raw.vegas/TGBT21A_501_11.raw.vegas.B.fits'))
# these are scans 167..196 on ScoX-1

for s in range(168,198,2):
    fn = f"online{s}.fits"
    scans=list(range(167,s+1))
    print(s,scans)
    sdf.write(fn,scan=list(range(167,s+1)), overwrite=True)   

#%% getps reloading?

sdf = GBTOffline(f1)

#%% play with changing files simulation growth

from dysh.fits.gbtfitsload import GBTOnline
from dysh.fits.gbtfitsload import GBTOnline


sdf = GBTOnline("online.fits")
sdf.summary()
sdf._reload()
sdf.summary()

sdf2 = GBTOnline("AGBT21B_024_01")

sdf1 = GBTOnline()
sdf1.summary()



sdf3 = GBTOffline("AGBT21B_024_01")

#%%

sdf =  GBTOnline()

#%%

f1 = dysh_data(example="getps") 
sdf1=GBTFITSLoad(f1)


#%%

os.environ["SDFITS_DATA"] = '/tmp/sdfits'
sdf = GBTOnline()
sdf.summary()
sdf


sdf = GBTOnline('online.fits')


#%%  issue 550

sdf1 = GBTOffline("TRFI_122024_U1")
print(sdf1.filename) # successfully returns file path
/home/sdfits/TRFI_122024_U1/TRFI_122024_U1.raw.vegas

sdf2 = GBTFITSLoad("/home/sdfits/TRFI_122024_U1/TRFI_122024_U1.raw.vegas")
print(sdf2.filename) # fails to return file path
AttributeError: 'GBTSDFITSLoad' object has no attribute '_filename'

# however printing the object works in both cases
print(sdf1)
print(sdf2)

#%% issue 550 my laptop
from dysh.fits.gbtfitsload import GBTFITSLoad
    
f1='AGBT15B_287_33'
sdf1 = GBTOffline(f1)    # ok
print(sdf1.filename) 
# /home/teuben/GBT/dysh_data/sdfits/AGBT15B_287_33/AGBT15B_287_33.raw.vegas

sdf2 = GBTFITSLoad(sdf1.filename)
print(sdf2.filename) 
# AttributeError: 'GBTFITSLoad' object has no attribute '_filename'
print(sdf2.filenames())

sdf3 = GBTFITSLoad(dysh_data(f1))
print(sdf3.filename)



#%% issue 550 with only one file kin the top level

f1 = 'test123'
sdf1 = GBTOffline(f1)    # ok
print(sdf1.filename) 

sdf2 = GBTFITSLoad(sdf1.filename)
print(sdf2.filename) 

sdf3 = GBTFITSLoad('/home/teuben/GBT/dysh_data/sdfits/test123/online168.fits')
print(sdf3.filename) 

sdf4 = SDFITSLoad('/home/teuben/GBT/dysh_data/sdfits/test123/online168.fits')
print(sdf4.filename)


#%%  354


file1 = dysh_data(test='TGBT21A_501_11/TGBT21A_501_11.raw.vegas.fits')
sdf1 = SDFITSLoad(file1)
print(sdf1)     # fails

sdf2= GBTFITSLoad(file1)
print(sdf2)      # works

#%%  20m data
from dysh.fits import gb20mfitsload

f1 = dysh_data(test="20m/Skynet_60476_DR21_118886_68343.cyb.fits")
sdf1 = gb20mfitsload.GB20MFITSLoad(f1)

sdf1.filename


