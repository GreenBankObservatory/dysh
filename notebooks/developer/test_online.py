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



#%%  debugging

import dysh
dysh.log.init_logging(3)   # 0=ERROR 1=WARNING 2=INFO 3=DEBUG

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


sdf1 = GBTOnline()
sdf1.summary()





