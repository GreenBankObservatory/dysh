#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 15 10:54:10 2025

@author: teuben
"""
import os

from dysh.fits.sdfitsload import SDFITSLoad
from dysh.fits.gbtfitsload import GBTFITSLoad
from dysh.fits.gbtfitsload import GBTOnline
from dysh.fits.gbtfitsload import GBTOffline
from dysh.util.files import dysh_data

#  some more liberal panda dataframe display options
import pandas as pd
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 1000)


#%%  issue 667
#    confusing Calibrating again message - now fixed - see below for the real 295 testing

f1 = dysh_data(example="getps")    # 'sdfits ver1.22' 
sdf1 = GBTFITSLoad(f1)
sdf1.get_summary()
sp1 = sdf1.getspec(0)  #  BUNIT='ct' now after the 295 fix
sp1.flux.unit

sb1 = sdf1.getps(fdnum=0, ifnum=0, plnum=0)
sb1.write('test1.fits',overwrite=True)

sdf2=GBTFITSLoad('test1.fits')
sdf2.get_summary()    # notice order is/was messed up,  known issue
sp2 = sdf2.getspec(0)   # K
sp2.flux.unit

sb2=sdf2.getps(fdnum=0, ifnum=0, plnum=0)   # cannot do, no OFFs anymore
# Exception: OFF scans not found in scan list [51, 53, 55, 57]
sb2=sdf2.gettp(fdnum=0, ifnum=0, plnum=0)
sb2.calibrate()
#Scan 51 was previously calibrated. Calibrating again.
#Scan 53 was previously calibrated. Calibrating again.
#Scan 55 was previously calibrated. Calibrating again.
#Scan 57 was previously calibrated. Calibrating again.

flux1 = sb1[0].calibrated(0)  #  in K
flux2 = sb2[0].calibrated(0)  # in ct
flux1.unit
flux2.unit

flux1-flux2
flux1.data - flux2.data
flux1.data - flux2.data

flux1 = sb1[0].calibrated(0).flux.value
flux2 = sb2[0].calibrated(0).flux.value
(flux1-flux2).std()


sdf1a = SDFITSLoad(f1)
sdf1a['TUNIT7']
sp1a = sdf1a.getspec(0)

sdf1b = SDFITSLoad('test1.fits')
sdf1b['TUNIT7']
sdf1b['BUNIT']    # this is odd
sp1b = sdf1b.getspec(0)

#%%  issue 295     
#    Units in sdfitsloader are not handled


f1 = dysh_data(test="getps")       # 'sdfits ver1.22'   FITSVER=1.9
sdf1 = GBTFITSLoad(f1)
sdf1.get_summary()
sdf1.getspec(0).flux.unit    # now "ct" while 295 fixed, before it was K

sdf1.write("getps0.fits", overwrite=True)   # DATA in col6, FLAGS added


sb1 = sdf1.getps(scan=51, ifnum=0, plnum=0, fdnum=0)
sb1[0].write("getps1.fits", overwrite=True) # 94 columns, data in col94    ScanBase.write()
sp1 = sb1.timeaverage()  
sp1.flux.unit
sp1.plot()
sp1.stats()


sb2 = sdf1.getps(ifnum=0, plnum=0, fdnum=0)
sb2.write("getps2.fits", overwrite=True)     # 94 columns, data in col94 - ScanBlock write

sp2 = sb2.timeaverage()     # WTF, why not working at one point
sp2.flux.unit
sp2.stats()
sp2.write("getps3.fits", overwrite=True)   # fix SINGLE DISH


sdf2 = GBTFITSLoad('getps2.fits')
sdf2.get_summary()
sdf2.getspec(0).flux.unit    # is now 'K', as it should be

#sb2_bad = sdf2.getps(scan=51, ifnum=0, plnum=0, fdnum=0)
sb2_bad = sdf2.gettp(scan=51, ifnum=0, plnum=0, fdnum=0)
sb2_bad[0].timeaverage().flux.unit    # is now 'ct', and that's wrong, but we did something bad  

#%% issue 639
#   scanblock.write() writes dubious SDFITS files 

# covered in the previous cell (getps1, getps2)

#%% old



f1 = dysh_data(example="test1")
sdf1 = GBTFITSLoad(f1)
sdf1.summary()

sb1 = sdf1.getps(fdnum=0, ifnum=0, plnum=0)
sb1.write('test1.fits',overwrite=True)

sdf2=GBTFITSLoad('test1.fits')
sdf2.summary()    # notixe order is messed up,  known issue

sb2=sdf2.getps(fdnum=0, ifnum=0, plnum=0)   # cannot do, no OFFs anymore
# Exception: OFF scans not found in scan list [51, 53, 55, 57]
sb2=sdf2.gettp(fdnum=0, ifnum=0, plnum=0)
sb2.calibrate()
#Scan 51 was previously calibrated. Calibrating again.
#Scan 53 was previously calibrated. Calibrating again.
#Scan 55 was previously calibrated. Calibrating again.
#Scan 57 was previously calibrated. Calibrating again.

flux1 = sb1[0].calibrated(0)
flux2 = sb2[0].calibrated(0)
flux1-flux2
#%%  issue 663
#    Wrong units in getspec

fn=dysh_data(test='AGBT05B_047_01/gbtidl/AGBT05B_047_01.getps.acs.fits')    # 'GBTIDL ver2.10.1' 
sdf=GBTFITSLoad(fn)
sdf['TUNIT7']         #  Ta
sp0=sdf.getspec(0)    # INFO Your data have no units, 'ct' was selected


#%% 596

# I calibrated some data to "Ta*", using bunit="Ta*" and then the plot shows antenna temperature as the y-axis label.
#   595 needs to be fixed first.

# marc: ScanBase needs the bunit attribute split into true FITS bunit and a scale string


#fnm = get_project_testdata() / "TGBT21A_501_11/TGBT21A_501_11.raw.vegas.fits"
fnm = dysh_data(test= "TGBT21A_501_11/TGBT21A_501_11.raw.vegas.fits")
sdf = GBTFITSLoad(fnm)
sba = sdf.getps(scan=152, bunit="ta*", zenith_opacity=0.05, ifnum=0, plnum=0, fdnum=0).timeaverage()
sba.plot()

sba.flux.unit

#%%   707

sdf1 = GBTFITSLoad(dysh_data(test="getps"))
sdf1.nrows(0)
# -> AttributeError: 'GBTFITSLoad' object has no attribute '_nrows'\

# same one for _binheader and _bintable
 
