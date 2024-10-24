#!/usr/bin/env python3

# fmt: off

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


def parr(data, n, w=1):
    """ print  values of an array +/-w around n
    """
    print(data[n-w:n+w+1])
    
    
#%%  testing simple spectra; even though we're not using specutils anymore

spec1 = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm, flux=np.random.default_rng(12345).random(49)*u.Jy)

spec1_bsmooth = box_smooth(spec1, width=3)
spec1_gsmooth = gaussian_smooth(spec1, stddev=3)
spec1_tsmooth = trapezoid_smooth(spec1, width=1)

gaussian_smooth(spec1, stddev=3) 


#%% testing using a PS



f1 = util.get_project_testdata() / 'TGBT21A_501_11/TGBT21A_501_11.raw.vegas.fits'
f1 = dysh_data(test="getps")

sdf1 = GBTFITSLoad(f1)
sdf1.info()
sdf1.summary(verbose=True)

p1 = sdf1.getps(scan=152, ifnum=0, plnum=0)
sp1 = p1[0].calibrated(0)
    
sp1b = sp1.smooth("boxcar",5,-1)
sp1h = sp1.smooth("hanning",1,-1)
sp1g = sp1.smooth("gaussian",5,-1)

d1=sp1.flux.value
d1b=sp1b.flux.value
d1h=sp1h.flux.value
d1g=sp1g.flux.value

# now the smoothref version; it needs to be odd, otherwise
# the IDL kernel isn't the same as astropy

p1s = sdf1.getps(scan=152, ifnum=0, plnum=0, smoothref=15)
sp1s = p1s[0].calibrated(0)
d1s = sp1s.flux.value

np.where(np.isnan(d1))   # 3072

#%% grab GBTIDL data
example1=\
    """
    filein,"TGBT21A_501_11.raw.vegas.fits"
    getps,152,ifnum=0,plnum=0,intnum=0   
    fileout,"TGBT21A_501_11_getps_scan_152_intnum_0_ifnum_0_plnum_0.fits"
    keep
    """

base =  'TGBT21A_501_11/TGBT21A_501_11_getps_scan_152_intnum_0_ifnum_0_plnum_0%s.fits'
base1  = base % ''
base1b = base % '_boxcar_5'
base1h = base % '_hanning'
base1g = base % '_gsmooth_5'
base1s = base % '_smthoff_15'

d2  = fits.open(util.get_project_testdata() / base1)[1].data['DATA'][0]
d2b = fits.open(util.get_project_testdata() / base1b)[1].data['DATA'][0]
d2h = fits.open(util.get_project_testdata() / base1h)[1].data['DATA'][0]
d2g = fits.open(util.get_project_testdata() / base1g)[1].data['DATA'][0]
d2s = fits.open(util.get_project_testdata() / base1s)[1].data['DATA'][0]

e2  = d1-d2
e2b = (d1b-d2b)
e2h = (d1h-d2h)
e2g = (d1g-d2g)
e2s = (d1s-d2s)

print("d1    :", np.nanstd(d1))
print("d1-d2 :", np.nanstd(d1-d2)) 
print("e box :", np.nanstd(e2b))  
print("e han :", np.nanstd(e2h))    
print("e gau :", np.nanstd(e2g))
print("e smth:", np.nanstd(e2s))

#%% issue 415


filename = dysh_data(test="AGBT05B_047_01/AGBT05B_047_01.raw.acs/AGBT05B_047_01.raw.acs.fits")

sdfits = GBTFITSLoad(filename)
sdfits.summary()
sdfits.flag_channel([[170,200],[2880,2980],[31000,32768]])   # note even 33000 cann be used with no warning
scan_block = sdfits.getps(ifnum=0, plnum=0)
ta = scan_block.timeaverage()
ta.plot(xaxis_unit="chan", yaxis_unit="mK", ymin=100, ymax=600, grid=True)
ta.baseline(model="chebyshev", degree=2, exclude=[(14000,18000)], remove=True)
ts = ta.smooth('gaussian',16)
ta.plot(xaxis_unit="chan", yaxis_unit="mK", ymin=100, ymax=600, grid=True)


