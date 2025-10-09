#!/usr/bin/env python3

# fmt: off

import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt

import astropy.units as u
from astropy.io import fits
from specutils import Spectrum1D
from specutils.manipulation import box_smooth, gaussian_smooth, trapezoid_smooth

import matplotlib.pyplot as plt

from dysh.fits.sdfitsload import SDFITSLoad
from dysh.fits.gbtfitsload import GBTFITSLoad
from dysh.fits.gbtfitsload import GBTOffline
import dysh.util as util
from dysh.util.selection import Selection
from dysh.util.files import dysh_data
from dysh.spectra.spectrum import Spectrum

def parr(data, n, w=1):
    """ print  values of an array +/-w around n
    """
    print(data[n-w:n+w+1])
    
def stats(sp):
    """ some more sstats on a spectrum
    """
    n = len(np.where(np.isnan(sp.flux))[0])
    print("nans:",n)
    
    
#%%  testing simple spectra; even though we're not using specutils anymore

spec1 = Spectrum1D(spectral_axis=np.arange(1, 50) * u.nm, flux=np.random.default_rng(12345).random(49)*u.Jy)

spec1_bsmooth = box_smooth(spec1, width=3)
spec1_gsmooth = gaussian_smooth(spec1, stddev=3)
spec1_tsmooth = trapezoid_smooth(spec1, width=1)

gaussian_smooth(spec1, stddev=3) 


#%% testing using a PS



f1 = util.get_project_testdata() / 'TGBT21A_501_11/TGBT21A_501_11.raw.vegas.fits'
# f1 = dysh_data(test="getps")

sdf1 = GBTFITSLoad(f1)
sdf1.info()
sdf1.get_summary(verbose=True)

p1 = sdf1.getps(scan=152, ifnum=0, plnum=0, fdnum=0)
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

p1s = sdf1.getps(scan=152, ifnum=0, plnum=0, fdnum=0, smoothref=15)
sp1s = p1s[0].calibrated(0)
d1s = sp1s.flux.value

np.where(np.isnan(d1))   # 3072
# after ...
# [3069, 3070, 3071, 3072, 3073, 3074, 3075]
    
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

#%% issue 415 : Masks not handled in smoothing
#               The smooth spectrum does not mask the flagged channels.

# from dysh.util import get_project_testdata
# filename = get_project_testdata() / "AGBT05B_047_01/AGBT05B_047_01.raw.acs/AGBT05B_047_01.raw.acs.fits"
filename = dysh_data(test="getps")

sdfits = GBTFITSLoad(filename)
sdfits.get_summary()
ta0 = sdfits.getps(ifnum=0, plnum=0, fdnum=0).timeaverage()
#ta0.plot(xaxis_unit="chan", yaxis_unit="mK", ymin=100, ymax=600, grid=True)
ta0.plot(xaxis_unit="chan",grid=True)

# notice some bad spikes and the upper edge

#      2..2      test
#    130..210    sync type behavior
#   1000..2000   test
#   2860..3010   galactic
#  31000..32768  edge effect
flags1 = [[2,3],[130,200],[1000,2000],[2860,3010],[31000,32768]]
flags2 = [14000,18000]
flags3 = [[130,200],[1000,2000],[2860,3010],[14000,18000],[31000,32768]]
flags3t = list(map(tuple,flags3))

sdfits.flags.clear()
sdfits.flag_channel(flags1)
sdfits.flags.show()
ta1 = sdfits.getps(ifnum=0, plnum=0, fdnum=0).timeaverage()
#ta1.plot(xaxis_unit="chan", yaxis_unit="mK", ymin=100, ymax=600, grid=True)
ta1.plot(xaxis_unit="chan",  grid=True)
# ok, this shows sections it skips where the flags were set


ta1.baseline(model="chebyshev", degree=2, exclude=[(14000,18000)])
   # this will show the baseline solution
   #   doing this first seems to fix the problem in the next one

ta1.baseline(model="chebyshev", degree=2, exclude=[(14000,18000)], remove=True)
# this now shows odd sections where baseline is -0.2 from the old continuum
# as if there are no more masks - should be an issue for Evan?

#ta1.baseline(model="chebyshev", degree=2, exclude=[flags2], remove=True)      #  IndexError: list index out of range
#ta1.plot(xaxis_unit="chan", yaxis_unit="mK", ymin=-200, ymax=300, grid=True)
ta1.plot(xaxis_unit="chan", grid=True)
# now the plot looks good


flags1.append(flags2)
print(flags1)
#ta0.baseline(model="chebyshev", degree=2, exclude=flags3, remove=True)
# ValueError: The truth value of an array with more than one element is ambiguous. Use a.any() or a.all()
# see issue 783
ta0.baseline(model="chebyshev", degree=2, exclude=flags3t, remove=True)



ta0.plot(xaxis_unit="chan", yaxis_unit="mK", grid=True)



ts = ta1.smooth('box',31)
ts.plot(xaxis_unit="chan",   grid=True)
# this looks ok as long as no baseline subtraction was done
# box looks ok 
ts = ta1.smooth('gaussian',31)
ts.plot(xaxis_unit="chan",   grid=True)
# did not skip masked sections


ts = ta1.smooth('hanning',31)
ts.plot(xaxis_unit="chan",   grid=True)
# ok

sdfits.flag(scan=20, 
            channel=[[2300,4096]], 
            intnum=[i for i in range(42,53)])



#%%

a=np.arange(20.0)
a1 = ma.masked_array(a, mask=a>10)
print(a.mean())
print(a1.mean())
a[5] = np.nan
print(a.mean())
print(np.nanmean(a))
print(np.nanmean(a1))

#%%

ta = GBTFITSLoad(dysh_data(test="getps")).getps(ifnum=0, plnum=0, fdnum=0).timeaverage()
print(ta.flux[1])
ta._data[1] = np.nan
print(ta.flux[1])

# note we need _data to patch manually, but _mask and mask are the same.
ta.mask[2] = True
print(ta.flux[2])

ta.mask[2] = False
print(ta.flux[2])git b

#%%


# from dysh.spectra.spectrum import Spectrum
f = Spectrum.fake_spectrum()
f.plot()
f.mask[100:200] = True
f.plot()
print(f.stats()["nan"])   # 100

f1 = f.smooth('gaussian',11)
f1.plot()
print(f1.stats()["nan"])   # 12
# this looks good now

#%%


sdfits = GBTFITSLoad(dysh_data(test="getps"))
ta0 = sdfits.getps(ifnum=0, plnum=0, fdnum=0).timeaverage()
ta0.plot(xaxis_unit="chan", grid=True)

flags3 = [[130,200],[1000,2000],[2860,3010],[14000,18000],[31000,32768]]
flags3t = list(map(tuple,flags3))

ta0.baseline(model="chebyshev", degree=2, exclude=flags3)
    # ValueError: The truth value of an array with more than one element is ambiguous. Use a.any() or a.all()
 
ta0.baseline(model="chebyshev", degree=2, exclude=flags3t)
    # ok
