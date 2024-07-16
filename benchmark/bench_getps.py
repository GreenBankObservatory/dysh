#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

filename = 'example_data/onoff-L/data/TGBT21A_501_11.raw.vegas.fits'
filename = None

#%% getps  --  dimensions taken from dysh_data(example='getps')
%%time

nscan=2
ntime=151
nif=5
npol=2
ncal=2
nchan=32768
nedge = nchan//10
    
# %time claims 3.0s, but 2.5s is the initialization
if filename == None:
    data = np.random.normal(1,0.1,(nscan,ntime,nif,npol,ncal,nchan))
    print("Size of PS data:",len(np.ravel(data)))
    tcal = np.ones( (nscan,ntime,nif,npol,ncal)  )
else:
    from astropy.io import fits
    hdu = fits.open(filename)  
    table = hdu[1].data
    data = table["DATA"]
    print(data.shape)      
    data = data.reshape(nscan,ntime,nif,npol,ncal,nchan)
    tcal = table["TCAL"]
    tcal = tcal.reshape(nscan,ntime,nif,npol,ncal)

#     raw data:  data[nscan][ntime][nif][npol][ncal][nchan]
#     raw tcal:  tcal[nscan][ntime][nif][npol][ncal]
# reductions
#     tsys:  tsys[nscan][ntime][nif][npol]
#     on/off:        ta0[ntime][nif][npol][nchan]
#     after t-aver:  ta1[nif][npol][nchan]
#     after p-aver:  ta2[nif][nchan]
#     pick an IF:    ta[nchan]

#%% bench_ps
%%time
# 0.45sec

m1 =  data[:,:,:,:,0,nedge:-nedge].mean(axis=4)
m2 = (data[:,:,:,:,1,nedge:-nedge]-data[:,:,:,:,0,nedge:-nedge]).mean(axis=4)
tsys = tcal[:,:,:,:,0] * m1/m2 + 0.5 * tcal[:,:,:,:,0]
on  = data[0,:,:,:,0,:]
off = data[1,:,:,:,0,:]
# cheat on tsys by picking one, and reshape to properly broadcast
tsys = tsys[0,:,:,:].reshape(ntime,nif,npol,1)
# Finally, Ta0
ta0 = tsys * (on-off)/off
# time aver
ta1 = ta0.mean(axis=0)
# pol aver
ta2 = ta1.mean(axis=1)
# first IF spectrum
ta = ta2[0]

#%%

from scipy.signal import convolve

box = np.ones(64)/64.0
tas = convolve(ta,box)
plt.plot(ta)
plt.plot(tas)