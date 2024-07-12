#!/usr/bin/env python3

"""
some PS and FS benchmarks using pure numpy, should compare well
to C code

Data dimensions taken from example='getps' and 'getfs'
Each almost 200M data points,  each took about  3sec   (%time in spyder)


"""

import numpy as np
import matplotlib.pyplot as plt



#%% getps  --  dimensions taken from dysh_data(example='getps')

# %time claims 3.0s

nscan=2
ntime=151
nif=5
npol=2
ncal=2
nchan=32768

nedge = nchan//10

data = np.random.normal(1,0.1,(nscan,ntime,nif,npol,ncal,nchan))
print("Size of PS data:",len(np.ravel(data)))

tcal = np.ones( (nscan,ntime,nif,npol,ncal)  )

# raw data:  data[nscan][ntime][nif][npol][ncal][nchan]
# raw tcal:  tcal[nscan][ntime][nif][npol][ncal]

# reduction-1
#     tsys:  tsys[nscan][ntime][nif][npol]


m1 =  data[:,:,:,:,0,nedge:-nedge].mean(axis=4)
m2 = (data[:,:,:,:,1,nedge:-nedge]-data[:,:,:,:,0,nedge:-nedge]).mean(axis=4)
tsys = tcal[:,:,:,:,0] * m1/m2 + 0.5 * tcal[:,:,:,:,0]

on = data[0,:,:,:,0,:]
off = data[1,:,:,:,0,:]

# cheat on tsys by picking one, and reshape to properly broadcast
tsys = tsys[0,:,:,:].reshape(ntime,nif,npol,1)

# Finally, Ta
ta0 = tsys * (on-off)/off

# time aver
ta1 = ta0.mean(axis=0)

# pol aver
ta2 = ta1.mean(axis=1)

# first spectrum
ta = ta2[0]

#%%
plt.plot(ta)


#%% getfs, dimensions taken from dysh_data(example='getfs')

# %time claims 3.3s

nscan = 34
nint  = 11
nif   = 4
npol  = 2
nsig  = 2
ncal  = 2
nchan = 16384


data = np.random.normal(1,0.1,(nscan,nint,nif,npol,nsig,ncal,nchan))
print("Size of FS data:",len(np.ravel(data)))
      
t1 = data[:,:,:,:,:,0,:]
t2 = data[:,:,:,:,:,1,:]
tsys = np.mean(t1,axis=5)/np.mean(t1-t2, axis=5)
data1 = 0.5*(data[:,:,:,:,:,0,:] + data[:,:,:,:,:,1,:])

# @todo  figure out how to do data1 * tsys

sig = data1[:,:,:,:,0,:]
ref = data1[:,:,:,:,1,:]

data2 = sig / (ref-sig)
data3 = np.mean(data2, axis=3)
data4 = np.mean(data3, axis=2)
data5 = np.mean(data4, axis=1)
data6 = np.mean(data5, axis=0)
