#!/usr/bin/env python3

"""
This script was developed in spyder during the vane calibration work.


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


def vcal(sdf, vane, sky, debug=False):
    """ find the two vane and sky and find tsys for each beam
    """
    kb=['DATE-OBS','SCAN', 'IFNUM', 'PLNUM', 'FDNUM', 'PROCSCAN', 'FEED', 'SRFEED', 'FEEDXOFF', 'FEEDEOFF']
    a = sdf._index[kb]
    feeds = a['FDNUM'].unique()
    feeds.sort()
    for f in feeds:
        b = a.loc[a['FDNUM']==f]
        c1 = b.loc[b['SCAN']==vane]
        c2 = b.loc[b['SCAN']==sky]
        print("FEED",f,c1,c2)
        
        # take aver of all the c1's (vane) and c2's (sky)
        c1_data = 1.0
        c2_data = 2.0
    return
    
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




#%% EXAMPLE-1   tp_nocal    NOD_BEAMS  10,1   (FEED 11,2)

f1 = dysh_data(accept='AGBT22A_325_15/AGBT22A_325_15.raw.vegas')  # accept='nod1'
sdf=GBTFITSLoad(f1)
# 8 files, 16 beams, each file has 2 beams - 4 scans, VANE/SKY/Nod/Nod
sdf.summary()    # 256 rows
# extract 290 and 289 (note order is odd in sdfits:   290 came before 289
mkdir("vane1")
sdf.write('vane1/file.fits',scan=[281,282], overwrite=True)   # 64 rows

vane1 = GBTFITSLoad('vane1')
vane1.summary()
vane1._index[k]    # 64 rows

vcal(vane1, 281, 282)


#%%   EDGE data

sdf2 = GBTFITSLoad(dysh_data('AGBT21B_024_01/AGBT21B_024_01.raw.vegas'))
a = sdf2.summary()  # 208 scans
b=a[a["OBJECT"] == "VANE"]
print(b)




