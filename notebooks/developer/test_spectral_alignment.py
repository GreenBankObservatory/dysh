#!/usr/bin/env python3

# fmt: off

import numpy as np
import matplotlib.pyplot as plt


from specutils import Spectrum1D
import astropy.units as u
import numpy as np
from specutils.manipulation import box_smooth, gaussian_smooth, trapezoid_smooth

"""
         0: True,
         4: True,
         5: True,
         6: True,
         7: True, # borderline
         8: True, # Low amp
         9: True,
        10: True,
        12: True,
        14: True,
        15: True,
        16: True, # Multiple lines?
        17: True, # Same frequency range as 16?
        18: True,
        19: True,
        20: True,
        21: True,
        22: True, # Same frequency range as 21?
        24: True,
        26: True,
        27: True,
        28: True,
        29: True, # Same frequency range as 28?
        30: True,
        31: True,
        38: True,
        42: True,
        45: True,
        49: True,
        53: True, # Baseline wobbly
        57: True, # Borderline
        63: True, # Borderline
"""


#%% W43G

f1 = dysh_data(accept='AGBT18B_014_02/AGBT18B_014_02.raw.vegas')
sdf1=GBTFITSLoad(f1)
# 8 files, 
sdf1.summary()

#   SCAN OBJECT VELOCITY   PROC  PROCSEQN  RESTFREQ  DOPFREQ # IF # POL # INT # FEED     AZIMUTH   ELEVATIO
#0     6   W43G    92.02  OffOn         1  5.764969  5.93154   64     2     3      1  117.575564  27.647278
#1     7   W43G    92.02  OffOn         2  5.764969  5.93154   64     2     3      1  117.172716  27.302879


s1 = sdf1.getps(ifnum=16,plnum=0).timeaverage()
s1.plot()
s1.smooth('box',20).plot()


sdf1.getps(ifnum=0,plnum=0).timeaverage().plot()
sdf1.getps(ifnum=4,plnum=0).timeaverage().plot()
sdf1.getps(ifnum=8,plnum=0).timeaverage().plot()     # weak
sdf1.getps(ifnum=21,plnum=0).timeaverage().plot()
sdf1.getps(ifnum=28,plnum=0).timeaverage().plot()
sdf1.getps(ifnum=29,plnum=0).timeaverage().plot()
sdf1.getps(ifnum=53,plnum=0).timeaverage().plot()



#%%

f2 = dysh_data(accept='AGBT18B_014_42/AGBT18B_014_42.raw.vegas')
sdf2=GBTFITSLoad(f2)
# 8 files, 
sdf2.summary()

#   SCAN OBJECT VELOCITY   PROC  PROCSEQN  RESTFREQ  DOPFREQ # IF # POL # INT # FEED     AZIMUTH   ELEVATIO
#0     6   W43G    92.02  OffOn         1  5.764969  5.93154   64     2     3      1  112.105292  22.643071
#1     7   W43G    92.02  OffOn         2  5.764969  5.93154   64     2     3      1   111.73205  22.280439

s2= sdf2.getps(ifnum=16,plnum=0).timeaverage()
s2.plot()
s2.smooth('box',20).plot()
