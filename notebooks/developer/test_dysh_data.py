#!/usr/bin/env python3

"""
Testing dysh_data  ( a full test is not possible in CI's pytest)

"""

import os
from astropy.io import fits
import astropy.units as u

from dysh.fits.sdfitsload import SDFITSLoad
from dysh.fits.gbtfitsload import GBTFITSLoad
from dysh.util.files import dysh_data

#  some more liberal panda dataframe display options
import pandas as pd
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 1000)

# pd.options.display.max_columns = None

dd = os.environ["DYSH_DATA"] 
print("Current value of DYSH_DATA:", dd)


#%%  debugging

import dysh
dysh.log.init_logging(3)   # 0=ERROR 1=WARNING 2=INFO 3=DEBUG

#%%

# A number of these could easily fail if you don't have your $DYSH_DATA set and populated
# but they should always work at GBO

print(dysh_data(example="nod-KFPA/data/TGBT22A_503_02.raw.vegas.trim.fits", verbose=True))
print(dysh_data(example="nod-KFPA/data/TGBT22A_503_02.raw.vegas"))
print(dysh_data(example="nod"))
print(dysh_data(example="getps"))
print(dysh_data(example="getps0"))   # bad url
print(dysh_data(example="getps1"))   # bad url
print(dysh_data(example="positionswitch/data/AGBT05B_047_01/AGBT05B_047_01.raw.acs/"))
print(dysh_data(example="test1"))

print(dysh_data(test="test1"))
print(dysh_data(test="getps"))
print(dysh_data(test="TRCO_230413_Ka/TRCO_230413_Ka_scan43.fits"))

print(dysh_data(accept='AGBT22A_325_15/AGBT22A_325_15.raw.vegas'))
print(dysh_data(accept='nod1'))

print(dysh_data("AGBT21B_024_01"))
print(dysh_data("AGBT21B_024_01/AGBT21B_024_01.raw.vegas"))
print(dysh_data("junk"))   # return None, file does not exist
