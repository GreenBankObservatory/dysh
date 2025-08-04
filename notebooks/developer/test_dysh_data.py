#!/usr/bin/env python3

"""
Testing dysh_data  ( a full test is not possible in CI's pytest)\
"""

import os

from dysh.fits.gbtfitsload import GBTFITSLoad
from dysh.fits.gbtfitsload import GBTOnline
from dysh.fits.gbtfitsload import GBTOffline
from dysh.util.files import dysh_data

#  some more liberal panda dataframe display options
import pandas as pd
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 1000)

# pd.options.display.max_columns = None

dd = os.environ["DYSH_DATA"] 
print("Current value of DYSH_DATA from the environment:", dd)

#%%   enable/disable DYSH_DATA

if "DYSH_DATA_ORIG" not in os.environ:
    print("Saving dd_orig")
    os.environ["DYSH_DATA_ORIG"] = os.environ["DYSH_DATA"] 
else:
    print("dd_orig", os.environ["DYSH_DATA_ORIG"])
    
if "DYSH_DATA" in os.environ:  
    print("removing dd")
    os.environ.pop("DYSH_DATA")
else:
    print("no dd")
    
#%% reset

os.environ["DYSH_DATA"] = os.environ["DYSH_DATA_ORIG"] 

#%% SDFITS_DATA

os.environ["SDFITS_DATA"] = '/tmp/sdfits'

#%%  debugging

import dysh
dysh.log.init_logging(3)   # 0=ERROR 1=WARNING 2=INFO 3=DEBUG
dysh.log.init_logging(0)

#%%  helper functions


def size_file(file_path):
    file_size_bytes = os.path.getsize(file_path)
    print(f"The size of file '{file_path}' is: {file_size_bytes/1e6} MB")
    
def size_dir(directory_path):
    total_size = 0
    for dirpath, dirnames, filenames in os.walk(directory_path):
        for f in filenames:
            fp = os.path.join(dirpath, f)
            # Avoid symbolic links to prevent double-counting or infinite loops
            if not os.path.islink(fp):
                total_size += os.path.getsize(fp)
    print(f"The size of dir '{directory_path}' is: {total_size/1e6} MB")


def my_size(fname):
    if os.path.isfile(fname):
        size_file(fname)
    elif os.path.isdir(fname):
        size_dir(fname)
    else:
        print(f"{fname} unkown type")



#%%

# A number of these could easily fail if you don't have your $DYSH_DATA set and populated
# but they should always work at GBO

print(dysh_data(example="nod-KFPA/data/TGBT22A_503_02.raw.vegas.trim.fits"))
print(dysh_data(example="nod-KFPA/data/TGBT22A_503_02.raw.vegas"))
print(dysh_data(example="nod"))
print(dysh_data(example="getps"))
print(dysh_data(example="getpslong"))
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

dysh_data('?')
dysh_data()

#
f = dysh_data(example="nod-KFPA/data/TGBT22A_503_02.raw.vegas.trim.fits")
cmd=f"cp {f} /tmp"
os.system(cmd)
dysh_data("/tmp/TGBT22A_503_02.raw.vegas.trim.fits")



#%%

# file
size_file(filename)

# dir
size_dir(filename)

#
my_size(filename)

#%%


filename = dysh_data(example="getps")
print(filename)
my_size(filename)

filename = dysh_data(example="getpslarge")
print(filename)
my_size(filename)


#%%
filename = dysh_data(example="positionswitch/data/AGBT05B_047_01/AGBT05B_047_01.raw.acs//AGBT05B_047_01.raw.acs.fits")
print(filename)
my_size(filename)
#%%
filename = dysh_data(example="positionswitch/data/AGBT05B_047_01/AGBT05B_047_01.raw.acs/")
print(filename)
my_size(filename)
#%%

filename = dysh_data('AGBT05B_047_01/AGBT05B_047_01.raw.acs')
print('F',filename)
my_size(filename)



#%%

filename = dysh_data('AGBT05B_047_01')
print('F',filename)
my_size(filename)
#%%

sdf = GBTOnline()
filename = sdf.filename
print('F',filename)
my_size(filename)
sdf.summary()

#%%

sdf = GBTOffline('AGBT05B_047_01')
filename = sdf.filename
print('F',filename)
my_size(filename)
sdf.summary()

#%%

filename = dysh_data('AGBT05B_047_01')
print('F',filename)
my_size(filename)
#%%

print(filename)
sdf = GBTFITSLoad(filename)
sdf.summary()

#%%

my_size(dysh_data(test="test1"))        # 46 MB --
my_size(dysh_data(example="test1"))     # 46 MB
my_size(dysh_data(test="getps"))        # 0.55 MB  
my_size(dysh_data(test="getpslarge"))   # 7.8 MB --

#%% i/o bench

# TGBT21A_501_11: pick something
scans = np.arange(152,154)   # NGC2415
scans = np.arange(156,160)   # NGC2782
scans = np.arange(171,197)   # Sco-X

sdf = GBTFITSLoad(dysh_data(example="getpslarge"))              # 7.26 s => 4.2s          2 x 7238 MB 
sb = sdf.getps(scan=scans, ifnum=0, plnum=0, fdnum=0)           # 1min 58s  => 4.5s
sdf.write('NGC2782.fits',scan=scans, overwrite=True, flags=False)           # 6min 38s - 2min 30s  2min 42s  => 3min 48s

sdf1 = GBTFITSLoad("NGC2782.fits")                              # 377 ms       976MB  - linear by size
sb1 = sdf1.getps(scan=scans, ifnum=0, plnum=0, fdnum=0)         # 6.78s   => 425ms    also about linear by size
sdf1.write("NGC2782a.fits", overwrite=True)                         # 13.9 s

#%%
sdf = GBTFITSLoad(dysh_data(example="getpslarge"))   
sdf.getps(scan=[152], ifnum=0, plnum=0, fdnum=0).timeaverage().plot()
sdf.getps(scan=[152,156,158], ifnum=0, plnum=0, fdnum=0).timeaverage().plot()



#%%
