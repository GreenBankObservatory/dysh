#!/usr/bin/env python3

"""

Author:   Peter Teuben
Created:  Wed Feb 21 20:50:40 2024
Modified: Wed Feb 21 20:50:40 2024

Description:
------------
This is a script developed in spyder, so it has lots of compute cells,
and lots of data to play with. Mostly with "getfs", but starting
out with a "getps".    Some of the data require you to extract certain
rows from certain fits files, for which we have a script "extract_rows.py"

There are also helper functions "qac()" and "peak()" which may be worthy
of adding to utils.
"""

import numpy as np
import numpy.ma as ma
from scipy.stats import norm
import astropy.units as u
from astropy.io import fits
import matplotlib.pyplot as plt
from dysh.fits.sdfitsload import SDFITSLoad
from dysh.fits.gbtfitsload import GBTFITSLoad
import dysh.util as util
from dysh.util.selection import Selection

def qac(data, label="", diff=True, prec=4, robust=False):
    """ simple regression stats on some data
        reports mean, std, min, max
        @todo   implement prec= and robust=
    """
    if robust:
        print("QAC_STATS: robust not implemented yet")
    else:
        if type(data) == 1:
            print("Need to convert")
            return
        flux = np.nansum(data)
        sump = data[data > 0.0].sum()
        sumn = data[data < 0.0].sum()
        sratio = (sump+sumn)/(sump-sumn)
        
        print('QAC_STATS',label,np.nanmean(data), np.nanstd(data), 
              np.nanmin(data), np.nanmax(data), flux, sratio, len(data))
        if diff:
            # ratio dd/data should be sqrt(2) for unrelated signal
            dd = data[:-1] - data[1:]
            print("QAC_STATS_DIFF", np.nanstd(dd)/np.nanstd(data), np.nanstd(dd),np.nanstd(data))
        
def peak(x,y,imax=-1):
    """ find or force a peak (when imax given)
        then do a 3-point polynomial fit and find the peak location
    """
    if imax < -1:
        raise ValueError("imax should be -1")
            
    if imax < 0:
        imax = y.argmax()
    y1 = y[imax-1]/y[imax]
    y2 = 1.0
    y3 = y[imax+1]/y[imax]
    dx = x[imax]-x[imax-1]
    loc = 0.5*(y1-y3)/(y1+y3-2*y2);
    xpeak = x[imax] + 0.5*(y1-y3)/(y1+y3-2*y2)*dx;
    print('PEAK @%d  %f %f  %g   %g %g' % (imax,x[imax],xpeak,y[imax],y1,y3)) 
 


def cmask(x, sections):
    """
    x = channel list
    sections = list of (c0,c1)


    """
    idx = None
    for i,s in enumerate(sections):
        c0,c1 = s 
        if i == 0:
            #idx = np.where((x>=c0) & (x<=c1))[0]
            m = (x >= c0) & (x <= c1)
        else:
            # idx = np.append(idx, np.where((x>=c0) & (x<=c1))[0])
            m = m |  (x >= c0) & (x <= c1)
    return ma.array(x, mask=m)

def minmatch(strings, s):
    n = len(strings)
    m = []
    for i in range(n):
        if strings[i].find(s) == 0:
            m.append(i)
    if len(m) == 1:
        return strings[m[0]]
    return None
      
def toy_minimum_string_match(s, valid_strings):
    """
    return the valid string given a minimum string input

    Parameters
    ----------
    s : TYPE
        DESCRIPTION.
    valid_strings : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """

    n = len(valid_strings)
    m = []
    for i in range(n):
        if valid_strings[i].find(s) == 0:
            m.append(i)
    if len(m) == 1:
        return valid_strings[m[0]]
    return None

def toggle_sections(nchan, s):
    """
    invert from exclude= to include= or vice versa
    assumes sections are ordered channels

    Parameters
    ----------
    nchan : TYPE
        DESCRIPTION.
    s : list of tuples
        sections that describe, e.g.
        [(100,200),(1000,2000)]

    Returns
    -------
    None.

    """
    ns = len(s)
    s1 = []
    e=0           #  set this to 1 if you want to be exact
    if s[0][0] == 0:
        print("edged")
        for i in range(ns-1):
            s1.append( (s[i][1]+e, s[i+1][0]-e) )    
    else:           
        print("internal")
        s1.append( (0,s[0][0]))
        for i in range(ns-1):
            s1.append( (s[i][1], s[i+1][0]) )
        s1.append( (s[ns-1][1], nchan-1))
    return s1
    
  
   
#%% new PS
   
f1 = util.get_project_testdata() / 'TGBT21A_501_11/TGBT21A_501_11.raw.vegas.fits'
s1 = SDFITSLoad(f1)

sdf1 = GBTFITSLoad(f1)
sdf1.info()
sdf1.summary(verbose=True)

#%%

s = Selection(sdf1)
s.select_range(row=[0,1])
s.show()

#%%

sdf1.select_range(row=[0,1])

p1 = sdf1.getps(scan=152, ifnum=0, plnum=0, calibrate=True, debug=True)

#p1 = sdf1.getps(scan=1, ifnum=0, plnum=0, calibrate=True, debug=True)    # ok
#p1 = sdf1.getps(scan=1, ifnum=0, plnum=1, calibrate=True, debug=True)    # warning
#p1 = sdf1.getps(scan=152, ifnum=0, plnum=1, calibrate=False) # warning
#   note:   calibrate=False still says it calibrated a spectrum

p1.timeaverage(weights=None).plot()

sp1 = p1.timeaverage(weights=None).flux.value
qac(sp1)
# 0.21853874407385052331 0.6796743343920357207 -3.7057502269744873047 4.3438787460327148438 0.0 0.3862831277239638603 32768
#%% new FS

f2 = util.get_project_testdata() / "TGBT21A_504_01/TGBT21A_504_01.raw.vegas/TGBT21A_504_01.raw.vegas.A.fits"
#   the small A.1.fits file only has the first 4 rows, i.e. plnum=1
f2 = util.get_project_testdata() / "TGBT21A_504_01/TGBT21A_504_01.raw.vegas/TGBT21A_504_01.raw.vegas.A.1.fits"
s2 = SDFITSLoad(f2)

sdf2 = GBTFITSLoad(f2)
if False:
    orows = []
    rows = sdf2.scan_rows([20],plnum=1)
    orows.append(np.arange(rows[0],rows[0]+4))
    hdu0  = sdf2._sdf[0]._hdu[0].copy()
    # @todo  this hdu0 isn't valid??? outhdu doesn't compute
    table = sdf2._sdf[0]._hdu[1].data[np.ravel(orows)]
    head  = sdf2._sdf[0]._hdu[1].header
    hdu1  = fits.BinTableHDU(table, header=head)
    outhdu= fits.HDUList([hdu0,hdu1])
    outhdu.writeto("fs_test.fits", overwrite=True)
if False:
    #s = Selection(sdf2)
    # s.select(plnum=1)     # first 4 rows are plnum=1
    sdf2.select_range(row=[0,3])
    # now select the first four "ROW" entries
    # NOTE: not working yet
    
sdf2.info()
sdf2.summary(verbose=True)     # selection applied here right?

# @todo make sure it fails for scan=21
p2 = sdf2.getfs(scan=20, ifnum=0, plnum=1, debug=False, fold=False)
#p2 = sdf2.getfs(scan=20, ifnum=0, plnum=0, debug=True)

p2a = sdf2.getfs(scan=20, ifnum=0, plnum=1, debug=False, fold=False)
#p2 = sdf2.getfs(scan=20, ifnum=0, plnum=1, debug=True, fold=True)
# p2 = sdf2.getfs(scan=20, ifnum=0, plnum=1, debug=True, fold=False)
# p2 = sdf2.getfs(scan=20, ifnum=0, plnum=0, debug=True)

#p2.calibrate()     # ->  print(f'FOLD={kwargs["fold"]}')   'fold' not present

p2.timeaverage(weights=None).plot()

sp1=p2[0].calibrated(0)
sp2=p2.timeaverage()

# peak(sp1.frequency.value, sp1.flux.value)

qac(sp1.flux[10000:13000].value,'sp1_3k')
qac(sp2.flux[10000:13000].value,'sp2_3k')

qac(sp2.flux.value,'sp2_full')
# plnum=0 first int only
# 0.16634342210884511815 6.268664315709928205 -12.429883003234863281 434.8529052734375 0.0 0.22693818093330696847 32768

# plnum=1 all
# 0.16815404608499103038 6.2927746314605035164 -11.961354477599949217 431.08564644969579865 5510.0717821129860834 0.3074105796267763535 32768
# plnum=0 all
# 0.24293653325963295009 10.193067871554597857 -12.8635846108866528596 924.8923895537704029 7960.5443218516525086 0.2987818028412268244 32768

#%%

s2_0 = s2.getspec(0)
s2_1 = s2.getspec(1)
s2_d = s2_1 - s2_0
qac(s2_d.flux.value[5000:10000])
#%%
#    fs_fold = standard
#    fs_fold1 = /linear
#    fs_fold2 = /linear,/wrap
#    fs_fold3 = /spline,/wrap
#    fs_nofold  [perfect match]

f2a = util.get_project_testdata() / "TGBT21A_504_01/TGBT21A_504_01.raw.vegas/fs_fold3.fits"
f2b = util.get_project_testdata() / "TGBT21A_504_01/TGBT21A_504_01.raw.vegas/fs_nofold.fits"
s2a = SDFITSLoad(f2a)
s2b = SDFITSLoad(f2b)
sp2a = s2a.getspec(0).flux.value
sp2b = s2b.getspec(0).flux.value    # /nofold

# only if 4 scans
p2a = sdf2.getfs(scan=20, ifnum=0, plnum=1, debug=False, fold=False)
#

sp1a = p2a[0].calibrated(0).flux.value   # dysh
sp1b = p2b[0].calibrated(0).flux.value

# sp2b
d = sp1a - sp2b


plt.clf()
plt.plot(sp2a)
plt.plot(sp1a)
d = sp2a-sp1a
plt.plot(d)
plt.title('dysh and dysh-gbtidl comparison of folded spectra')
          
qac(d[100:5000])
qac(sp2a[100:5000])
qac(sp1a[100:5000])

qac(d[10000:13000])
qac(sp2a[10000:13000])
qac(sp1a[10000:13000])
#   the nofold spectra are identical !!
plt.plot(sp2b)
plt.plot(sp1b)
qac(sp2b-sp1b)
plt.plot(sp2b)
qac(sp2b[100:5000])
qac(sp2b[10000:13001])

#%% Example of reducing size to one scan one plnum
sp=p2[0].calibrated(0)

#%% grabbing a  GBTIDL calibrated file
pjt0 = 'pjt0.fits'   # plnum=0
s2a = SDFITSLoad(pjt0)
s2a.getspec(0).plot()

# or this
sdf2a = GBTFITSLoad(pjt0)
sdf2a.getspec(0).plot()

#%% GBTIDL example 1: PS

f3 = util.get_project_testdata() / 'ngc5291.fits'
s3 = SDFITSLoad(f3)
s3.info()
sdf3 = GBTFITSLoad(f3)

#%% GBTIDL example 2: FS
f4 = util.get_project_testdata() / 'W3OH.fits'
s4 = SDFITSLoad(f4)
s4.info()

sdf4 = GBTFITSLoad(f4)

#%% GBTIDL example 3: TPnod

f5 = util.get_project_testdata() / 'IC1481.fits'
s5 = SDFITSLoad(f5)
s5.info()
sdf5 = GBTFITSLoad(f5)

#%% my EDGE data: this is Track in an OnOffOn
f6 = util.get_project_testdata() / 'AGBT15B_287_33.raw.vegas/AGBT15B_287_33.raw.vegas.A.fits'
s6 = SDFITSLoad(f6)

sdf6 = GBTFITSLoad(f6)

#%% longer FS data on M33S:   scans 6..39    11840 rows

# f7 = util.get_project_testdata() / 'AGBT20B_014_03.raw.vegas.A6.fits'
# f7 = util.get_project_testdata() / 'AGBT20B_014_03.raw.vegas.A7.fits'
f7 = util.get_project_testdata() / 'AGBT20B_014_03.raw.vegas.A67.fits'
# f7 = util.get_project_testdata() / 'AGBT20B_014_03.raw.vegas/AGBT20B_014_03.raw.vegas.A.fits'
# f7 = util.get_project_testdata() / 'junk.fits'
# f7 = util.get_project_testdata() / 'AGBT20B_014_03.raw.vegas.A1.fits'
# f7 = util.get_project_testdata() / 'AGBT20B_014_03.raw.vegas.A6.fits'


#%%
s7 = SDFITSLoad(f7)
sdf7 = GBTFITSLoad(f7)
sdf7.summary()

if False:
    s = Selection(sdf7)
    # s.select(plnum=1)     # first 4 rows are plnum=1
    s.select_range(row=[0,3])
    # sdf7.select_range(row=[0,3])
    
if False:
    c1=s7.getspec(0)
    c2=s7.getspec(1)
    d12=c2-c1
    d12.plot()
    c1.plot()

#%%

ifnum=0
plnum=0
scans=[6,7]

#p7 = sdf7.getfs(ifnum=ifnum, plnum=plnum, debug=False, fold=True)   # no scans
#p7 = sdf7.getfs(plnum=plnum, debug=False, fold=True)                # no scans
#p7 = sdf7.getfs(scan=scans, plnum=plnum, debug=False, fold=True)    #4 ok
#p7 = sdf7.getfs(scan=scans, ifnum=ifnum, debug=False, fold=True)    #1 should be 2
p7 = sdf7.getfs(scan=scans, ifnum=ifnum, plnum=plnum, debug=False, fold=True, use_sig=True)
# p7 = sdf7.getfs(scan=6, ifnum=0, plnum=1, debug=True,  fold=False)
# p7 = sdf7.getfs(scan=6, ifnum=0, debug=True,  fold=False)

# p7.timeaverage(weights=None).plot()

sp1=p7[0].calibrated(0)
sp2=p7[0].timeaverage()

nchan = len(sp1.flux)

peak(sp2.frequency.value, sp2.flux.value)

if True:
    sp1.plot(xaxis_unit="chan",title=f'plnum={plnum} ifnum={ifnum}  scans={scans}')
    sp2.plot(title='timeaveraged')
else:
    print('no plot')

#%%
s1=(3000,5000)
s2=(7000,9000)
f1 = np.ones(len(f1))
c1 = np.arange(len(f1))
f2 = cmask(c1,[s1,s2])
f1 = ma.array(f1, mask=ma.getmask(f2))
for i in range(len(f1)):
    sp2.flux[i] *= f1[i]
sp2.plot()

#%%     patch spectrum sp1 coords
print("Patching spectrum with a polynomial noise") 
np.random.seed(123)
nchan = len(sp1.flux)
x = np.linspace(0,1,nchan)
sp1._data = np.random.normal(0.0,0.1,nchan)
sp1._data = sp1._data - x/4 + 4*x**2

#%%


model = 'chebyshev'    # will normalize internally 
model = 'legendre'     # will normalize internally 
model = 'hermite'      # will normalize internally 
model = 'polynomial'   # the worst, it needs norm=True

norm = True
kms = u.km/u.s
ghz = u.GHz

sp1.plot(xaxis_unit="km/s",yaxis_unit="K",xmin=-600,xmax=600, ymin=-1,ymax=1)
qac(sp1.flux[6600:7700])

#   the examples with units don't work
if False:
    ex1 = [(0,4200),(5300,6600),(7600,16383)]
    ex1 = toggle_sections(nchan, [(4200,5300),(6600,7600)])
    ex2 = [(-1500*kms,-400*kms),(-200*kms,50*kms),(600*kms,2200*kms)]
    ex3 = [(1.41*ghz,1.42*ghz),(1.422*ghz,1.426*ghz)]
    sp1.baseline(degree=2, exclude=ex1, remove=True, model=model, normalize=norm)
else:
    ex1 = [(4200,5300),(6600,7600)]
    sp1.baseline(degree=2, include=ex1, remove=True, model=model, normalize=norm)

sp1.plot(xaxis_unit="km/s",yaxis_unit="K",xmin=-600,xmax=600, ymin=-1,ymax=1)
qac(sp1.flux[6600:7700])
print(sp1.baseline_model)

#%%

nsp = len(p7)
for i in range(nsp):
    #sp_i = p7[i].timeaverage()
    #qac(sp_i.flux.value)
    sp_i = p7[i].calibrated(0)
    qac(sp_i.flux.value)
print(f"Found {len(p7)} scanblocks")

#%%  classic 2005 W3OH in scans 79..83,  6 integrations each
f8 = util.get_project_testdata() / 'TREG_050627.raw.acs.79.fits'
f8 = util.get_project_testdata() / 'TREG_050627.raw.acs.fits'
s8 = SDFITSLoad(f8)
sdf8 = GBTFITSLoad(f8)
sdf8.summary()
if False:
    # extract scan 79 only
    orows = []
    
    orows.append(np.arange(rows[0],rows[0]+96))
    hdu0  = sdf2._sdf[0]._hdu[0].copy()
    # @todo  this hdu0 isn't valid??? outhdu doesn't compute
    table = sdf2._sdf[0]._hdu[1].data[np.ravel(orows)]
    head  = sdf2._sdf[0]._hdu[1].header
    thdu  = fits.BinTableHDU(table, header=head)
    outhdu= fits.HDUList([hdu0,hdu1])
    outhdu.writeto("fs_test.fits", overwrite=True)
print("%s %s" % (sdf8._index['OBJECT'][0],sdf8._index['DATE-OBS'][0]))

p8 = sdf8.getfs(scan=79, ifnum=0, plnum=1, debug=False, fold=False)
p8.timeaverage(weights=None).plot()


# show the ripple again, recall these older data store calon/off differently
sp8a = s8.getspec(1)
sp8a.plot()
sp8b = s8.getspec(3)
sp8b.plot()
sp8_diff = sp8a-sp8b
sp8_diff.plot()



#%%

dd = util.get_project_testdata() / "TGBT21A_504_01"

f1 = dd / "TGBT21A_504_01.raw.vegas.A.1.fits"
f2 = dd / "TGBT21A_504_01.raw.vegas.A.2.fits"
sdf1 = GBTFITSLoad(f1)
sdf2 = GBTFITSLoad(f2)
sdf1c = sdf1.getfs(scan=20, plnum=1, fold=False)
sdf2c = sdf2.getfs(scan=20, plnum=0, fold=False)
#sdf1c.timeaverage().plot()
#sdf2c.timeaverage().plot()
spf1 = sdf1c.timeaverage().flux.value
spf2 = sdf2c.timeaverage().flux.value    

n1 = dd / "TGBT21A_504_01.raw.vegas.A.1.nofold.fits"
n2 = dd / "TGBT21A_504_01.raw.vegas.A.2.nofold.fits"
n3 = dd / "TGBT21A_504_01.raw.vegas.A.00_nofold.fits"

sdn1 = GBTFITSLoad(n1)
sdn2 = GBTFITSLoad(n2)
sdn3 = GBTFITSLoad(n3)

spn1 = sdn1.getspec(0).flux.value 
spn2 = sdn2.getspec(0).flux.value 
spn3 = sdn3.getspec(0).flux.value 

plt.figure()
plt.clf()
#plt.plot(spn1-spf1)
#plt.plot(spn2-spf2)
plt.plot(spn1)
plt.title(f1)

plt.figure()
plt.clf()
#plt.plot(spn1-spf1)
#plt.plot(spn2-spf2)
plt.plot(spn2)
plt.title(f2)

qac(spn1[3000:5000])





#%%

dd = util.get_project_testdata() / "TGBT21A_504_01/old2"

f1 = dd / "TGBT21A_504_01.raw.vegas.A.11.fits"
f2 = dd / "TGBT21A_504_01.raw.vegas.A.22.fits"
sdf1 = GBTFITSLoad(f1)
sdf2 = GBTFITSLoad(f2)
sdf1c = sdf1.getfs(scan=20, plnum=1, fold=False)
sdf2c = sdf2.getfs(scan=20, plnum=0, fold=False)
#sdf1c.timeaverage().plot()
#sdf2c.timeaverage().plot()
spf1 = sdf1c.timeaverage().flux.value
spf2 = sdf2c.timeaverage().flux.value    

n1 = dd / "TGBT21A_504_01.raw.vegas.A.11_nofold.fits"
n2 = dd / "TGBT21A_504_01.raw.vegas.A.22_nofold.fits"

sdn1 = GBTFITSLoad(n1)
sdn2 = GBTFITSLoad(n2)

spn1 = sdn1.getspec(0).flux.value 
spn2 = sdn2.getspec(0).flux.value 
    

plt.figure()
plt.clf()
#plt.plot(spn1-spf1)
#plt.plot(spn2-spf2)
plt.plot(spn1)
plt.title(f1)

plt.figure()
plt.clf()
#plt.plot(spn1-spf1)
#plt.plot(spn2-spf2)
plt.plot(spn2)
plt.title(f2)

qac(spn1[3000:5000])

# QAC_STATS -4.4523509292574095904e-09 2.7577887449653572036e-08
qac((spn1-spf1)[3000:5000])

#%%



dd = util.get_project_testdata() / "TGBT21A_504_01/"

f1 = dd / "TGBT21A_504_01.raw.vegas.A.10.fits"
n1 = dd / "TGBT21A_504_01.raw.vegas.A.10_nofold.fits"

f1 = dd / "TGBT21A_504_01.raw.vegas"
n1 = dd / "TGBT21A_504_01.vegas.nofold.fits"

sdf1 = GBTFITSLoad(f1)

sdf1c = sdf1.getfs(scan=20, plnum=1, fold=False)
sdf2c = sdf1.getfs(scan=20, plnum=0, fold=False)

spf1 = sdf1c.timeaverage(weights='tsys').flux.value
spf2 = sdf2c.timeaverage(weights='tsys').flux.value    


sdn1 = GBTFITSLoad(n1)

spn1 = sdn1.getspec(0).flux.value 
spn2 = sdn1.getspec(1).flux.value 
    

plt.figure()
plt.clf()
plt.plot(spn1-spf1)
#plt.plot(spn2-spf2)
#plt.plot(spn1)
#plt.plot(spf1)
plt.title(f1)



qac(spn1[3000:5000],'nofold')

# QAC_STATS -4.4523509292574095904e-09 2.7577887449653572036e-08
qac((spn1-spf1)[3000:5000],'diff1')
qac((spn2-spf2)[3000:5000],'diff2')
