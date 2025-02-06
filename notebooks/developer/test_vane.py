#!/usr/bin/env python3

"""
This script was developed in spyder during the vane calibration work. 

Argus and the W (3mm) band use a vane to determine Tsys
Q (7mm) band has a vane as well as noise diodes.   Cross check?



The following examples are covered here:

From nodding we have two cases:
getps - review the standard test='getps'
nod1  - VANE/SKY     fdnum=[10,1]           tp nocal;  TBD
edge1 - VANE/SKY
edge2 - VANE/SKY   at 112,114 GHz
nod3  - CALSEQ         fdnum=[0,1]            tp nocal;  TBD
vane4 - ?
vane5 - VANE/SKY
vane6 - CALSEQ    at 4 freqs

Related issues:

https://github.com/GreenBankObservatory/dysh/issues/457     flagging problem

for Peter:    cd ~/GBT/dysh_data/nodding
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
from dysh.spectra.core import mean_tsys
    
#  useful keys for a mult-beam observation listing

k=['DATE-OBS','SCAN', 'IFNUM', 'PLNUM', 'FDNUM', 'INTNUM', 'PROCSCAN','FEED', 'SRFEED', 'FEEDXOFF', 'FEEDEOFF', 'SIG', 'CAL', 'PROCSEQN', 'PROCSIZE']
ks=['DATE-OBS','SCAN', 'SUBOBSMODE', 'IFNUM', 'PLNUM', 'FDNUM', 'INTNUM', 'CAL', 'PROCSEQN']
kw=['DATE-OBS','SCAN', 'SUBOBSMODE', 'IFNUM', 'PLNUM', 'FDNUM', 'INTNUM', 'CAL', 'PROCSEQN', 'CALPOSITION']

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

def mkdir(name, clean=True):
    """ simpler frontend for making a directory that might also already exist
        clean=True:    also remove files inside
    """
    os.makedirs(name, exist_ok = True)
    if clean:
        fns = os.listdir(name)
        for fn in fns:
            print(f"Removing {fn} from {name}")
            os.unlink(os.path.join(name,fn))


def vcal(sdf, vane, sky, debug=False):
    """ find the two vane and sky and find tsys for each beam
    """
    kb=['DATE-OBS','SCAN', 'IFNUM', 'PLNUM', 'FDNUM', 'PROCSCAN', 'FEED', 'SRFEED', 'FEEDXOFF', 'FEEDEOFF']
    a = sdf._index[kb]
    feeds = a['FDNUM'].unique()
    feeds.sort()
    if debug:
        print("Found feeds", feeds)
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

def getbeam(sdf, debug=False):
    """ find the two nodding beams based on FDNUM, FEED 
        needs PROCSCAN='BEAM1' or 'BEAM2'
    """
    kb=['DATE-OBS','SCAN', 'IFNUM', 'PLNUM', 'FDNUM', 'PROCSCAN', 'FEED', 'SRFEED', 'FEEDXOFF', 'FEEDEOFF']
    a = sdf._index[kb]
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
    
def mean_data(data, fedge=0.1, nedge=None):
    """ special mean to exclude the edges like mean_tsys()
    """
    nchan = len(data)
    if nedge is None:
        nedge = int(nchan * fedge)
    chrng = slice(nedge, -(nedge - 1), 1)
    meandata = np.nanmean(data[chrng])
    return meandata
    
def vanecal(sdf, vane, sky,  feeds=[], tcal=None, tmpfile=False):
    """ loop over feeds to get tsys factor
        for efficiency of large data, it's better to sdf.write() and use only the
        vane and sky scans'
        Example EDGE:  (1.3GB)
             all 163 scans:     6m33s    3m6s
             write small sdf:     25s    
             just vane/cal:       20s
    """
    if len(feeds) == 0:
        print("Warning, no feeds= given")
        return None
    tsys = np.zeros(len(feeds), dtype=float)
    if tmpfile:
        # needs testing, what to do about ifnum, plnum - force one
        print(f"Writing scans {vane} and {sky} to _vanecal")
        mkdir('_vanecal')
        sdf.write('_vanecal/file.fits',scan=[vane,sky],overwrite=True)
        sdf1 = GBTFITSLoad('_vanecal')
    else:
        sdf1 = sdf
       
    #  for VANE/CAL data usually tcal=1
    if tcal is None:
        tcal = sdf._index['TCAL'].mean()

    
    i = 0
    for f in feeds:
        v = sdf1.gettp(scan=vane, fdnum=f, calibrate=True, cal=False).timeaverage()
        s = sdf1.gettp(scan=sky,  fdnum=f, calibrate=True, cal=False).timeaverage()
        if True:
            mean_off = mean_data(s.data)
            mean_dif = mean_data(v.data - s.data)
            tsys[i] = tcal * mean_dif/mean_off
        else:   
            tsys[i] = mean_tsys(v.flux, s.flux, 1.0) * tcal    # wrong for vane/cal
            # optional?    += tcal/2.0
        i = i + 1
    print("TCAL=",tcal)
     
    return tsys

    
def calseq(sdf, scan, tcold=54, ifnum=0, plnum=0, fdnum=0, freq=None):
    """ W-band receivers use a CALSEQ
        This routine returns the gain and Tsys for W-band channel
        
        Tcold = 54 - 0.6*(FREQ-77)      FREQ in GHz
    """
    if freq is not None:
        # see eq.(13) in GBT memo 302
        tcold = 54 - 0.6*(freq-77)
        print(f"Warning: calseq using freq={freq} GHz and setting tcold={tcold} K")
        
    twarm = sdf._index['TAMBIENT'].mean()
        
    if False:
        a = sdf._index[['CALPOSITION']]
        o = list(a.loc[a['CALPOSITION']=='Observing'].index)
        c1 = list(a.loc[a['CALPOSITION']=='Cold1'].index)
        c2 = list(a.loc[a['CALPOSITION']=='Cold2'].index)
        vsky = sdf.gettp(scan=scan,ifnum=ifnum,plnum=plnum,fdnum=fdnum,intnum=o,calibrate=True,cal=False).timeaverage()
        vcold1  = sdf.gettp(scan=scan,ifnum=ifnum,plnum=plnum,fdnum=fdnum,intnum=c1,calibrate=True,cal=False).timeaverage()
        vcold2  = sdf.gettp(scan=scan,ifnum=ifnum,plnum=plnum,fdnum=fdnum,intnum=c2,calibrate=True,cal=False).timeaverage()
    else:
        tp_args = {"scan":scan,"ifnum":ifnum,"plnum":plnum,"fdnum":fdnum,"calibrate":True,"cal":False}
        vsky = sdf.gettp(CALPOSITION="Observing", **tp_args).timeaverage()
        vcold1  = sdf.gettp(CALPOSITION="Cold1", **tp_args).timeaverage()
        vcold2  = sdf.gettp(CALPOSITION="Cold2", **tp_args).timeaverage()
    
    if fdnum == 0:
        g = (twarm-tcold)/mean_data(vcold2.data-vcold1.data)
    elif fdnum == 1:
        g = (twarm-tcold)/mean_data(vcold1.data-vcold2.data)
    else:   
        print(f"Illegal fdnum={fdnum} for a CALSEQ")
        return None
    tsys = mean_data(g*vsky.data)
   
    print(f"Twarm={twarm} Tcold={tcold}")
    print(f"IFNUM {ifnum} PLNUM {plnum} FDNUM {fdnum}")
    print(f"Tsys = {tsys}")
    print(f"Gain [K/counts] = {g}")
    return tsys, g

# scan = auto calseq scan
# tcold = effective temperature of cold load (e.g., 50K)
# ifnum = IFnum of spectral window
# plnum = pol-number
# fdnum = beam-number

"""
;;Output:
;;Prints Tsys and gain and returns OUTgain
;;OUTgain= gain = (Twarm-Tcold)/(warm-cold) [K/counts]
;;Tsys=median(gain*sky)

if (n_elements(ifnum) eq 0) then ifnum = 0
if (n_elements(fdnum) eq 0) then fdnum = 0
if (n_elements(plnum) eq 0) then plnum = 0
if (n_elements(tcold) eq 0) then tcold = 54.

# CALPOSITION contains the wcalpos

gettp,scan,plnum=plnum,fdnum=fdnum,ifnum=ifnum,quiet=1,wcalpos='Observing'
vsky=getdata(0)
twarm=!g.s[0].twarm
gettp,scan,plnum=plnum,fdnum=fdnum,ifnum=ifnum,quiet=1,wcalpos='Cold1'
vcold1=getdata(0)
gettp,scan,plnum=plnum,fdnum=fdnum,ifnum=ifnum,quiet=1,wcalpos='Cold2'
vcold2=getdata(0)


;;Feed =1 or 2 for the two possible receiver beams
feed=!g.s[0].feed
gain=0.0
if (feed eq 1) then gain=(twarm-tcold)/median(vcold2-vcold1)
if (feed eq 2) then gain=(twarm-tcold)/median(vcold1-vcold2)
tsys=median(gain*vsky)

print,'Twarm, Tcold:',twarm,tcold
print,'IFNUM, FDNUM, PLNUM:',ifnum,fdnum,plnum
print,'Tsys =',tsys
print,'Gain [K/counts]=',gain
OUTgain=gain

"""

def getnod(sdf, scan1, scan2, beam1, beam2):
    ps1_on = sdf.gettp(scan=scan1, fdnum=beam1, calibrate=True, cal=False).timeaverage()
    ps1_off = sdf.gettp(scan=scan2, fdnum=beam1, calibrate=True, cal=False).timeaverage()
    sp1 = (ps1_on - ps1_off)/ps1_off

    ps2_on = sdf.gettp(scan=scan2, fdnum=beam2, calibrate=True, cal=False).timeaverage()
    ps2_off = sdf.gettp(scan=scan1, fdnum=beam2, calibrate=True, cal=False).timeaverage()
    sp2 = (ps2_on - ps2_off)/ps2_off

    return (sp1,sp2)

#%% calc_etamb

def calc_etamb(freq, Jupiter=False):

    """
    For a given input frequency, calculate and return the eta_mb value
    for the GBT using David's equation from GBT Memo 302. 

    The equation is not included in the memo and provided by private 
    communication:
        eta_jupiter = 1.23 * eta_aperture + 0.005 * (nu - 60) - 0.00003*(nu-60)**2
    where nu is in GHz.

    It comes from a polynominal fit to nu and eta_jupiter using eta_aperture
    as a function of frequency.

    This correction assumes that our sources are approximately the
    size of Jupiter (43" diameter), which isn't a bad assumption for
    extended molecular gas.

    """

    import math
    from astropy import constants as c

    # surface error for GBT with optimal surface and excellent weather
    esurf = 0.0235*u.cm # cm; GBT memo 302 says 230micron = 0.0230

    # make the input frequency into a quantity if it isn't already
    if not isinstance(freq,u.Quantity):
        if freq > 1.0e9: # freq likely in Hz
            freq = freq * u.Hz
        elif freq < 200.0: # freq likely in GHz
            freq = freq * u.GHz
        else:
            print("check units on input frequency. Should be either in GHz or Hz")
    # convert frequency to GHz to use in fit
    freq = freq.to(u.GHz)

    # calculate equivalent wavelength
    wave = freq.to(u.cm, equivalencies=u.spectral()) 

    # aperture efficiency for GBT using Ruze equation.
    # 0.71 is the aperture efficiency of the GBT at lower frequency.
    eta_a = 0.71 * math.exp(- ( 4 * math.pi * esurf / wave)**2)

    if Jupiter:
        # calculate eta_mb for a jupiter sized source via 
        # polynominal fit from GBT Memo 302.  The equation is not
        # included in the memo and provided by private communication:
        #        eta_jupiter = 1.23 * eta_aperture + 0.005 * (nu - 60) -
        #        0.00003*(nu-60)**2 
        # where nu is in GHz. This correction assumes that our sources 
        # are approximately the size of Jupiter (43" diameter)
        eta_mb = 1.23 * eta_a + 0.005*(freq.value-60) - 0.00003 * (freq.value - 60)**2
    # else calculate "small source" eta_mb.
    elif freq > 110.0*u.GHz:
        # corrected for higher freq (D.Frayer) 
        eta_mb = 1.24 * (freq.value/113)**0.8 * math.exp ( -(1.3*freq.value/113)**2 )
    elif freq > 100.0*u.GHz:
        # GBT memo 302 finds that at high frequencies the eta_mb/eta_b ratio is more like 1.45 due to a slightly larger beam size factor (1.28 instead of 1.2).
        eta_mb = 1.45 * eta_a
    else: 
        # GBT memo 302 finds that for frequencies of 86-90GHz you can use the expected theoretical ratio of 1.274
        eta_mb = 1.274 * eta_a

    print("ETA A,MB: ",eta_a, eta_mb, "at freq ", freq)
        
    return (eta_a, eta_mb)



#%%  classic tp/ps

f1 = dysh_data(test="getps")        # OnOff, small one (1 IF, 1 PL, 1 INT)
#f1 = dysh_data(example="getps")    # OnOff, long one (4 IFs, 2 PLs, 151 INTs)
#f1 = dysh_data(test="AGBT18B_354_03/AGBT18B_354_03.raw.vegas") # OffOn
fn1 = f1.parts[-1]
print(f"Using {fn1}")


sdf1 = GBTFITSLoad(f1)
sdf1.summary()
sdf1.summary(verbose=True)
sdf1._index[k]

sp1 = sdf1.getps().timeaverage()
sp1.plot(title='getps: %s' % fn1)
sp1.stats(qac=True)
# 0.21853874407385052 0.6796743343920357 -3.7057502269744873 4.343878746032715

tp1a = sdf1.gettp(scan=152).timeaverage()
tp1b = sdf1.gettp(scan=153).timeaverage()

tp1a.plot(title='gettp scan=152')
tp1b.plot(title='gettp scan=153')

tp1a.stats(qac=True)  # '490015338.1141468 152969155.2487644 3311173.5 805602304.0'

#%% NOD EXAMPLE-1   tp_nocal    NOD_BEAMS  10,1   (FEED 11,2)

_help = """
   SCAN    OBJECT VELOCITY   PROC  PROCSEQN    RESTFREQ     DOPFREQ # IF # POL # INT # FEED     AZIMUTH   ELEVATIO
0   281      VANE      0.0  Track         1  111.711282  111.711282    1     1     2     16  310.707976  56.549105
1   282       SKY      0.0  Track         1  111.711282  111.711282    1     1     2     16  310.708028  56.549066
2   290  1-631680      0.0    Nod         2  111.711282  111.711282    1     1     6     16    310.3087  55.619947
3   289  1-631680      0.0    Nod         1  111.711282  111.711282    1     1     6     16  310.344245  55.705345
"""


f1 = dysh_data(accept='AGBT22A_325_15/AGBT22A_325_15.raw.vegas')  # accept='nod1'
sdf1=GBTFITSLoad(f1, skipflags=True)
# 8 files, 16 beams, each file has 2 beams - 4 scans, VANE/SKY/Nod/Nod
sdf1.summary()    # 256 rows
# extract 290 and 289 (note order is odd in sdfits:   290 came before 289
mkdir("vane1")
sdf1.write('vane1/file.fits',scan=[281,282], overwrite=True)   # 64 rows

vane1 = GBTFITSLoad('vane1')
vane1.summary()
vane1._index[ks]    # 64 rows

vcal(vane1, 281, 282)

v1 = vane1.getspec(0).flux.value
s1 = vane1.getspec(4).flux.value
t1 = s1/(v1-s1)                          # 0.53 +/0 0.02
np.nanmean(t1[100:900])  # 0.52613854
np.nanstd(t1[100:900])   # 0.016047847

v1 = vane1.gettp(scan=281, fdnum=8, calibrate=True, cal=False).timeaverage()
# ok   Message: 'fitsindex='
s1 = vane1.gettp(scan=282, fdnum=8, calibrate=True, cal=False).timeaverage()
# ok   Message: 'fitsindex='

t1 = s1/(v1-s1) 
t1 = t1[100:900]
# UnitTypeError: Longitude instances require units equivalent to 'rad', so cannot set it to ''


v1.plot(title='VANE')
s1.plot(title='SKY')

print(vanecal(vane1, 281, 282, feeds=range(16)))
# 1.39763406 1.73117692 0.27474769 1.51735319 2.0114504  1.60521632
# 2.03984903 1.94033088 1.92176731 2.05082117 1.98364723 1.7910409
# 2.01003515 1.89671295 1.98234    1.99136881]

sp1,sp2 = getnod(sdf1, 289, 290, 1, 10)
sp1.plot()
sp2.plot()

sp3 = 0.5*(sp1+sp2)
sp3.plot()


#%%   EDGE1 data  (1363MB)

_help = """
     SCAN    OBJECT VELOCITY       PROC  PROCSEQN    RESTFREQ     DOPFREQ # IF # POL # INT # FEED     AZIMUTH   ELEVATIO
0      17      VANE      0.0      Track         1  113.571858  113.571858    1     1    21     16  165.049197  79.052319
1      18       SKY      0.0      Track         1  113.571858  113.571858    1     1    21     16  165.050181  79.052262
2      19   NGC0001      0.0        Nod         1  113.571858  113.571858    1     1   122     16  166.640475   79.14949
3      20   NGC0001      0.0        Nod         2  113.571858  113.571858    1     1   122     16  168.125608  79.199539
4      21      VANE      0.0      Track         1  113.571858  113.571858    1     1    21     16  178.544679  79.395538
5      22       SKY      0.0      Track         1  113.571858  113.571858    1     1    21     16  178.545661  79.395501
6      23   NGC0001      0.0  DecLatMap         1  113.571858  113.571858    1     1    68     16  177.410654  78.233486
"""

sdf2 = GBTFITSLoad(dysh_data('AGBT21B_024_01/AGBT21B_024_01.raw.vegas'), skipflags=True)
a = sdf2.summary()  # 208 scans   1363MB
print(a)
b=a[a["OBJECT"] == "VANE"]
print(b)   # there are 13 vane's in here

# make a smaller dataset for testing
mkdir("edge1")
# sdf2.write("edge1/file.fits", fdnum=8, scan=[17,18,19,20], intnum=[0,1], overwrite=True)
# sdf2.write("edge1/file.fits", fdnum=8, scan=[17,18,19,20], overwrite=True)
sdf2.write("edge1/file.fits", scan=[17,18,19,20], overwrite=True)                # 4576 rows
#  CPU times: user 23.5 s, sys: 3.02 s, total: 26.5 s   Wall time: 25.2 s

edge1 = GBTFITSLoad("edge1")
edge1.summary()
edge1._index[ks]

tsys = vanecal(edge1, 17, 18, feeds=range(16))
#  sdf2   CPU times: user 6min 29s, sys: 10.5 s, total: 6min 40s   Wall time: 6min 33s
#  edge1  CPU times: user   20.1 s, sys: 271 ms, total:   20.4 s   Wall time:   20.1 s
print(tsys)

# 17,18
# [1.3001234  1.55490619 1.57251803 1.51060869 1.63612382 1.42200443
#  1.59264637 1.58702561 1.65487421 1.60720362 1.52214009 1.52169199
#  1.72276893 1.50905505 1.83133043 1.68868181]

# 21,22
# [1.28974196, 1.53318127, 1.55064672, 1.50331784, 1.61864391,
#  1.37614122, 1.56854187, 1.56631664, 1.63607151, 1.56983139,
#  1.4890203 , 1.50031436, 1.71650273, 1.47633432, 1.80116151,
#  1.679226  ])

# plotting a passband
edge1.gettp(scan=17, fdnum=8, calibrate=True, cal=False).timeaverage().plot()
# ok !!

# not working yet, since it wants ON/OFF, needs the GETTP hack
sp2 = sdf2.getnod(scan=19)
# IndexError: index 0 is out of bounds for axis 0 with size 0


# if getnod() is not converted, we can use the old trick of doing two gettp()
scan1 = 19
scan2 = 20
beam1 = 4
beam2 = 7
sp1,sp2 = getnod(edge1, scan1, scan2, beam1, beam2)

sp3 = 0.5*(sp1+sp2)
sp3.plot(xaxis_unit="km/s")
# center is at 0 km/s, that seems wrong, should be VLSR 4500 km/s at 113.57


#%% EDGE2 nod example on NGC5908
sdf1=GBTFITSLoad(dysh_data('AGBT21B_024_14/AGBT21B_024_14.raw.vegas'), skipflags=True)

mkdir("edge2")
sdf1.write("edge2/file.fits", scan=[327,328,329,330,331,332,333,334], overwrite=True)                # 4576 rows
#  CPU times: user 23.5 s, sys: 3.02 s, total: 26.5 s   Wall time: 25.2 s

edge2 = GBTFITSLoad("edge2")
edge2.summary()
edge2._index[ks] # 5248 rows

_s = """
   SCAN   OBJECT VELOCITY   PROC  PROCSEQN    RESTFREQ     DOPFREQ # IF # POL # INT # FEED     AZIMUTH   ELEVATIO
0   327     VANE      0.0  Track         1  112.365952  112.365952    1     1    21     16  314.014676  69.054141
1   328      SKY      0.0  Track         1  112.365952  112.365952    1     1    21     16  314.015393  69.054142
2   329     VANE      0.0  Track         1   114.03043   114.03043    1     1    21     16  333.855801  70.214517
3   330      SKY      0.0  Track         1   114.03043   114.03043    1     1    21     16  333.856727  70.214422
4   331  NGC5908      0.0    Nod         1   114.03043   114.03043    1     1    61     16  333.089151  70.104129
5   332  NGC5908      0.0    Nod         2   114.03043   114.03043    1     1    61     16  332.870869   70.02778
6   333  NGC5908      0.0    Nod         1   114.03043   114.03043    1     1    61     16  332.654981  69.984783
7   334  NGC5908      0.0    Nod         2   114.03043   114.03043    1     1    61     16   332.44083  69.907597
"""
tsys1 = vanecal(edge2, 327, 328, feeds=range(16))
print(tsys1)
# @112.3 GHz

#[1.41488748 1.58898855 1.68263663 1.45543941 1.69865351 1.13087422
# 1.71964229 1.75642407 1.73828958 1.78190857 1.7309088  1.66905382
# 1.81711721 1.67381457 1.77568321 1.8291556 ]


tsys2 = vanecal(edge2, 329, 330, feeds=range(16))
print(tsys2)
# @114.0 GHz

#[1.10748701 1.25370611 1.27460674 1.5011671  1.39648119 1.14553602
# 1.28962367 1.37215013 1.36604623 1.36293573 1.32261247 1.35472034
# 1.40061496 1.25941423 1.47254888 1.42668759]

# ratio:        1.24 +/- 0.11

sp1,sp2 = getnod(edge2, 331, 332, 4, 7)

sp3 = 0.5*(sp1+sp2)
sp3.plot()
sp3.stats(roll=1, qac=True)   # 4.923401847859579e-05 0.0016115564993204495 -0.0007704968985869036 0.051103618263309156'


sp1,sp2 = getnod(edge2, 333, 334, 4, 7)

sp3 = 0.5*(sp1+sp2)
sp3.plot()
sp3.stats(roll=1,qac=True)   # 3.4933664104756355e-06 0.00023323680056886874 -0.0008156976082961529 0.0034146710055591675'

# no obvious signal


#%% NOD EXAMPLE-3 tp_nocal   NOD_BEAMS 0,1

_help = """
    SCAN OBJECT VELOCITY    PROC  PROCSEQN RESTFREQ DOPFREQ # IF # POL # INT # FEED     AZIMUTH   ELEVATIO
0    130    M82      0.0  CALSEQ         1   87.645    86.4    4     2    41      2  334.375845  46.552969
1    131    M82      0.0     Nod         1   87.645    86.4    4     2    61      2  334.340511  46.455481
2    132    M82      0.0     Nod         2   87.645    86.4    4     2    61      2  334.420688  46.356491
3    133    M82      0.0     Nod         1   87.645    86.4    4     2    61      2  334.270943  46.255965
4    134    M82      0.0     Nod         2   87.645    86.4    4     2    61      2  334.352313  46.156723
5    135    M82      0.0     Nod         1   87.645    86.4    4     2    61      2  334.205169   46.05736
6    136    M82      0.0     Nod         2   87.645    86.4    4     2    61      2  334.287644  45.957891
7    137    M82      0.0     Nod         1   87.645    86.4    4     2    61      2  334.142569  45.858239
8    138    M82      0.0     Nod         2   87.645    86.4    4     2    61      2  334.225833  45.757117
9    139    M82      0.0     Nod         1   87.645    86.4    4     2    61      2  334.082225  45.655873
10   140    M82      0.0     Nod         2   87.645    86.4    4     2    61      2  334.167154  45.555955
11   141    M82      0.0  CALSEQ         1   87.645    86.4    4     2    41      2  334.027735  45.461716
"""

f3 = dysh_data(accept='AGBT15B_244_07/AGBT15B_244_07.raw.vegas')
sdf3=GBTFITSLoad(f3, skipflags=True)
# 8 fits files,   2 for beams, 4 for IF  - 12 scans (2 CALSEQ - W-band receiver at 87 GHz)
sdf3.summary()
# 11072 rows

mkdir("nod3cal")
sdf3.write("nod3cal/file.fits", scan=130, ifnum=0, plnum=0, overwrite=True)
nod3cal = GBTFITSLoad("nod3cal")
nod3cal.summary()
nod3cal._index[kw]   # 82 rows  (Observing, Cold1, Cold2 and a few Unknown in state transition)

calseq(nod3cal, 130)
# calseq(sdf3, 130)
#   -> RecursionError: maximum recursion depth exceeded


v1 = nod3cal.gettp(scan=130, fdnum=0, calibrate=True, cal=False).timeaverage()
s1 = nod3cal.gettp(scan=130, fdnum=1, calibrate=True, cal=False).timeaverage()
v1.plot()
s1.plot()
t1 = s1/(v1-s1) 

tsys = vanecal(nod3cal, 130, 130, feeds=[0,1])

mkdir("nod3")
sdf3.write('nod3/file.fits', scan=[131,132], ifnum=0, plnum=0, overwrite=True)  #244

 
nod3 = GBTFITSLoad('nod3')
nod3.summary()
nod3._index[k]    # 244 rows
getbeam(nod3)     # [0,1]

sp0 = nod3.getnod(scan=[131,132],ifnum=0,plnum=0)

# .timeaverage()


nod3.gettp(scan=131,fdnum=0,calibrate=True, cal=False).timeaverage().plot()
nod3.gettp(scan=132,fdnum=0,calibrate=True, cal=False).timeaverage().plot()
nod3.gettp(scan=131,fdnum=1,calibrate=True, cal=False).timeaverage().plot()
# ok

#%% VANE4    issue 257

_help = """
   SCAN    OBJECT VELOCITY   PROC  PROCSEQN RESTFREQ DOPFREQ # IF # POL # INT # FEED AZIMUTH ELEVATIO
0     1  B0329+54      0.0  Track         1      1.5     1.5    1     2    16      1     0.0      0.0
1     2  B0329+54      0.0  Track         1      1.5     1.5    1     2    16      1     0.0      0.0
2     3  B0329+54      0.0  Track         1      1.5     1.5    1     2    16      1     0.0      0.0
3     4  B0329+54      0.0  Track         1      1.5     1.5    1     2    16      1     0.0      0.0
4     5  B0329+54      0.0  Track         1      1.5     1.5    1     2    16      1     0.0      0.0
"""

f4 = dysh_data("TSCI_RYAN_16/TSCI_RYAN_16.raw.vegas")
sdf4 = GBTFITSLoad(f4)
sdf4.summary()

kw=['DATE-OBS','SCAN', 'SUBOBSMODE', 'IFNUM', 'PLNUM', 'FDNUM', 'INTNUM', 'TCOLD', 'CALPOSITION']
sdf4._index[kw]   # 192 rows

sdf4.getps()
# Exception: Multiple SUBOBSMODE present, cannot deal with this yet ['TPWCAL', 'TPNOCAL']

# scan 1 is the TPWCAL
sp4a = sdf4.gettp(scan=1, plnum=0)
sp4b = sdf4.gettp(scan=1, plnum=1)               
    
sp4a.timeaverage().plot(ymin=-1e8, ymax=1e9, title="TSCI_RYAN_16 scan=1 plnum=0", xaxis_unit="chan")
sp4b.timeaverage().plot(ymin=-1e8, ymax=1e9, title="TSCI_RYAN_16 scan=1 plnum=1", xaxis_unit="chan")
     
sdf4.gettp(scan=2, plnum=0, calibrate=True, cal=False).timeaverage().plot(ymin=-1e8, ymax=1e9)
# ok

#%%  VANE5   argus TGBT22A_605_05 (26 GB big dataset)

_help = """
    SCAN OBJECT VELOCITY       PROC  PROCSEQN   RESTFREQ    DOPFREQ # IF # POL # INT # FEED    AZIMUTH   ELEVATIO
0     10   VANE      0.0      Track         1  93.173704  93.173704    1     1    10     16  67.295835  76.953031
1     11    SKY      0.0      Track         1  93.173704  93.173704    1     1    10     16   67.29655  76.953058
2     12   DR21      0.0      Track         1  93.173704  93.173704    1     1    60     16  66.498004  77.259536
3     13   DR21      0.0      Track         1  93.173704  93.173704    1     1    60     16  66.270661  77.4695
...
24    34   DR21      0.0      Track         1  93.173704  93.173704    1     1    60     16  57.709781  81.763408
25    35   VANE      0.0      Track         1  93.173704  93.173704    1     1    25     16    57.4161  82.059614
26    36    SKY      0.0      Track         1  93.173704  93.173704    1     1    25     16  57.417741  82.0596
27    37   DR21      0.0  RALongMap         1  93.173704  93.173704    1     1    75     16  55.839773  82.313428
...
50    60   DR21      0.0  RALongMap        24  93.173704  93.173704    1     1    75     16  39.876154  84.616655
51    61   VANE      0.0      Track         1  93.173704  93.173704    1     1    25     16  40.375997  84.666108
52    62    SKY      0.0      Track         1  93.173704  93.173704    1     1    25     16  40.377915   84.66616
"""

f5 = dysh_data(example="mapping-Argus/data/TGBT22A_603_05.raw.vegas")
sdf5 = GBTFITSLoad(f5, skipflags=True)
sdf5.summary()   # 53 scans (10..62)
# lots memory

mkdir("vane5")
sdf5.write("vane5/file.fits", scan=[10,11,12], overwrite=True)
# takes long time

vane5 = GBTFITSLoad("vane5", skipflags=True)
vane5.summary()


tsys_10 = vanecal(vane5, 10, 11, feeds=range(16))
#      1.45318402, 1.25949678, 1.27948539, 1.41625876, 1.26606387,
#      1.42277304, 1.2289571 , 1.29762214, 1.27845341, 1.24670253,
#      1.24174926, 1.32368279, 1.2663306 , 1.2432676 , 1.29513699,
#      1.23222654])
 
tsys_34 = vanecal(sdf5, 34, 35, feeds=range(16))
#      [1.46858929, 1.27776463, 1.29033069, 1.43358742, 1.28317839,
#       1.42513274, 1.24034782, 1.31302693, 1.29617406, 1.26008328,
#       1.25831389, 1.34137986, 1.06626227, 1.25749316, 1.32069702,
#       1.252414  ])

# mean and std of diff:    -0.0020865396509081036  0.05241421768622772


#%%  VANE6   W-band rxco-W/data/TSCAL_220105_W.raw.vegas  (10 MB)

_help = """
    SCAN     OBJECT VELOCITY    PROC  PROCSEQN RESTFREQ DOPFREQ # IF # POL # INT # FEED     AZIMUTH   ELEVATIO
0     32  2253+1608      0.0  CALSEQ         1     86.0    86.0    1     2    46      2  206.811054  65.771037
1     33  2253+1608      0.0     Nod         1     86.0    86.0    1     2    31      2  207.340151  65.684269
2     34  2253+1608      0.0     Nod         2     86.0    86.0    1     2    31      2  208.036343  65.599144
3     23  2253+1608      0.0  CALSEQ         1     77.0    77.0    1     2    46      2   201.75799  66.502127
4     24  2253+1608      0.0     Nod         1     77.0    77.0    1     2    31      2  202.306842  66.431865
5     25  2253+1608      0.0     Nod         2     77.0    77.0    1     2    31      2  202.895056  66.380331
6     26  2253+1608      0.0  CALSEQ         1     72.0    72.0    1     2    46      2  203.469145  66.274739
7     27  2253+1608      0.0     Nod         1     72.0    72.0    1     2    31      2  203.999448  66.200655
8     28  2253+1608      0.0     Nod         2     72.0    72.0    1     2    31      2  204.579033  66.145504
9     29  2253+1608      0.0  CALSEQ         1     82.0    82.0    1     2    46      2  205.141539  66.032763
10    30  2253+1608      0.0     Nod         1     82.0    82.0    1     2    31      2  205.671661  65.952387
11    31  2253+1608      0.0     Nod         2     82.0    82.0    1     2    31      2  206.242727  65.893704
"""

f6 = dysh_data(example="rxco-W/data/TSCAL_220105_W.raw.vegas")
sdf6 = GBTFITSLoad(f6, skipflags=True)
sdf6.summary()
sdf6._index[kw]    # 1728

mkdir("vane6")
scans = [32,33,34]
scans = [32]
scans = [32,23,26,29]
scans = [23,26,29,32]
sdf6.write("vane6/file.fits", plnum=0, scan=scans, overwrite=True)
#sdf6.write("vane6/file.fits", scan=scans, overwrite=True)

vane6 = GBTFITSLoad("vane6")
vane6.summary()
vane6._index[kw]    # 92 in just scan 32      368 if all 4 CALSEQ

calseq(vane6, 23, plnum=0, ifnum=0, fdnum=1)
# syntaxError: too many nested parentheses    <--- only if all plnum's selected to into vane6
calseq(vane6, 26, fdnum=1)
calseq(vane6, 29, fdnum=1)
calseq(vane6, 32, fdnum=1)
#    fdnum=0    608K    (612 when adding the Nod)
#          1    179K    (189 when adding the Nod)
#    something odd:   when doing [32,33,34] i get different answers from using just [32]

# getbeam fails if only scan 32 is used
getbeam(vane6)

vane6.gettp(scan=32,fdnum=0,calibrate=True, cal=False).timeaverage().plot()
vane6.gettp(scan=32,fdnum=1,calibrate=True, cal=False).timeaverage().plot()

sp1,sp2 = getnod(vane6, 33, 34, 0, 1)

sp3 = 0.5*(sp1+sp2)
sp3.plot()
sp3.stats(roll=1, qac=True)   # 4.923401847859579e-05 0.0016115564993204495 -0.0007704968985869036 0.051103618263309156'



#%%

