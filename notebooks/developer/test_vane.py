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
edge2 - VANE/SKY   at 112,114 GHz                                               some weak signal!
nod3  - CALSEQ         fdnum=[0,1]          tp nocal;  TBD
vane4 - "ISSUE 257"
vane5 - VANE/SKY
vane6 - CALSEQ    at 4 freqs

Related issues:

https://github.com/GreenBankObservatory/dysh/issues/257  
https://github.com/GreenBankObservatory/dysh/issues/457     flagging problem
https://github.com/GreenBankObservatory/dysh/issues/484     tp_nocal

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
# new
from dysh.fits.core import mean_data
from dysh.fits.core import getbeam
from dysh.fits.core import calseq
from dysh.fits.core import vanecal
from dysh.fits.core import plot_vegas
 
#  useful keys for a mult-beam observation listing

k=['DATE-OBS','SCAN', 'IFNUM', 'PLNUM', 'FDNUM', 'INTNUM', 'PROCSCAN','FEED', 'SRFEED', 'FEEDXOFF', 'FEEDEOFF', 'SIG', 'CAL', 'PROCSEQN', 'PROCSIZE']
ks=['DATE-OBS','SCAN', 'SUBOBSMODE', 'IFNUM', 'PLNUM', 'FDNUM', 'INTNUM', 'CAL', 'PROCSEQN']
kw=['DATE-OBS','SCAN', 'SUBOBSMODE', 'IFNUM', 'PLNUM', 'FDNUM', 'INTNUM', 'CAL', 'PROCSEQN', 'CALPOSITION','TCAL','TSYS']

#  some more liberal panda dataframe display options
import pandas as pd
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 1000)

# pd.options.display.max_columns = None

#%%  debugging

# ?? can one only do this once per kernel ??

import dysh
dysh.log.init_logging(3)   # 0=ERROR 1=WARNING 2=INFO 3=DEBUG
dysh.log.init_logging(1)    # 0 didn't work

#%%  helper classes

class Tsys(object):
    """ 
    Toolkit to compute the TSYS 
    
    Parameters
    ----------
    ifnum: int
        list (or single value) of the IFs
    fdnum: int
        list (or single value) of the feeds
    plnum: int
        list (or single value) of the polarizations

    In situations where computing the Tsys is cumbersome, a special Tsys class
    is available to contruct a container that can be passed to any of the get*
    routines to apply a Tsys calibration.
    
    1. Setting up a Tsys container needs the number of IF/FD/PL:
        
        mytsys = Tsys(fdnum=1, ifnum=[0,1], plnum=1)
        mydata.data[0,0,0] = 200    # ifnum=0
        mydata.data[0,1,0] = 210    # ifnum=1
        getps(ifnum=[0,1], plnum=1, fdnum=1, cal_tsys = mytsys)
  
    """
    def __init__(self, fdnum=0, ifnum=0, plnum=0, **kawrgs):
        # set the feeds
        if type(fdnum) == list:   
            self.fdnum = fdnum
        else:
            self.fdnum = [fdnum]
            
        # set the IF's
        if type(ifnum) == list:
            self.ifnum = ifnum
        else:
            self.ifnum = [ifnum]

        # set the polarizations
        if type(plnum) == list:   
            self.plnum = plnum
        else:
            self.plnum = [plnum]
        print(f"TSYS for fdnum={self.fdnum} ifnum={self.ifnum} plnum={self.plnum}")
        
        self.data = np.ones((len(self.fdnum),len(self.ifnum),len(self.plnum)))
        
    def get(self, fdnum, ifnum, plnum):
        """ get a single Tsys
        """
        # @todo check index
        return self.data[self.fdnum.index(fdnum), self.ifnum.index(ifnum), self.fdnum.index(fdnum)]
    
    def set(self, tsys, fdnum, ifnum, plnum):
       """ set a single Tsys
       """
       # @todo check index
       self.data[self.fdnum.index(fdnum), self.ifnum.index(ifnum), self.fdnum.index(fdnum)]= tsys 
        
junk = Tsys([0,1],[0,1,2,3],[0,1])
print(junk.data)

mytsys = Tsys(fdnum=1, ifnum=[0,1], plnum=1)
mytsys.data[0,0,0] = 200    # ifnum=0
mytsys.data[0,1,0] = 210    # ifnum=1
print(mytsys.data)
mytsys.set(220,1,1,1)
print(mytsys.get(1,1,1))


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

def getbeam1(sdf, debug=False):
    """ find the two nodding beams based on FDNUM, FEED 
        needs PROCSCAN='BEAM1' or 'BEAM2'
    """
    kb=['DATE-OBS','SCAN', 'IFNUM', 'PLNUM', 'FDNUM', 'PROCSCAN', 'FEED', 'SRFEED', 'FEEDXOFF', 'FEEDEOFF']
    kb=['FEEDXOFF','FEEDEOFF','PROCSCAN','FDNUM','FEED']
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
    
def mean_data1(data, fedge=0.1, nedge=None):
    """ special mean to exclude the edges like mean_tsys()
    """
    nchan = len(data)
    if nedge is None:
        nedge = int(nchan * fedge)
    chrng = slice(nedge, -(nedge - 1), 1)
    meandata = np.nanmean(data[chrng])
    return meandata
    
def vanecal1(sdf, vane_sky,  feeds=[], tcal=None, tmpfile=False):
    """ loop over feeds to get tsys factor
        for efficiency of large data, it's better to sdf.write() and use only the
        vane and sky scans'
        Example EDGE:  (1.3GB)
             all 163 scans:     6m33s    3m6s
             write small sdf:     25s    
             just vane/cal:       20s
    """
    vane = vane_sky[0]
    sky = vane_sky[1]
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
        
    if True:
        # until we figure this out via getatmos
        tcal = 100.0    # set to 100K for now

    
    i = 0
    for f in feeds:
        v = sdf1.gettp(scan=vane, fdnum=f, calibrate=True, cal=False).timeaverage()
        s = sdf1.gettp(scan=sky,  fdnum=f, calibrate=True, cal=False).timeaverage()
        mean_off = mean_data(s.data)
        mean_dif = mean_data(v.data - s.data)
        tsys[i] = tcal * mean_off/mean_dif
        #  vanecal.pro seems to do    tcal / median( (v-s)/s)
        i = i + 1
    print("Warning, assumed TCAL=",tcal)
    # print("TSYS=",tsys)
     
    return tsys

    
def calseq1(sdf, scan, tcold=54, fdnum=0, ifnum=0, plnum=0, freq=None):
    """ W-band receivers use a CALSEQ
        This routine returns the gain and Tsys for W-band channel
        
        Tcold = 54 - 0.6*(FREQ-77)      FREQ in GHz
    """
    if freq is not None:
        # see eq.(13) in GBT memo 302
        tcold = 54 - 0.6*(freq-77)
        print(f"Warning: calseq using freq={freq} GHz and setting tcold={tcold} K")
        
    twarm = sdf._index['TAMBIENT'].mean()
        

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

def getnod(sdf, scans, beams, ifnum=0, plnum=0, tsys=None):
    """ fake getnod() based on alternating gettp() with averaging done internally
        use the real sdf.getnod() for final analysis
        new API:    scans[], feeds[], tsys[]
        """
    if tsys is None:
        tsys = [1.0, 1.0]
    else:
        print(f"@todo use tsys={tsys} for proper weighting")
    ps1_on = sdf.gettp(scan=scans[0], fdnum=beams[0], ifnum=ifnum, plnum=plnum, calibrate=True, cal=False).timeaverage()
    ps1_off = sdf.gettp(scan=scans[1], fdnum=beams[0], ifnum=ifnum, plnum=plnum, calibrate=True, cal=False).timeaverage()
    sp1 = (ps1_on - ps1_off)/ps1_off * tsys[0]

    ps2_on = sdf.gettp(scan=scans[1], fdnum=beams[1], ifnum=ifnum, plnum=plnum, calibrate=True, cal=False).timeaverage()
    ps2_off = sdf.gettp(scan=scans[0], fdnum=beams[1], ifnum=ifnum, plnum=plnum, calibrate=True, cal=False).timeaverage()
    sp2 = (ps2_on - ps2_off)/ps2_off * tsys[1]

    return (sp1,sp2)


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


def plot_vegas1(sdf, scans, title = None, tsys = False, inverse=False, edge=50, ylim=None):
    """   plot vegas 16 beams
    """
    fig,ax = plt.subplots(4,4, sharex='col', sharey='row',
                          gridspec_kw={'hspace': 0, 'wspace': 0})
    
    for r in range(4):
        for c in range(4):
            p = r*4 + c
            ax[r,c].plot()
            #ax[r,c].set_xlabel(f"{r} {c} {p}")
            for s in scans:
                v1 = sdf.gettp(scan=s, fdnum=p, calibrate=True, cal=False).timeaverage()
                if tsys:
                    s1 = sdf.gettp(scan=s+1, fdnum=p, calibrate=True, cal=False).timeaverage()        
                    vs = s1.data/(v1.data - s1.data)
                    if inverse:
                        vs = 1/vs
                else:
                    vs = v1.data
                ax[r,c].plot(vs[edge:-edge], label=f"{p}")
                if tsys:
                    # note wwe haven't put a proper channel / freq axis
                    nc = len(vs[edge:-edge])
                    fix = np.ones(nc)
                    ax[r,c].plot(fix, color='black')
                #ax[r,c].scatter([0,900],[vs[edge],vs[-edge]], label=f"{p}")
                
                if ylim is not None:
                    ax[r,c].set_ylim(ylim)
            ax[r,c].legend()
    if title is None:
        plt.suptitle(sdf.filenames()[0])
    else:
        plt.suptitle(title)
    plt.tight_layout()
    
    if tsys:
        print(f"Showing sky/(vane-sky) for scans={scans},scans+1")
    else:
        print(f"Showing total power for scans={scans}")
              


#%%  classic tp/ps

f1 = dysh_data(test="getps")        # OnOff, small one (1 IF, 1 PL, 1 INT)
#f1 = dysh_data(example="getps")    # OnOff, long one (4 IFs, 2 PLs, 151 INTs)
#f1 = dysh_data(test="AGBT18B_354_03/AGBT18B_354_03.raw.vegas") # OffOn
fn1 = f1.parts[-1]
print(f"Using {fn1}")


sdf0 = GBTFITSLoad(f1)
sdf0.summary()
sdf0.summary(verbose=True)
sdf0._index[k]

sp1 = sdf0.getps().timeaverage()
sp1.plot(title='getps: %s' % fn1)
sp1.stats(qac=True)
# 0.21853874407385052 0.6796743343920357 -3.7057502269744873 4.343878746032715

tp1a = sdf0.gettp(scan=152).timeaverage()
tp1b = sdf0.gettp(scan=153).timeaverage()

# why can't I do   len(tp1a), but I need len(tp1a.data)

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


f1 = dysh_data(accept='AGBT22A_325_15/AGBT22A_325_15.raw.vegas')                  # accept='nod1'
sdf1=GBTFITSLoad(f1)   # skipflags=True)
# 8 files, 16 beams, each file has 2 beams - 4 scans, VANE/SKY/Nod/Nod
sdf1.summary()    # 256 rows
sdf1.flags.show() # 64 rows

# extract 290 and 289 (note order is odd in sdfits:   290 came before 289 (odd))
mkdir("vane1")
scans=[281,282]
sdf1.write('vane1/file.fits',scan=scans, overwrite=True)   # 64 rows


# our EDGE notes claim that 1,9 are the nodding beams
feeds1 = getbeam(sdf1)    # [10,1]
print("Nodding feeds",feeds1)

# load smaller dataset
vane1 = GBTFITSLoad('vane1')
vane1.summary()
vane1._index[kw]    # 64 rows  -- TPNOCAL
getbeam(vane1)   # -> notice not defined here in just the vane/sky

#  in these, feed=2 is bad
plot_vegas(vane1,[281,282])
#plot_vegas(vane1,[281],tsys=True,ylim=[0.45,0.7])
plot_vegas(vane1,[281],tsys=True)

#  do a single vane/sky comparison for fdnum=8 (the first entry, but not a nodding beam)
v1 = vane1.getspec(0).flux.value      #  vane
s1 = vane1.getspec(4).flux.value      #  sky
t1 = s1/(v1-s1)    
t1 = t1[100:900]                      # 0.53 +/0 0.02
np.nanmean(t1)  # 0.52613854
np.nanstd(t1)   # 0.016047847
plt.plot(t1)

# mean_tsys adds tcal/2
tsys1 = mean_tsys(v1, s1, 1, mode=1) - 0.5
tsys0 = mean_tsys(v1, s1, 1, mode=0) - 0.5
print("mean_tsys mode0,1=",tsys0,tsys1)

t2 = (v1-s1)/s1
t2 = t2[100:900]                      
np.nanmean(t2) 
np.nanstd(t2)  
plt.plot(t2)

#  @todo request:   what's the variation in time? why is this so hard
v1 = vane1.gettp(scan=281, fdnum=8, calibrate=True, cal=False).timeaverage()
s1 = vane1.gettp(scan=282, fdnum=8, calibrate=True, cal=False).timeaverage()

# nchan is 1024
t1 = s1/(v1-s1) 
# t1 = t1[100:900]
# ISSUE: UnitTypeError: Longitude instances require units equivalent to 'rad', so cannot set it to ''
t1 = t1.data[100:900]
plt.plot(t1)

v1.plot(title='VANE')
s1.plot(title='SKY')

tsys=vanecal(vane1, [281, 282], feeds=feeds1)
print(f"For feeds={feeds1} Tsys={tsys}")       #  [50.41218948 57.76417121]

# now process the NOD assuming the just aqcuired Tsys
sp1,sp2 = getnod(sdf1, [289, 290], feeds1, tsys)
sp1.plot(title='Feed1')
sp2.plot(title='Feed2')

sp3 = 0.5*(sp1+sp2)
sp3 *= tsys.mean()
object = sp3.meta['OBJECT']
sp3.plot(title='Averaged Feeds')
sp3.plot(title=f"Source: {object}",xaxis_unit="km/s")               # no obvious signal - still inherits the title (bug??)

sp3.stats(qac=True)
sp3.stats(qac=True, roll=1)



sdf1.getnod(scan=[289,290],fdnum=[10,1])
# ISSUE#484: IndexError: index 0 is out of bounds for axis 0 with size 0


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

sdf1 = GBTFITSLoad(dysh_data('AGBT21B_024_01/AGBT21B_024_01.raw.vegas'), skipflags=True)
#  ISSUE: CPU=9s with skipflags,  but 9min 5 sec with skipflags
a = sdf1.summary()  # 208 scans   1363MB
print(a)
b=a[a["OBJECT"] == "VANE"]
print(b)   # there are 13 vane's in here (17,21,58,65,102,111,115,152,196,198,212,249)

# comparing the SKY tp at two different times
plot_vegas(sdf1,[17,21])    # this is SLOOOW, see using it frome edge1 below

# make a smaller dataset for testing with just the sky/vane's
mkdir("edge1")
scans = [17,18,21,22]
# sdf1.write("edge1/file.fits", fdnum=8, scan=scans, intnum=[0,1], overwrite=True)
# sdf1.write("edge1/file.fits", fdnum=8, scan=scans, overwrite=True)
# sdf1.write("edge1/file.fits", scan=scans, overwrite=True)                # 4576 rows  -  also fails if skipflags=True
#   ISSUE:    TypeError: object of type 'Column' has no len()
#  CPU times: user 23.5 s, sys: 3.02 s, total: 26.5 s   Wall time: 25.2 s
sdf1.write("edge1/file.fits", scan=scans, overwrite=True, flags=False)     # now ok
#    temp fix

edge1 = GBTFITSLoad("edge1")
edge1.summary()
edge1._index[ks]

# in these data feed=2 is ok
# but feed=8 changed a lot from 17 to 21
plot_vegas(edge1,[17,21])      # super fast now
# looking at the spectral Tsys, feed 8 didn't change like their TP did
plot_vegas(edge1,[17,21],tsys=True)

beams1=getbeam(sdf1)   # 4,7
print(beams1)

tsys = vanecal(sdf1, [17, 18], feeds=beams1)    # 12 sec
print(tsys)
# [61.12006835 63.01095552]


tsys = vanecal(edge1, [17, 18], feeds=beams1)   #  0.8 sec
print(tsys)
# [61.12006835 63.01095552].


# plotting a TP spectrum
sdf1.gettp(scan=17, fdnum=8, calibrate=True, cal=False).timeaverage().plot()
edge1.gettp(scan=17, fdnum=8, calibrate=True, cal=False).timeaverage().plot()
# ok !!

# not working yet, since it wants ON/OFF, needs the GETTP hack
sp2 = sdf1.getnod(scan=19)
# IndexError: index 0 is out of bounds for axis 0 with size 0


# if getnod() is not converted, we can use the old trick of doing two gettp()
scans = [19,20]
sp1,sp2 = getnod(sdf1, scans, beams1)

sp3 = 0.5*(sp1+sp2)
sp3 *= tsys.mean()
sp3 = sp3[100:900]
object = sp3.meta['OBJECT']
sp3.plot(title=f"Source: {object}",xaxis_unit="km/s")     # no obvious signal, pretty bad baseline
# center is at 0 km/s, where VLSR 4500 km/s at 113.57 is

sp3.stats(qac=True)
sp3.stats(qac=True, roll=1)

#%% EDGE2 nod example on NGC5908       1640 MB

_help = """
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


sdf2=GBTFITSLoad(dysh_data('AGBT21B_024_14/AGBT21B_024_14.raw.vegas'), skipflags=True)
sdf2.summary()

mkdir("edge2")
scans = [327,328,329,330,331,332,333,334]
sdf2.write("edge2/file.fits", scan=scans, overwrite=True)                # 4576 rows
#  CPU times: user 23.5 s, sys: 3.02 s, total: 26.5 s   Wall time: 25.2 s

beam2 = getbeam(sdf2)   # 1,9
print("feeds",beam2)

edge2 = GBTFITSLoad("edge2")
edge2.summary()
edge2._index[kw] # 5248 rows

plot_vegas(edge2,[327,329],"edge2 TP at 112 and 114 GHz")
# this shows 327 is generally better tsys, except feed=5 about the same (but high)
plot_vegas(edge2,[327,329],"edge2 Tsys at 112 and 114 GHz",tsys=True, ylim=[0.5,1.1])

# this was at 112 GHz for NGC570
tsys1 = vanecal(edge2, [327, 328], feeds=beam2)
print(tsys1)
# @112.3 GHz   [62.93311553 56.119602  ]
# [0.70676998 0.62933116 0.59430538 0.68707773 0.58870158 0.88427164
#  0.58151629 0.56933859 0.57527814 0.56119602 0.57773119 0.59914186
#  0.55032223 0.59743774 0.56316352 0.54670035]

# this was at 114 GHz for NGC5908
tsys2 = vanecal(edge2, [329, 330], feeds=beam2)
print(tsys2)
# @114.0 GHz[   79.76351041 73.37103122]
# [0.90294513 0.7976351  0.78455571 0.66614836 0.71608555 0.87295378
#  0.77542001 0.72878323 0.73203965 0.73371031 0.75607937 0.73815973
#  0.7139721  0.79401993 0.67909461 0.7009243 ]


# ratio:        1.24 +/- 0.11

# sp1,sp2 = getnod(edge2, 331, 332, 4, 7)
sp1,sp2 = getnod(edge2, [331, 332], beam2)

#sp3a= 0.5*(sp1+sp2)
sp3a = sp1.average(sp2)
sp3a *= tsys1.mean()     # @todo how to properly weigh these spectra    weight = delta-T * delta-F / Tsys
object = sp3a.meta['OBJECT']
sp3a.plot(title=f"Source: {object}",xaxis_unit="km/s")   # no obvious signal - huge amp wave 0.5
sp3a.stats(qac=True)
sp3a.stats(roll=1, qac=True)   # -0.0007032170649519567 0.02473309930926793 -0.6811288707087332 0.038270612438221545


sp1,sp2 = getnod(edge2, [333, 334], beam2)

#sp3b = 0.5*(sp1+sp2)
sp3b = sp1.average(sp2)
sp3b *= tsys2.mean()
object = sp3b.meta['OBJECT']
sp3b.plot(title=f"Source: {object}",xaxis_unit="km/s")  # some signal !!!   amp wave now 0.1
sp3b.stats(qac=True)
sp3b.stats(roll=1,qac=True)   # 1.159984571846731e-05 0.016863803397811985 -0.05766586065002078 0.04746087356417471'


if False:
    sp1,sp2 = getnod(edge2, [331, 332], beam2, tsys=tsys1)
    sp3,sp4 = getnod(edge2, [333, 334], beam2, tsys=tsys2)
    sp3 = sp1.average([sp2,sp3,sp4])
else:    
    sp3 = (sp3a + sp3b)/2
    sp3 = sp3a.average(sp3b)
sp3.plot(title=f"Source: {object}",xaxis_unit="km/s")

# do a baseline subtraction
kms = u.km/u.s
sp4=sp3[400:600]    # @todo need to figure this out in km/s
sp4.baseline(model="poly", degree=5, exclude=[-150*kms,150*kms], remove=True)
#sp4.baseline(model="cheby", degree=5, exclude=[-150*kms,150*kms], remove=True)
# chebyshev', 'legendre', or 'hermite'
print(sp4.baseline_model)
sp4.plot(title=f"Source: {object}",xaxis_unit="km/s")

sp5 = sp4.smooth('box', 3)
sp5.plot(title=f"Source: {object}",xaxis_unit="km/s")

#%% NOD EXAMPLE-3 tp_nocal   NOD_BEAMS 0,1       729 MB

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

# note:   VHEL = 269 km/s   VLSR=275   -   but why is dopfreq so much lower ???

f3 = dysh_data(accept='AGBT15B_244_07/AGBT15B_244_07.raw.vegas')
sdf3=GBTFITSLoad(f3, skipflags=True)
# 8 fits files,   2 for beams, 4 for IF  - 12 scans (2 CALSEQ - W-band receiver at 87 GHz)
sdf3.summary()
# 11072 rows

tsys3,g = calseq(sdf3, 130) 
print(tsys3,g)                  # 100.13203834626455 9.115908926574802e-07

mkdir("nod3cal")
sdf3.write("nod3cal/file.fits", scan=130, ifnum=0, plnum=0, overwrite=True)


nod3cal = GBTFITSLoad("nod3cal")
nod3cal.summary()
nod3cal._index[kw]   # 82 rows  (Observing, Cold1, Cold2 and a few Unknown in state transition)  -- TPNOCAL

tsys3,g = calseq(nod3cal, 130)     # 100.27992034259859, 9.129371938341321e-07
print(tsys3,g)
# calseq(sdf3, 130)      # 100.13203834626455, 9.115908926574802e-07   - if/pol multivalued....
#   -> RecursionError: maximum recursion depth exceeded
#   now ok


#   131, 133, 135, 137, 139
s = 131
s = 133
s = 135
s = 137
s = 139
mkdir("nod3")
sdf3.write('nod3/file.fits', scan=[s,s+1], ifnum=0, plnum=0, overwrite=True)  #244

nod3 = GBTFITSLoad('nod3')
nod3.summary()
nod3._index[k]    # 244 rows

beams3 = getbeam(nod3)     # [0,1]
print(beams3)

sp1,sp2 = getnod(nod3, [s,s+1],beams3)
sp3 = 0.5*(sp1+sp2)
sp3 *= tsys

object = sp3.meta['OBJECT']
sp3.plot()
sp3.plot(title=f"Source: {object}",xaxis_unit="km/s")

sp3.stats(qac=True)
sp3.stats(roll=1,qac=True)   

sp3.smooth('box',32).plot(title=f"Source: {object}",xaxis_unit="km/s")


sp0 = nod3.getnod(scan=[s,s+1],ifnum=0,plnum=0)
#   IndexError: index 0 is out of bounds for axis 0 with size 0
# .timeaverage()

nod3.gettp(scan=131,fdnum=0,calibrate=True, cal=False).timeaverage().plot()
nod3.gettp(scan=132,fdnum=0,calibrate=True, cal=False).timeaverage().plot()
nod3.gettp(scan=131,fdnum=1,calibrate=True, cal=False).timeaverage().plot()
# ok

sp = []
for s in range(131,140,2):
    print(s)
    sp1,sp2 = getnod(sdf3, [s,s+1], [0,1], plnum=0, ifnum=1)
    sp3,sp4 = getnod(sdf3, [s,s+1], [0,1], plnum=1, ifnum=1)
    sp.append(sp1)
    sp.append(sp2)
    sp.append(sp3)
    sp.append(sp4)
    
sp_final = sp[0].average(sp[1:]).smooth('gauss',16)
sp_final.plot()
    
    

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
sdf4._index[kw]   # 192 rows   TCAL ~ 1.5   TPNOCAL

sdf4.getps()
# Exception: Multiple SUBOBSMODE present, cannot deal with this yet ['TPWCAL', 'TPNOCAL']

# scan 1 is the TPWCAL
sp4a = sdf4.gettp(scan=1, plnum=0)
sp4b = sdf4.gettp(scan=1, plnum=1)               
    
sp4a.timeaverage().plot(ymin=-1e8, ymax=1e9, title="TSCI_RYAN_16 scan=1 plnum=0", xaxis_unit="chan")
sp4b.timeaverage().plot(ymin=-1e8, ymax=1e9, title="TSCI_RYAN_16 scan=1 plnum=1", xaxis_unit="chan")
     
sdf4.gettp(scan=2, plnum=0, calibrate=True, cal=False).timeaverage().plot(ymin=-1e8, ymax=1e9)
# ok

#%%  VANE5   argus TGBT22A_605_05 (26 GB big dataset)   N2H+

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

# needs monitoring for performance
f5 = dysh_data(example="mapping-Argus/data/TGBT22A_603_05.raw.vegas")
sdf5 = GBTFITSLoad(f5, skipflags=True)
# CPU times: user 4.74 s, sys: 4.81 s, total: 9.55 s Wall time: 1min 20s     34.8g  13.7g   6.4g 
# CPU times: user 47.7 s, sys: 1.53 s, total: 49.3 s Wall time: 48.3 s
sdf5.summary()   # 53 scans (10..62)
# lots memory

mkdir("vane5")
sdf5.write("vane5/file.fits", scan=[10,11,12], overwrite=True)
# takes long time, also seemed to use 130GB on e.g. lma, so we keep a local copy as well (820MB)

vane5 = GBTFITSLoad("vane5", skipflags=True)
vane5.summary()
vane5._index[kw]      # 2560 ;    FSW12NOCAL

getbeam(vane5)   # no beams, no nodding here

tsys_10 = vanecal(vane5, [10, 11], feeds=range(16))
print(tsys_10)
# [0.95318402 0.75949678 0.77948539 0.91625876 0.76606387 0.92277304
#  0.7289571  0.79762214 0.77845341 0.74670253 0.74174926 0.82368279
#  0.7663306  0.7432676  0.79513699 0.73222654]

# tsys_10a  = vanecal(sdf5, 10, 11, feeds=range(16))

# on vane5:   CPU times: user 9.28 s, sys: 166 ms, total: 9.44 s        Wall time: 9.48     (on lma;  
#             CPU times: user 6.33 s, sys: 36 ms, total: 6.36 s         Wall time: 6.33 
# on sdf5:    CPU times: user 1min 33s, sys: 3.73 s, total: 1min 37s    Wall time: 1min 37s


# mean and std of diff:    -0.0020865396509081036  0.05241421768622772

tsys_35 = vanecal(vane5, [35, 36], feeds=range(16))
print(tsys_35)

plot_vegas(vane5, [10,11])
plot_vegas(vane5, [10], tsys=True, edge=2000)

plot_vegas(vane5, [35,36])
plot_vegas(vane5, [35], tsys=True, edge=2000)

plot_vegas(vane5, [10,35], tsys=True, edge=2000)


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

beams6 = getbeam(sdf6)    # 0,1
print(beams6)

mkdir("vane6")
scans = [32,33,34]
scans = [32]
scans = [32,23,26,29]   # why is order scans not the same as input?  here i see  32,26,29,23
scans = [23,26,29,32]
sdf6.write("vane6/file.fits", plnum=0, scan=scans, overwrite=True)
#sdf6.write("vane6/file.fits", scan=scans, overwrite=True)

vane6 = GBTFITSLoad("vane6")
vane6.summary()
vane6._index[kw]    # 92 in just scan 32      368 if all 4 CALSEQ  --  TPNOCAL

calseq(vane6, 23, plnum=0, ifnum=0, fdnum=1)
# syntaxError: too many nested parentheses    <--- only if all plnum's selected to into vane6
calseq(vane6, 26, fdnum=1)
calseq(vane6, 29, fdnum=1)
calseq(vane6, 32, fdnum=0)
calseq(vane6, 32, fdnum=1)
#    fdnum=0    608K    (612 when adding the Nod)
#          1    179K    (189 when adding the Nod)
#    something odd:   when doing [32,33,34] i get different answers from using just [32]


vane6.gettp(scan=32,fdnum=0,calibrate=True, cal=False).timeaverage().plot()
vane6.gettp(scan=32,fdnum=1,calibrate=True, cal=False).timeaverage().plot()

sp1,sp2 = getnod(vane6, [33, 34], beams6)

sp3 = 0.5*(sp1+sp2)
sp3 = sp3[50:-50]
sp3.plot(xaxis_unit="km/s")              #   no obvious signal, very bad baselines

sp3.stats(qac=True)
sp3.stats(roll=1, qac=True)   

#%%
sdf8=GBTFITSLoad(dysh_data('AGBT21B_024_20/AGBT21B_024_20.raw.vegas'), skipflags=True)
sdf8.summary()

#%%
sdf8=GBTFITSLoad(dysh_data('AGBT21B_024_21/AGBT21B_024_21.raw.vegas'), skipflags=True)
sdf8.summary()

"""
   207             VANE      0.0   Track   1  114.030   1   21  16  346.4  72.4
   208              SKY      0.0   Track   1  114.030   1   21  16  346.4  72.4
   209          NGC5908      0.0     Nod   1  114.030   1   61  16  345.5  72.4
   210          NGC5908      0.0     Nod   2  114.030   1   61  16  345.2  72.3
   211          NGC5908      0.0     Nod   1  114.030   1   61  16  344.9  72.3
   212          NGC5908      0.0     Nod   2  114.030   1   61  16  344.6  72.2
"""

mkdir("edge4")
scans = [207,208,209,210,211,212]
sdf8.write("edge4/file.fits", scan=scans, overwrite=True)

beams8 = getbeam(sdf8)
print(beams8)                     # [1, 9]

edge4 = GBTFITSLoad("edge4")
edge4.summary()
edge4._index[ks]

# in these data feed=2 is ok
# but feed=8 changed a lot from 17 to 21
plot_vegas(edge4,[207,208])      # super fast now
# looking at the spectral Tsys, feed 8 didn't change like their TP did
plot_vegas(edge4,[207],tsys=True, ylim=[0.5,0.9])


tsys4 = vanecal(edge4, [207,208], feeds=beams8)    
print(tsys4)
# [74.92613003 68.55878443]




sp1,sp2 = getnod(edge4, [209, 210], beams8)

sp3a= 0.5*(sp1+sp2)
sp3a *= tsys4.mean()     # @todo how to properly weigh these spectra    weight = delta-T * delta-F / Tsys
object = sp3a.meta['OBJECT']
sp3a.plot(title=f"Source: {object}",xaxis_unit="km/s")   # no obvious signal - huge amp wave 0.5
sp3a.stats(qac=True)
sp3a.stats(roll=1, qac=True)   # -0.0007032170649519567 0.02473309930926793 -0.6811288707087332 0.038270612438221545



sp1,sp2 = getnod(edge4, [211, 212], beams8)

sp3b= 0.5*(sp1+sp2)
sp3b *= tsys4.mean()     # @todo how to properly weigh these spectra    weight = delta-T * delta-F / Tsys
object = sp3b.meta['OBJECT']
sp3b.plot(title=f"Source: {object}",xaxis_unit="km/s")   # no obvious signal - huge amp wave 0.5
sp3b.stats(qac=True)
sp3b.stats(roll=1, qac=True)   # -0.0007032170649519567 0.02473309930926793 -0.6811288707087332 0.038270612438221545


#%%
sdf61=GBTFITSLoad(dysh_data('AGBT21B_024_61/AGBT21B_024_61.raw.vegas'), skipflags=True)
sdf61.summary()
# scans 13..163

mkdir("edge61")
scans = [13,14]
sdf61.write("edge61/file.fits", scan=scans, overwrite=True)


edge61 = GBTFITSLoad("edge61")
edge61.summary()
edge61._index[ks]


vanecal(edge61, [13,14], range(16), tcal=451.59348)

plot_vegas(edge61,[13,14])
plot_vegas(edge61,[13], tsys=True)

