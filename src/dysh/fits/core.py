"""
Core functions for FITS/SDFITS
"""

import numpy as np
import matplotlib.pyplot as plt

def default_sdfits_columns():
    """The default column names for GBT SDFITS.

    Returns
    -------
        colnames : list
        A list of the GBT SDFITS column names
    """
    colnames = [
        "OBJECT",
        "BANDWID",
        "DATE-OBS",
        "DURATION",
        "EXPOSURE",
        "TSYS",
        "TDIM7",
        "TUNIT7",
        "CTYPE1",
        "CRVAL1",
        "CRPIX1",
        "CDELT1",
        "CTYPE2",
        "CRVAL2",
        "CTYPE3",
        "CRVAL3",
        "CRVAL4",
        "OBSERVER",
        "OBSID",
        "SCAN",
        "OBSMODE",
        "FRONTEND",
        "TCAL",
        "VELDEF",
        "VFRAME",
        "RVSYS",
        "OBSFREQ",
        "LST",
        "AZIMUTH",
        "ELEVATIO",
        "TAMBIENT",
        "PRESSURE",
        "HUMIDITY",
        "RESTFREQ",
        "DOPFREQ",
        "FREQRES",
        "EQUINOX",
        "RADESYS",
        "TRGTLONG",
        "TRGTLAT",
        "SAMPLER",
        "FEED",
        "SRFEED",
        "FEEDXOFF",
        "FEEDEOFF",
        "SUBREF_STATE",
        "SIDEBAND",
        "PROCSEQN",
        "PROCSIZE",
        "PROCSCAN",
        "PROCTYPE",
        "LASTON",
        "LASTOFF",
        "TIMESTAMP",
        "QD_XEL",
        "QD_EL",
        "QD_BAD",
        "QD_METHOD",
        "VELOCITY",
        "FOFFREF1",
        "ZEROCHAN",
        "ADCSAMPF",
        "VSPDELT",
        "VSPRVAL",
        "VSPRPIX",
        "SIG",
        "CAL",
        "CALTYPE",
        "TWARM",
        "TCOLD",
        "CALPOSITION",
        "BACKEND",
        "PROJID",
        "TELESCOP",
        "SITELONG",
        "SITELAT",
        "SITEELEV",
        "IFNUM",
        "PLNUM",
        "FDNUM",
        "INT",
        "INTNUM",  # not all SDFITS files have INT, so we always create INTNUM
        "NSAVE",
        # The following are added by the GBTFITSLoad constructor.
        # Arguable whether they should be included or not.
        "HDU",
        "BINTABLE",
        "ROW",
        "PROC",
        "OBSTYPE",
        "SUBOBSMODE",
    ]
    return colnames

# check w/ mean_tsys
def mean_data(data, fedge=0.1, nedge=None, median=False):
    """ special mean to exclude the edges like mean_tsys() which is like dcmeantsys()
    """
    nchan = len(data)
    if nedge is None:
        nedge = int(nchan * fedge)
    chrng = slice(nedge, -(nedge - 1), 1)
    if median:
        meandata = np.nanmedian(data[chrng])
    else:
        meandata = np.nanmean(data[chrng])
    return meandata

def getbeam(sdf, debug=False):
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
    
def calseq(sdf, scan, tcold=54, fdnum=0, ifnum=0, plnum=0, freq=None):
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


def vanecal(sdf, vane_sky,  feeds=range(16), mode=2, tcal=None, verbose=True):
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
    sdf1 = sdf
       
    #  for VANE/CAL data usually tcal=1
    if tcal is None:
        tcal = sdf._index['TCAL'].mean()
        if tcal == 1.0:
            # until we figure this out via getatmos
            tcal = 100.0    # set to 100K for now

    i = 0
    for f in feeds:
        v = sdf1.gettp(scan=vane, fdnum=f, calibrate=True, cal=False).timeaverage()
        s = sdf1.gettp(scan=sky,  fdnum=f, calibrate=True, cal=False).timeaverage()
        if mode == 0:
            mean_off = mean_data(s.data)
            mean_dif = mean_data(v.data - s.data)
            tsys[i] = tcal * mean_off/mean_dif
        elif mode == 1:
            tsys[i] = tcal / mean_data( (v.data-s.data)/s.data )
        elif mode == 2:
            tsys[i] = tcal / np.nanmedian( (v.data-s.data)/s.data )
        #  vanecal.pro seems to do    tcal / median( (v-s)/s)
        #  as well as not take off the edges
        i = i + 1
    if verbose:
        for i in range(len(feeds)):
            print(f"fdnum,Tsys   {feeds[i]:2d}  {tsys[i]:10.5f}")
        print(f"<Tsys>  {np.nanmean(tsys):.5f} +/- {np.nanstd(tsys):.5f}")
        print(f"mode={mode}")
        print("TCAL=",tcal)
     
    return tsys

def plot_vegas(sdf, scans, title = None, tsys = False, inverse=False, edge=50, ylim=None):
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
              
