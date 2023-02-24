#!/usr/bin/env python
# coding: utf-8

# # Revised GBT classes
# 
# update the classes from gbtoy as we redesign the package. 
# 
# 
# This should eventually reproduce Example 1 (position switching) from the GBTIDL manual. The datafile **ngc5291.fits** you need is [here](http://safe.nrao.edu/wiki/pub/GB/Data/GBTIDLExampleAndSampleData/ngc5291.fits) or locally on **/n/chara/teuben/GBT**.  We will also show regressions values with the GBTIDL spectra.

# ## Example
# 
# We start off with a session as a GBT user might see this:
# 
#         ps = GBTLoadPS('a.fits')        # load the SDFITS file
#         ps.summary()                    # overview times, sources, scans etc.
#         ps.finalspectrum()              # calibrate and time/pol/scan average
#         ps.plot()                       # review the final plot 
#         ps.save('a1.fits')              # save the spectrum (also SDFITS format)
#         
# This is an example of a well behaved spectrum. No masking, no baseline fitting, just simple averaging. 
# 
# All the way at the end of this notebook these 5 lines will be reviewed again.

import cProfile, pstats,io
from pstats import SortKey
import sys
import os
import time
import copy

from astropy.io import fits
import numpy as np
import matplotlib 
import pandas as pd
#from astropy.io.table import Table
#import mpl_animators
from matplotlib import pyplot as plt
from astropy.modeling import models, fitting
from specutils import Spectrum1D, SpectrumList
import astropy.units as u
from astropy.wcs import WCS


# # Helper functions
# 
# * **get_size**:     get the memory footprint of an object, recursively includes contained objects
# * **dcmeantsys**:   calibration routine to get Tsys from calON/calOFF noise diode
# * **uniq**:         returns unique values in the order from the list (unlike np.unique)
# * **sonoff**:       helper to find the On and Off scan numbers

# In[87]:


def get_size(obj, seen=None):
    #https://goshippo.com/blog/measure-real-size-any-python-object/
    """Recursively finds size of objects"""
    size = sys.getsizeof(obj)
    if seen is None:
        seen = set()
    obj_id = id(obj)
    if obj_id in seen:
        return 0
    # Important mark as seen *before* entering recursion to gracefully handle
    # self-referential objects
    seen.add(obj_id)
    if isinstance(obj, dict):
        size += sum([get_size(v, seen) for v in obj.values()])
        size += sum([get_size(k, seen) for k in obj.keys()])
    elif hasattr(obj, '__dict__'):
        size += get_size(obj.__dict__, seen)
    elif hasattr(obj, '__iter__') and not isinstance(obj, (str, bytes, bytearray)):
        size += sum([get_size(i, seen) for i in obj])
    return size


# In[88]:


def dcmeantsys(calon, caloff, tcal, mode=0, fedge=10, nedge=None):
    """
    following the GBTIDL routine with same name, get the tsys from 
    the neighboring calon and caloff we define an extra way to set 
    the edge size, nedge, if you prefer to use number of edge channels
    instead of the inverse fraction
    
    calon/caloff is meant to reflect the state of the noise diode
    
    mode=0     do the mean before the division
    mode=1     do the mean after the division
    """
    nchan = len(calon)
    if nedge == None:
        nedge = nchan // fedge     # 10 %
    if mode == 0:
        meanoff = np.mean(caloff[nedge:-nedge])
        meandiff = np.mean(calon[nedge:-nedge] - caloff[nedge:-nedge])
        meanTsys = ( meanoff / meandiff * tcal + tcal/2.0 )
    else:
        meanTsys = np.mean( caloff[nedge:-nedge] / (calon[nedge:-nedge] - caloff[nedge:-nedge]) )
        meanTsys = meanTsys * tcal + tcal/2.0
    return meanTsys


# In[89]:


def uniq(seq):
    """ from http://stackoverflow.com/questions/480214/how-do-you-remove-duplicates-from-a-list-in-python-whilst-preserving-order """
    seen = set()
    seen_add = seen.add
    return [ x for x in seq if x not in seen and not seen_add(x)]


# In[90]:


def sonoff(scan, procseqn):
    """
    return the list of On and Off scan numbers
    there must be a more elegant python way to do this....
    """
    sp = {}
    for (i,j) in zip(scan, procseqn):
        sp[i] = j
    
    us1 = uniq(scan)
    up1 = uniq(procseqn)
    
    sd = {}
    for i in up1:
        sd[i] = []
        
    for s in us1:
        sd[sp[s]].append(s)

    return sd


# # Spectrum (not clear this is even needed)
# 
# The class that contains one spectrum, something like Spectrum1D in specutils or Spectrum in pyspeckit. For now we are hardcoding the spectral axis in km/s, assuming the RESTFREQ is good. Doppler tracking ignored here for now.
# 

# In[91]:

class Spectrum(Spectrum1D):
    """
    Spectrum class built on top of specutils.Spectrum1D
    """
    def __init__(self, flux=None, spectral_axis=None, wcs=None,
             velocity_convention=None, rest_value=None, redshift=None,
             radial_velocity=None, bin_specification=None, **kwargs):
        super().__init__(flux,spectral_axis,wcs,velocity_convention,rest_value,redshift,radial_velocity,bin_specification,**kwargs)
        self._baseline = np.zeros(np.shape(flux))
        self._weights = np.ones(np.shape(flux))
        self._meta = None
            
    def stats(self, chans=None, edge=0, label=""):
        """
        show Mean,RMS,Min,Max
        """
        if chans==None:
            if edge == 0:
                c0 = 0
                c1 = len(self.data)
            else:
                c0 = edge
                c1 = -edge
        else:
            c0 = chans[0]
            c1 = chans[1]
        mean = self.data[c0:c1].mean()
        rms  = self.data[c0:c1].std()
        dmin = self.data[c0:c1].min()
        dmax = self.data[c0:c1].max()
        ndat = c1-c0
        print("%s  %s %s %s %s %d" %  (label,repr(mean),repr(rms),repr(dmin),repr(dmax),ndat))
        return (mean,rms,dmin,dmax,ndat)
    
    def _flux(self, xrange=None, chans=None):
        """
        """
        dx = self.xvals[1]-self.xvals[0]
        if chans==None and xrange==None:
            flux = self.data.sum() * dx
        elif chans != None:
            c0 = chans[0]
            c1 = chans[1]
            flux = self.data[c0:c1].sum() * dx
        elif xrange != None:
            flux = 0
        else:
            flux = -1
        print("Flux %g" % flux)
        
    
    def _plot(self, xrange=None, yrange=None, chans=None, label=None):
        """
        simple spectrum plot
        xrange and yrange are in physical units, e.g. xrange=[4500,5500]
        chans= are in channel numbers, e.g. chans=[2000,3000]
        """
        if xrange != None:   plt.xlim(xrange[0], xrange[1])
        if yrange != None:   plt.ylim(yrange[0], yrange[1])
        if chans == None:
            x = self.xvals
            y = self.data
        else:
            x = self.xvals[chans[0]:chans[1]]
            y = self.data[chans[0]:chans[1]]
        plt.plot(x,y)
        if 'bl' in self.gbt:
            if chans == None:
                bl = self.gbt['bl']
            else:
                bl = self.gbt['bl'][chans[0]:chans[1]]
            #print(type(x),type(bl))
            plt.plot(x,bl,'r-')          # plot full spectrum in red
            # plot the ranges...
            bl_chans = self.gbt['baseline']
            #print(bl_chans)
            for i in range(len(bl_chans)):    # loop over segments to plot in thick black
                c0 = bl_chans[i][0]
                c1 = bl_chans[i][1] 
                plt.plot(x[c0:c1],bl[c0:c1],'k-',linewidth=3)
            
        plt.xlabel("Velocity (km/s)")
        plt.ylabel("$T_A$ (K)")
        title = "Source: %s" % self.meta['OBJECT']
        if label != None:  title = title + " " + label
        plt.title(title)
        
    def _smooth(self, win=11, method='hanning'):
        """
        win needs to be odd !!!
        method = ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']
        """
        data = np.copy(self.data)
        s=np.r_[data[win-1:0:-1],data,data[-2:-win-1:-1]]
        if method == 'flat':     # moving average
            w = np.ones(win,'d')
        else:
            w = eval('np.'+method+'(win)')
        data = np.convolve(w/w.sum(),s,mode='valid')
        w2 = (win-1)//2
        self.data = data[w2:-w2]
        
    def _baseline(self, degree=2, chans=None, replace=False):
        """
        example      chans=[(1000,2000),(5000,6000)]
        """
        poly = models.Polynomial1D(degree=degree)
        fit = fitting.LinearLSQFitter()
        self.gbt['poly_degree'] = degree
        self.gbt['baseline'] = chans
        x = self.xvals
        y = np.copy(self.data)
        if chans != None:
            for i in range(len(chans)):
                c0 = chans[i][0]
                c1 = chans[i][1]
                #print("Segment",c0,c1)
                if i==0:
                    x0 = x[c0:c1]
                    y0 = y[c0:c1]
                else:
                    x0 = np.append(x0,x[c0:c1])
                    y0 = np.append(y0,y[c0:c1])
        else:
            x0 = x
            y0 = y          
        
        bl0 = fit(poly,x0,y0)
        if True:
            # residuals
            bl = bl0(x0) - y0
            print("Residuals[%d] %g %g" % (degree,bl.mean(),bl.std()))
        
        bl = bl0(x)
        #print(type(bl),type(x),type(bl))
        self.gbt['bl'] = bl
        if replace:
            self.data = self.data - bl
        
    def _xy(self):
        """
        return spectrum
        """
        return (self.xvals, self.data)
         


# # SDFITSLoad
# 
# This is the class that loads an SDFITS file. Normally not called by users, but by classes such as GBTLoadPS()
# 
# 

# In[92]:


class SDFITSLoad(object):
    """
    container for a bintable from a selected HDU
    normally not used by users
    """
    def __init__(self, filename, src=None, hdu=None, **kwargs):
        """
        """
        print("==SDFITSLoad %s" % filename)
        kwargs_opts = {'fix':None}
        kwargs_opts.update(kwargs)
        self._filename = filename
        self._bintable = []
        self._binheader = []
        self._data = []
        self._spectra = [] # list of SpecList
        self._hdu = fits.open(filename)  
        self._header = self._hdu[0].header
        self.load(src,hdu,kwargs_opts['fix'])
    @property
    def filename(self):
        return self._filename
    
    def reset(self,hdu=None):
        self._bintable = []
        self._binheader = []
        self._data = []
        self._spectra = []
        self._hdu = fits.open(self._filename)  
        self._header = self._hdu[0].header
        #if hdu is not None:
         #   self.load(src,hdu)
    
    def load(self, src=None, hdu=None, fix=False):
        """
        for given hdu make this bintable available
        """
        self._nrows = []
        if hdu is not None:
            ldu = list([hdu])
        else:
            ldu = range(1,len(self._hdu))
        for i in ldu:
            j=i-1
            self._bintable.append(self._hdu[i])
            self._binheader.append(self._hdu[i].header)
            if src is None:
                self._data.append(self._bintable[j].data)
            else:
                wh1 = self._bintable[j].data[:]['OBJECT'] == src
                if wh1.sum() == 0:
                    srcs = np.unique(self._bintable[j].data[:]['OBJECT'])
                    raise Exception(f"Source name {src} not found in HDU {i}. Sources present are {srcs}")
                self._data.append(self._bintable[j].data[wh1])                                 
            self._nrows.append(len(self._data[j]))
            
    def _loadlists(self,fix=True,dometa=True):
        for i in range(len(self._bintable)):
            sl = SpectrumList()
            for j in range(self.nrows(i)):
                sp = self.rawspectra(i)[j]*u.K
                crval1  = s._data[i]['CRVAL1'][j]
                cdelt1  = s._data[i]['CDELT1'][j]
                crpix1  = s._data[i]['CRPIX1'][j]
                ctype1  = s._data[i]['CTYPE1'][j]
                # Ensure rest frequency is in Hertz
                # CUNIT1 is not always present
                restfrq = s._data[i]['RESTFREQ'][j]
                if "CUNIT1" in s._data[i].names:
                    cunit1  = s._data[i]['CUNIT1'][j]
                    rfq = restfrq * u.Unit(cunit1)
                    restfreq = rfq.to("Hz").value
                cunit1 = "Hz"
                crval2  = s._data[i]['CRVAL2'][j]
                crval3  = s._data[i]['CRVAL3'][j]
                ctype2  = s._data[i]['CTYPE2'][j]
                ctype3  = s._data[i]['CTYPE3'][j]
                # 'FREQ-OBS' to 'FREQ'; assuming SPECSYS='TOPOCENT'
                #if ctype1 == 'FREQ-OBS': ctype1  = 'FREQ'
                # only axis1 needs a full description, axis2,3,4 are all single points
                wcs = WCS(header={'CDELT1': cdelt1, 'CRVAL1': crval1, 'CUNIT1': cunit1,
                                  'CTYPE1': ctype1, 'CRPIX1': crpix1, 'RESTFRQ': restfrq,
                                  'CTYPE2': ctype2, 'CRVAL2': crval2,
                                  'CTYPE3': ctype3, 'CRVAL3': crval3},
                         fix=fix)
                # GBT really fucks up FREQ/VELDEF/VELFRAME
                if False:
                    if "VELFRAME" in s._data[i].names:
                        vframe = s._data[i]['VELFRAME'][j]
                    elif "VFRAME" in s._data[i].names:
                        vframe = s._data[i]['VFRAME'][j]
                    else:
                        vframe = None
                    if "VELDEF" in s._data[i].names:
                        vdef = s._data[i]['VELDEF'][j]
                    else:
                        vdef = None

                    convention = self.velocity_convention(vdef,vframe)
                meta = {}
                convention="doppler_radio"
                if dometa:
                    q = set(s._data[i].names) - set(['DATA'])
                    for n in q:
# FitsRec.__getitem__ [] is VERY EXPENSIVE!
                        meta[n] = s._data[i][n][j]
                    if fix:
                        self.fix_meta(meta)
                    #try:
                    #    convention = self.velocity_convention(meta['VELDEF'],meta['VELFRAME'])
                    #except Exception:
                    #    #print("WARNING: insufficient veldef/velframe, assuming convention is 'doppler_radio'")
                    #    convention="doppler_radio"
                #s1d = Spectrum1D(flux=sp, wcs=wcs, meta=meta, velocity_convention=convention)
                s1d = Spectrum(flux=sp, wcs=wcs, meta=meta, velocity_convention=convention)
                sl.append(s1d)
            s._spectra.append(sl)
    
    def fix_meta(self,meta):
        """Do any repair to the meta/header for peculariaties in definitions from a particular observatory
        The passed-in dictionary will be repaired in place.
        At minimum this method must populate meta['VELDEF'] and meta['VELFRAME']
        """
        pass
    
    def velocity_convention(self,veldef,velframe):
        # sub-classes must implement this so I can keep this class generic.
        # GBT uses VELDEF and VELFRAME incorrectly. 
        return "doppler_radio"
    
    def udata(self,bintable,key):
        return uniq(self._data[bintable][key])
    def ushow(self,bintable,key):
        print(f'{bintable} {key}: {self.udata(bintable,key)}')
        
    def naxis(self,bintable,naxis):
        nax = f'NAXIS{naxis}'
        return self._binheader[bintable][nax]
    
    def nintegrations(self,bintable,source=None):
        if source is not None:
            nint = np.shape(self._data[bintable]['OBJECT'] == source)[0]//self.npol(bintable)
        else:
            nint = np.shape(self._data[bintable])[0]//self.npol(bintable)
        return nint

    def rawspectra(self,bintable):
        return self._data[bintable]['DATA']
    
    def nrows(self,bintable):
        return self._nrows[bintable]

    def nchan(self,bintable):
        return np.shape(self.rawspectra(bintable))[1]
    
    def npol(self,bintable):
        return len(self.udata(bintable,'CRVAL4'))
    
    def sources(self,bintable):
        return self.udata(bintable,'OBJECT')
    
    def scans(self,bintable):
        return self.ushow(bintable,'SCAN')
    
    def __len__(self):
        return self.nrows
    
    def _summary(self,bintable):
        j=bintable
        nrows = self.naxis(j,2)
        nflds = self._binheader[j]['TFIELDS']
        restfreq = np.unique(self._data[j]['RESTFREQ'])/1.0E9
    #
        print("HDU       %d" %  (j+1))
        print("BINTABLE: %d rows x %d cols with %d chans" % (self._nrows[j],nflds,self.nchan(j)))
        print("Selected  %d/%d rows" % (self._nrows[j],nrows))
        print("Sources: ",self.sources(j))
        print("RESTFREQ:",restfreq,'GHz')
        print("Scans:   ",self.scans(j))
        print("Npol:    ",self.npol(j))   
        print("Nint:    ",self.nintegrations(j))
        
    def summary(self):
        print("File:     %s"%self._filename)
        for i in range(len(self._bintable)):
            self._summary(i)
    
    def __repr__(self):
        return self._filename    

class Obsblock():
    def __init__(self,speclist):
        self._speclist = speclist
        #self._summary = 

    def __op__(self,opname):
        pass
    

#class Summary():
#    def __init__(self,bintable):
        


##########
## notes
# Summary can be from an astropy Table or pandas. 
# Same table can be used for data selection.  Can get quite sophisticated if pandas
# SpecList should be encapuslated in ObsBlock class, similar to pyspeckit.
# Need baseline spectraum as part of obsblock class for each Spectrum1d
# this argues for a gbt.Spectrum class that contains the baseline, which is I think what Peter did
# Also need 'smooth spectrum'.  undo('baseline') or undo('smooth') resets these to None respectively print('undid %s',*arg) hehehe. 
# Does spec1d have weights? _default_weigts
#######

# In[93]:


if False:
    s = SDFITSLoad('ngc5291.fits',fix=False)
    s.summary()
    get_size(s)
    s._loadlists(fix=False)
    get_size(s)


    fname="/bigdisk/data/gbt/examples/mapping-L/data/TGBT17A_506_11.raw.vegas/TGBT17A_506_11.raw.vegas.A.fits"
    s = SDFITSLoad(fname,src="3C286",hdu=1)

    s.summary()
    get_size(s)

    s = SDFITSLoad(fname)
    s.summary()
    get_size(s)


# # GBTLoad
# 
# This is the base class from which we derive all GBTLoad* subclasses that can load and calibrate spectra. It can also be used to load all the spectra, but there is no structure defined, e.g. to check and guide calibration manually.


class GBTLoad(SDFITSLoad):
    def __init__(self, filename, src=None,hdu=None):
        """
        Holds a raw "unstructured" series of scans, normally not used by users
        """       
        SDFITSLoad.__init__(self,filename,src,hdu,fix=False)
        print("==GBTLoad %s" % filename)

        self.ushow(0,'OBJECT')
        self.ushow(0,'SCAN')
        self.ushow(0,'SAMPLER')
        #ushow('PLNUM')
        #ushow('IFNUM')
        self.ushow(0,'SIG')
        self.ushow(0,'CAL')
        self.ushow(0,'PROCSEQN')
        self.ushow(0,'PROCSIZE')
        self.ushow(0,'OBSMODE')  
        self.ushow(0,'SIDEBAND')

# PS example      
if False:
    ex1 = GBTLoad('ngc5291.fits')
    ex1.summary()
    get_size(ex1)

if __name__ == "__main__":
    examples = "/data/gbt/examples/"
    files = [
    "mapping-L/data/TGBT17A_506_11.raw.vegas/TGBT17A_506_11.raw.vegas.A.fits", 
    "onoff-L/data/TGBT21A_501_11.raw.vegas.fits", 
    "rxco-W/data/TSCAL_220105_W.raw.vegas/TSCAL_220105_W.raw.vegas.A.fits", 
    "rxco-W/data/TSCAL_220105_W.raw.vegas/TSCAL_220105_W.raw.vegas.E.fits",
    "misc/ngc5291.fits",
    "misc/W3OH.fits",
    "misc/IC1481.fits",
    ]
    timestr = ""
    for fn in files[3:]:
        pr = cProfile.Profile()
        inf = f'{examples}{fn}'
        size = os.path.getsize(inf)/1048576
        #print(f'\n{inf} size: {size:.1f} MB')
        t0 = time.perf_counter_ns()
        pr.enable()
        s = SDFITSLoad(inf)
        t1 = time.perf_counter_ns()
        s._loadlists(fix=False,dometa=True)
        pr.disable()
        t2 = time.perf_counter_ns()
        #s.summary()
        nhdu = len(s._hdu)
        nrows = np.sum(s._nrows)

        #s = io.StringIO()
        sortby = SortKey.CUMULATIVE
        ps = pstats.Stats(pr).sort_stats(sortby)
        ps.print_stats(10)

        timestr+=f'   {size:.0f}\t{nhdu}\t{nrows}\t{(t1-t0)/1E6:.1f}\t{(t2-t1)/1E6:.1f}\t{(t2-t0)/1E6:.1f}\n'
    print('\n# Timing (ms)')
    print('# Size\tNhdu\tNrows\tLoad\tLoad_Lists\tTotal')
    print(timestr)

    #print(s.getvalue())

