"""Load generic SDFITS files
    - Not typically used directly.  Sub-class for specific telescope SDFITS flavors.
"""
import sys
import copy
from astropy.wcs import WCS
from astropy.units import cds
from astropy.io import fits
from astropy.modeling import models, fitting
import astropy.units as u
from astropy.table import Table
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from ..spectra.spectrum import Spectrum
from ..spectra.obsblock import Obsblock
from ..spectra import dcmeantsys
from ..util import uniq


class SDFITSLoad(object):
    '''
    Generic Container for a bintable(s) from selected HDU(s)

    Parameters
    -----------
    filename: - str 
        input file name
    source: - str
        target source to select from input file. Default: all sources
    hdu: - int or list
        Header Data Unit to select from input file. Default: all HDUs
    '''
    def __init__(self, filename, source=None, hdu=None, **kwargs):
        print("==SDFITSLoad %s" % filename)
        cds.enable()  # to get mmHg
        kwargs_opts = {'fix':False}
        kwargs_opts = {'wcs':False}
        kwargs_opts.update(kwargs)
        self._filename = filename
        self._bintable = []
        self._ptable = []
        self._binheader = []
        self._data = []
        self._obsblock = [] # list of SpecList
        self._hdu = fits.open(filename)  
        self._header = self._hdu[0].header
        self.load(hdu,**kwargs_opts)
        self.create_index()
        #self._hdu.close()  # can't access hdu[i].data member of you do this.

    @property 
    def filename(self):
        return self._filename
    
    @property
    def index(self,hdu):
        return self._ptable[hdu]

    def reset(self,hdu=None):
        self._bintable = []
        self._binheader = []
        self._ptable = []
        self._data = []
        self._spectra = []
        self._hdu = fits.open(self._filename)  
        self._header = self._hdu[0].header
        #if hdu is not None:
         #   self.load(src,hdu)
    
    def create_index(self,hdu=None):
        if hdu is not None:
            ldu = list([hdu])
        else:
            ldu = range(1,len(self._hdu))
        self._ptable = []
        for i in ldu:
            t = Table.read(self._hdu[i]) 
            t.remove_column('DATA')
            self.stripT(t)
            print(f"doing pandas for HDU {i}")
            self._ptable.append(t.to_pandas())
            del t
            
    def load(self, hdu=None, **kwargs):
        """
        for given hdu make this bintable available
        Note mmHg and UTC are unrecognized units.  mmHg is in astropy.units.cds but UTC is just wrong.
        """

        self._bintable = []
        self._ptable = []
        self._binheader = []
        self._data = []
        self._nrows = []
        source = kwargs.get('source',None)
        fix = kwargs.get('fix',False)
        wcs = kwargs.get('wcs',False)

        if hdu is not None:
            ldu = list([hdu])
        else:
            ldu = range(1,len(self._hdu))
        for i in ldu:
            j=i-1
            self._bintable.append(self._hdu[i]) 
            self._binheader.append(self._hdu[i].header)
            # TODO: don't allow preselection here, it just screws bookkeepingu p.
            # All selection must happen in obsblocks.
            #if source is None:
            #self._data.append(self._bintable[j].data[:]["DATA"])
            #else:
            #    wh1 = np.char.strip(self._bintable[j]['OBJECT']) == source # true/false array
            #    if wh1.sum() == 0: #all False none found
            #        srcs = np.unique(self._bintable[j]['OBJECT'])
            #        raise Exception(f"Source name {source} not found in HDU {i}. Sources present are {srcs}")
            #    self._data.append(self._bintable[j]["DATA"][wh1])
            #    mask = np.where(wh1)[0]
            #    self._bintable[j] = self._bintable[j][mask]  # header will be wrong?
            self._nrows.append(self._binheader[j]["NAXIS2"])

        if kwargs.get("index",False):
            self.create_index(hdu)

    def stripT(self,b):
        # remove leading and trailing chars from all strings in table 
        for n in b.colnames:
            if np.issubdtype(b.dtype[n],str):
                b[n] = np.char.strip(b[n])

    def strip(self):
        # remove leading and trailing chars from all strings in table 
        for b in self._ptable:
            for n in b.colnames:
                if np.issubdtype(b.dtype[n],str):
                    b[n] = np.char.strip(b[n])

    def _loadlists(self,hdu,fix=False,wcs=False,maxspect=1E16):
        self._obsblock = []
        i=0
        k = -1
        print("HDU = ",hdu)
        if hdu is not None:
            b = self._ptable[hdu-1]
            rawspect = self._bintable[i].data["DATA"]
            sl = SpectrumList()
            maxload = int(np.min([maxspect,self.nrows(i)]))
            print(f"Creating {maxload} Spectrum in bintable {i} HDU {hdu}",file=sys.stderr)
            for j in range(maxload):
                k = k+1
                # need extra [[]] because we have 1x1 spatial NAXIS
                # otherwise, slicing the spectrum won't work.
                if wcs: 
                    sp = np.array([[self.rawspectrum(i,j)]])
                else:
                    #sp = self.rawspectrum(i,j)*u.K
                    sp = np.copy(rawspect[j])#*u.K
                    #sp = np.random.rand(32768)*u.K
                naxis1 =  sp.shape[0]#self.nchan(i)
                printme = int(0.1*len(b))
                if (k%printme) == 0: 
                    print(f"Row {k} nchan {naxis1} {type(sp)}", file=sys.stderr)
                    #print(f"NAXIS1 is {naxis1}",file=sys.stderr)
                crval1  = b['CRVAL1'][j]
                cdelt1  = b['CDELT1'][j]
                crpix1  = b['CRPIX1'][j]
                ctype1  = b['CTYPE1'][j]
                # Ensure rest frequency is in Hertz
                # CUNIT1 is not always present
                restfrq = b['RESTFREQ'][j]
                if "CUNIT1" in b.columns:
                    cunit1  = b['CUNIT1'][j]
                    rfq = restfrq * u.Unit(cunit1)
                    restfreq = rfq.to("Hz").value
                cunit1 = "Hz"
                crval2  = b['CRVAL2'][j]
                crval3  = b['CRVAL3'][j]
                ctype2  = b['CTYPE2'][j]
                ctype3  = b['CTYPE3'][j]
                # 'FREQ-OBS' to 'FREQ'; assuming SPECSYS='TOPOCENT'
                #if ctype1 == 'FREQ-OBS': ctype1  = 'FREQ'
                # only axis1 needs a full description, axis2,3,4 are all single points
                if wcs:
                    wcs = WCS(header={'CDELT1': cdelt1, 'CRVAL1': crval1, 'CUNIT1': cunit1,
                                      'CTYPE1': 'FREQ', 'CRPIX1': crpix1, 'RESTFRQ': restfrq,
                                      'CTYPE2': ctype2, 'CRVAL2': crval2, 'CRPIX2': 1,
                                      'CTYPE3': ctype3, 'CRVAL3': crval3, 'CRPIX3': 1,
                                      'CUNIT2': 'deg', 'CUNIT3':'deg',
                                      'NAXIS1': naxis1, 'NAXIS2':1, 'NAXIS3':1
                                     },
                             fix=fix)
                else:
                    wcs = None
                # GBT really fucks up FREQ/VELDEF/VELFRAME
                #if False:
                if "VELFRAME" in b.columns:
                    vframe = b['VELFRAME'][j]
                elif "VFRAME" in b.columns:
                    vframe = b['VFRAME'][j]
                else:
                    vframe = None
                if "VELDEF" in b.columns:
                    vdef = b['VELDEF'][j]
                else:
                    vdef = None

                meta={'CDELT1': cdelt1, 'CRVAL1': crval1, 'CUNIT1': cunit1,
                      'CTYPE1': 'FREQ', 'CRPIX1': crpix1, 'RESTFRQ': restfrq,
                      'CTYPE2': ctype2, 'CRVAL2': crval2, 'CRPIX2': 1,
                      'CTYPE3': ctype3, 'CRVAL3': crval3, 'CRPIX3': 1,
                      'CUNIT2': 'deg', 'CUNIT3':'deg',
                      'NAXIS1': naxis1, 'NAXIS2':1, 'NAXIS3':1,
                      'VELDEF': vdef, 'VELFRAME': vframe
                     }

                convention = self.velocity_convention(vdef,vframe)
                #meta = dict(b.loc[j])# Necessary? Since we are sending whole pandas table to Obsblock
                if fix:
                    self.fix_meta(meta)
                if False:
                    try:
                        convention = self.velocity_convention(meta['VELDEF'],meta['VELFRAME'])
                    except Exception:
                        #print("WARNING: insufficient veldef/velframe, assuming convention is 'doppler_radio'")
                        convention="doppler_radio"
                meta = {}
                sl.append(Spectrum(flux=sp*u.K,wcs=wcs,  meta=meta, velocity_convention=convention))
            self._obsblock.append(Obsblock(sl,self._ptable[i]))
            i=i+1


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
        return uniq(self._ptable[bintable][key])
                                  
    def ushow(self,bintable,key):
        print(f'{bintable} {key}: {self.udata(bintable,key)}')
        
    def naxis(self,bintable,naxis):
        nax = f'NAXIS{naxis}'
        return self._binheader[bintable][nax]
    
    def nintegrations(self,bintable,source=None):
        if source is not None:
            nint = np.shape(np.char.strip(self._data[bintable]['OBJECT']) == source)[0]//self.npol(bintable)
        else:
            nint = np.shape(self._data[bintable])[0]//self.npol(bintable)
        return nint

    def rawspectra(self,bintable):
        return self._bintable[bintable].data[:]["DATA"]

    def rawspectrum(self,bintable,i):
        return self._bintable[bintable].data[:]["DATA"][i]
    
    def nrows(self,bintable):
        return self._nrows[bintable]

    def nchan(self,bintable):
        return np.shape(self.rawspectrum(bintable,1))[0]
    
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
        restfreq = np.unique(self._ptable['RESTFREQ'])/1.0E9
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
        """Print a summary of each record of the data"""
        print("File:     %s"%self._filename)
        for i in range(len(self._bintable)):
            print("i=",i)
            self._summary(i)
    
    def __repr__(self):
        return self._filename    

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

