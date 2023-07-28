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
from ..spectra import veldef_to_convention
from ..util import uniq, stripTable
#from .. import version


class SDFITSLoad(object):
    '''
    Generic Container for a bintable(s) from selected HDU(s)

    Parameters
    ----------
        filename : str 
            input file name
        source  : str
            target source to select from input file. Default: all sources
        hdu : int or list
            Header Data Unit to select from input file. Default: all HDUs
    '''
    def __init__(self, filename, source=None, hdu=None, **kwargs):
        if kwargs.get("verbose",None):
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
        self._hdu = fits.open(filename)  
        self._header = self._hdu[0].header
        self.load(hdu,**kwargs_opts)
        self.create_index()
        #self._hdu.close()  # can't access hdu[i].data member of you do this.

    def info(self):
        """Return the `~astropy.HDUList` info()"""
        return self._hdu.info()

    @property 
    def bintable(self):
        """The list of bintables"""
        return self._bintable

    def binheader(self):
        """The list of bintable headers"""
        return self._binheader

    @property 
    def filename(self):
        """The input SDFITS filename"""
        return self._filename
    
    def index(self,hdu):
        """The index table"""
        return self._ptable[hdu]

    def reset(self,hdu=None):
        """Reset all attributes"""
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
        """Create the index of the SDFITS file.

        Parameters
        ----------
            hdu : int or list
                Header Data Unit to select from input file. Default: all HDUs
        """
        if hdu is not None:
            ldu = list([hdu])
        else:
            ldu = range(1,len(self._hdu))
        self._ptable = []
        for i in ldu:
            t = Table.read(self._hdu[i]) 
            t.remove_column('DATA')
            stripTable(t)
            self._ptable.append(t.to_pandas())
            del t
            
    def load(self, hdu=None, **kwargs):
        """
        Load the bintable for given hdu.
        Note mmHg and UTC are unrecognized units.  mmHg is in astropy.units.cds but UTC is just wrong.

        Parameters
        ----------
            hdu : int or list
                Header Data Unit to select from input file. Default: all HDUs
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

    def _loadlists(self,hdu,fix=False,wcs=False,maxspect=1E16):
        '''Create an obsblock from all rows in bintable.  For debug/performance testing only'''
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
        """Do any repair to the meta/header for peculariaties in definitions 
        from a particular observatory
        The passed-in dictionary will be repaired in place.
        At minimum this method must populate meta['VELDEF'] and meta['VELFRAME']

        Parameters
        ----------
            meta : dict
                The header of the `~Spectrum` to be fixed, corresponding to the `meta` attribute of the Spectrum.
        """
        pass
    
    def velocity_convention(self,veldef,velframe):
        '''Compute the velocity convention string use for velocity conversions, 
        given the VELDEF and VELFRAME values. 
        Return value must be a recognized string of `~specutils.Spectrum1D`, one of
        {"doppler_relativistic", "doppler_optical", "doppler_radio"}
        Sub-classes should implement, because different observatories use VELDEF and 
        VELFRAME inconsistently. This base class method hard-coded to return "doppler_radio."
            
        Parameters
        ----------
            veldef : str
                The velocity definition string (`VELDEF` FITS keyword)
            velframe : str
                The velocity frame string (`VELFRAME` FITS keyword)
        '''
        return "doppler_radio"
    
    def udata(self,bintable,key):
        """The unique list of values of a given header keyword

        Parameters
        ----------
            bintable :  int
                The index of the `bintable` attribute
            key : str
                The keyword to retrieve
        Returns
        -------
            udata : list
                The unique set of values for the input keyword.
        """
        return uniq(self._ptable[bintable][key])
                                  
    def ushow(self,bintable,key):
        """Print the unique list of values of a given header keyword

        Parameters
        ----------
            bintable :  int
                The index of the `bintable` attribute
            key : str
                The keyword to retrieve
        """
        print(f'{bintable} {key}: {self.udata(bintable,key)}')
        
    def naxis(self,bintable,naxis):
        '''The NAXISn value of the input bintable.

        Parameters
        ----------
            bintable :  int
                The index of the `bintable` attribute
            naxis : int
                The NAXIS whose length is requested

        Returns
        -------
            naxis : the length of the NAXIS
        '''
        nax = f'NAXIS{naxis}'
        return self._binheader[bintable][nax]
    
    def nintegrations(self,bintable,source=None):
        '''The number of integrations on a given source divided by the number of polarizations

        Parameters
        ----------
            bintable :  int
                The index of the `bintable` attribute
            source: str
                The source name (OBJECT keyword) or None for all sources. Default: None

        Returns
        -------
            nintegrations : the number of integrations 
        '''
        
        data = self.rawspectra(bintable)
        if source is not None:
            df = self.select('OBJECT',source,self._ptable[bintable])
            #nfeed = df["FEED"].nunique()
            numsources = len(df)
            #nint = numsources//(self.npol(bintable)*nfeed)
            nint = numsources//self.npol(bintable)
        else:
            nint = self.nrows(bintable)//self.npol(bintable)
        return nint

    def rawspectra(self,bintable):
        '''Get the raw (unprocessed) spectra from the input bintable.

        Parameters
        ----------
            bintable :  int
                The index of the `bintable` attribute

        Returns
        -------
            rawspectra : ~numpy.ndarray
                The DATA column of the input bintable
        '''
        return self._bintable[bintable].data[:]["DATA"]

    def rawspectrum(self,bintable,i):
        '''Get a single raw (unprocessed) spectrum from the input bintable.
    TODO: arguments are backwards from getrow(), getspec()

        Parameters
        ----------
            bintable :  int
                The index of the `bintable` attribute
            i :  int
                The row index to retrieve

        Returns
        -------
            rawspectrum : ~numpy.ndarray
                The i-th row of DATA column of the input bintable
        '''
        return self._bintable[bintable].data[:]["DATA"][i]

    def getrow(self,i,bintable=0):
        return self._bintable[bintable].data[i]

    def getspec(self,i,bintable=0):
        """get a row (record) as a Spectrum"""
        meta = self._ptable[bintable].iloc[i]
        data = self.rawspectrum(bintable,i)
        naxis1 = len(data)
        ctype1 = meta['CTYPE1']
        ctype2 = meta['CTYPE2']
        ctype3 = meta['CTYPE3']
        crval1 = meta['CRVAL1']
        crval2 = meta['CRVAL2']
        crval3 = meta['CRVAL3']
        crpix1 = meta['CRPIX1']
        cdelt1 = meta['CDELT1']
        restfrq = meta['RESTFREQ']
        if 'CUNIT1' in meta:
            cunit1 = meta['CUNIT1']
        else:
            cunit1 = "Hz" #@TODO this is in gbtfits.hdu[0].header['TUNIT11'] but is it always TUNIT11?
        rfq = restfrq * u.Unit(cunit1)
        restfreq = rfq.to("Hz").value

        #@TODO WCS is expensive.  Figure how to calculate spectral_axis instead.
        wcs = WCS(header={'CDELT1': cdelt1, 'CRVAL1': crval1, 'CUNIT1': cunit1,
                                      'CTYPE1': 'FREQ', 'CRPIX1': crpix1, 'RESTFRQ': restfreq,
                                      'CTYPE2': ctype2, 'CRVAL2': crval2, 'CRPIX2': 1,
                                      'CTYPE3': ctype3, 'CRVAL3': crval3, 'CRPIX3': 1,
                                      'CUNIT2': 'deg', 'CUNIT3':'deg',
                                      'NAXIS1': naxis1, 'NAXIS2':1, 'NAXIS3':1
                                     },
                             )
        vc = veldef_to_convention(meta['VELDEF'])
        
        # raw data are in counts
        return Spectrum(data*u.count,wcs=wcs,meta=meta.to_dict(),velocity_convention=vc)

    
    def nrows(self,bintable):
        '''The number of rows of the input bintable

        Parameters
        ----------
            bintable :  int
                The index of the `bintable` attribute

        Returns
        -------
            nrows : int
                Number of rows, i.e., the length of the input bintable
        '''
        return self._nrows[bintable]

    def nchan(self,bintable):
        '''The number of channels per row of the input bintable. Assumes all rows have same length.

        Parameters
        ----------
            bintable :  int
                The index of the `bintable` attribute

        Returns
        -------
            nchan : int
                Number channels in the first spectrum of the input bintbale
        '''
        return np.shape(self.rawspectrum(bintable,1))[0]
    
    def npol(self,bintable):
        '''The number of polarizations present in the input bintable. 

        Parameters
        ----------
            bintable :  int
                The index of the `bintable` attribute

        Returns
        -------
            npol: int
                Number of polarizations as given by `CRVAL4` FITS header keyword.
        '''
        return len(self.udata(bintable,'CRVAL4'))
    
    def sources(self,bintable):
        '''The number of sources present in the input bintable. 

        Parameters
        ----------
            bintable :  int
                The index of the `bintable` attribute

        Returns
        -------
            sources: int
                Number of sources as given by `OBJECT` FITS header keyword.
        '''
        return self.udata(bintable,'OBJECT')
    
    def scans(self,bintable):
        #@TODO move this to GBTFISLoad?
        '''The number of scans resent in the input bintable. 

        Parameters
        ----------
            bintable :  int
                The index of the `bintable` attribute

        Returns
        -------
            scans: int
                Number of scans as given by `SCAN` FITS header keyword.
        '''
        return self.udata(bintable,'SCAN') #self.ushow(bintable,'SCAN')
    
    def _summary(self,bintable):
        j=bintable
        nrows = self.naxis(j,2)
        nflds = self._binheader[j]['TFIELDS']
        restfreq = np.unique(self._ptable[j]['RESTFREQ'])/1.0E9
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
    
    def __len__(self):
        return self.nrows
    
    def __repr__(self):
        return self._filename    

    def _bintable_from_rows(self,rows=None,bintable=None):
        """extract a bintable from an existing 
           bintable in this SDFITSLoad object
            
        Parameters
        ----------

            rows: int or list-like
                Range of rows in the bintable(s) to write out. e.g. 0, [14,25,32]. 
            bintable :  int
                The index of the `bintable` attribute

        Returns
        -------
            outbintable : ~astropy.iofits.BinTableHDU
                The binary table HDU created containing the selected rows

        """
        # ensure it is a list if int was given 
        if type(rows) == int:
            rows = [rows]
        outbintable = self._bintable[bintable].copy()
        #print(f"bintable copy data length {len(outbintable.data)}")
        outbintable.data = outbintable.data[rows]
        #print(f"bintable rows data length {len(outbintable.data)}")
        outbintable.update()
        return outbintable

    def write(self,fileobj,rows=None,bintable=None,output_verify="exception",overwrite=False,checksum=False):
        """
        Write the `SDFITSLoad` to a new file, potentially sub-selecting rows or bintables.

        Parameters
        ----------
        fileobj : str, file-like or `pathlib.Path`
            File to write to.  If a file object, must be opened in a
            writeable mode.

        rows: int or list-like
            Range of rows in the bintable(s) to write out. e.g. 0, [14,25,32]. Default: None, meaning all rows
            Note: Currently `rows`, if given, must be contained in a single bintable and bintable must be given

        bintable :  int
            The index of the `bintable` attribute or None for all bintables. Default: None

        output_verify : str
            Output verification option.  Must be one of ``"fix"``,
            ``"silentfix"``, ``"ignore"``, ``"warn"``, or
            ``"exception"``.  May also be any combination of ``"fix"`` or
            ``"silentfix"`` with ``"+ignore"``, ``+warn``, or ``+exception"
            (e.g. ``"fix+warn"``).  See https://docs.astropy.org/en/latest/io/fits/api/verification.html for more info

        overwrite : bool, optional
            If ``True``, overwrite the output file if it exists. Raises an
            ``OSError`` if ``False`` and the output file exists. Default is
            ``False``.

        checksum : bool
            When `True` adds both ``DATASUM`` and ``CHECKSUM`` cards
            to the headers of all HDU's written to the file.
        """
        if bintable is None:
            if rows is None:
                # write out everything 
                self._hdu.writeto(fileobj,
                    output_verify=output_verify,
                    overwrite=overwrite,checksum=checksum)
            else:
                raise ValueError("You must specify bintable if you specify rows")
        else:
            if rows is None:
                # bin table index counts from 0 and starts at the 2nd HDU (hdu index 1), so add 2
                self._hdu[0:bintable+2]._hdu.writeto(fileobj,
                    output_verify=output_verify, 
                    overwrite=overwrite, checksum=checksum)
            else:
                hdu0 = self._hdu[0].copy()
                # need to get imports correct first
                #hdu0.header["DYSHVER"] = ('dysh '+version(), "This file was created by dysh")
                outhdu = fits.HDUList(hdu0)
                outbintable = self._bintable_from_rows(rows,bintable)
                outhdu.append(outbintable)
                outhdu.writeto(fileobj,
                    output_verify=output_verify,
                    overwrite=overwrite, checksum=checksum)
