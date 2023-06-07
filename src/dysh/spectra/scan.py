import astropy.units as u
from copy import deepcopy
import numpy as np
from astropy.wcs import WCS
from .spectrum import Spectrum
from . import dcmeantsys

class PSScan(object):
    """
    Holds a position scan pair

    Parameters
    ----------
        sdfits : ~SDFITSLoad
            input SDFITSLoad object (or derivative)
        scans : dict
            dictionary with keys 'ON' and 'OFF' containing unique list of ON (sig) and OFF (ref) scan numbers
        scanrows : dict
            dictionary with keys 'ON' and 'OFF' containing the list of rows in `sdfits` corresponding to ON (sig) and OFF (ref) integrations
        bintable : int
            the index for BINTABLE in `sdfits` containing the scans
    """
    def __init__(self, sdfits, scans, scanrows, bintable):
        self._sdfits = sdfits # parent class
        self._status = 0 #@TODO make these an enumeration, possibly dict
        #                           # ex1:
        self._nint = 0               # 11
        self._npol = 0               #  2
        self._on = None              # 44
        self._off = None             # 44
        self._calibrated = None      # 22
        self._timeaveraged = None    #  2
        self._polaveraged = None     #  1
        self._bintable_index = bintable
        self._nrows = len(scanrows['ON'])
        print(f"PSSCAN nrows = {self.nrows}")
        
    @property
    def status(self):
        """Status flag, will be used later for undo"""
        return self._status

    @property
    def nrows(self):
        """The number of rows in this Scan"""
        return self._nrows

    @property
    def npol(self):
        """The number of polarizations in this Scan"""
        return self._npol

    def __len__(self):
        return self._nrows
    
class GBTPSScan(PSScan):
    """GBT specific version of Position Switch Scan (PSScan)
    Parameters
    ----------
        sdfits : ~GBFITSLoad
            input GBFITSLoad object 
        scans : dict
            dictionary with keys 'ON' and 'OFF' containing unique list of ON (sig) and OFF (ref) scan numbers
        scanrows : dict
            dictionary with keys 'ON' and 'OFF' containing the list of rows in `sdfits` corresponding to ON (sig) and OFF (ref) integrations
        calrows : dict
            dictionary containing with keys 'ON' and 'OFF' containing list of rows in `sdfits` corresponding to cal=T (ON) and cal=F (OFF) integrations. 
        bintable : int
            the index for BINTABLE in `sdfits` containing the scans
    """
    def __init__(self, gbtfits, scans, scanrows, calrows, bintable=0):
        PSScan.__init__(self,gbtfits,scans,scanrows,bintable)
        # The rows of the original bintable corresponding to ON (sig) and OFF (reg)
        self._scanrows = scanrows
        #print(f"scanrows ON {self._scanrows['ON']}")
        #print(f"scanrows OFF {self._scanrows['OFF']}")
        
        # calrows perhaps not needed as input since we can get it from gbtfits object?
        #calrows['ON'] are rows with noise diode was on, regardless of sig or ref
        #calrows['OFF'] are rows with noise diode was off, regardless of sig or ref
        self._calrows = calrows
        self._npol =  gbtfits.npol(bintable) #TODO deal with bintable
        self._nint = gbtfits.nintegrations(bintable)
        # todo use gbtfits.velocity_convention(veldef,velframe)
        vc = "doppler_radio"
        # so quick with slicing!
        self._sigonrows = list(set(self._calrows["ON"]).intersection(set(self._scanrows["ON"])))
        self._sigoffrows = list(set(self._calrows["OFF"]).intersection(set(self._scanrows["ON"])))
        self._refonrows = list(set(self._calrows["ON"]).intersection(set(self._scanrows["OFF"])))
        self._refoffrows = list(set(self._calrows["OFF"]).intersection(set(self._scanrows["OFF"])))
        self._sigcalon = gbtfits.rawspectra(bintable)[self._sigonrows]
        self._sigcaloff = gbtfits.rawspectra(bintable)[self._sigoffrows]
        self._refcalon = gbtfits.rawspectra(bintable)[self._refonrows]
        self._refcaloff = gbtfits.rawspectra(bintable)[self._refoffrows]
        self._tsys = None

    @property
    def tsys(self):
        """The system temperature array. This will be `None` until 
           calibration is done

        Returns
        -------
            tsys : ~numpy.ndarray 
                System temperature values in K
        """
        return self._tsys

    #TODO something clever 
    # self._calibrated_spectrum = Spectrum(self._calibrated,...) [assuming same spectral axis]
    def calibrated(self,i):
        """Return the calibrated Spectrum
        
        Parameters
        ---------
            i : int
                The index into the calibrated array

        Returns
        -------
            spectrum : ~Spectrum
        """
        meta = dict(self._sdfits.index(self._bintable_index).iloc[self._scanrows["ON"][i]])
        meta['TSYS'] = self._tsys[i]
        naxis1 = len(self._calibrated[i])
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


        
        vc = "doppler_radio"
        
        return Spectrum(self._calibrated[i]*u.K,wcs=wcs,meta=meta,velocity_convention=vc)

    def calibrate(self, **kwargs):
        """
        Position switch calibration, following equations 1 and 2 in the GBTIDL clibration manual
        """
        kwargs_opts = {
            'verbose': False
        }
        kwargs_opts.update(kwargs)

        self._status = 1
        #npolint = self.npol * self.nint
        nspect = self.nrows//2
        self._calibrated = np.empty(nspect, dtype=np.ndarray)
        self._tsys= np.empty(nspect, dtype=float)

        tcal = list(self._sdfits.index(self._bintable_index).iloc[self._refonrows]["TCAL"])
        #@Todo  this loop could be replaces with clever numpy
        if len(tcal) != nspect: 
            raise Exception(f"TCAL length {len(tcal)} and number of spectra {nspect} don't match")
        for i in range(nspect):
            #tcal = self.off[2*i].meta['TCAL']
            #tcal2= self.on[2*i].meta['TCAL']
            tsys = dcmeantsys(calon=self._refcalon[i],  caloff=self._refcaloff[i],tcal=tcal[i])
            #tsys2= dcmeantsys(self.sigcalon[i],  self.sigcalon[i],tcal2)
            #if kwargs_opts['verbose']: 
            #    print(i,tcal,tsys,tsys2)
            #                 2*i is the CalON     2*i+1 the CalOFF
            #
            #calon = 2*i
            #caloff = calon+1
            sig = 0.5*(self._sigcalon[i] + self._sigcaloff[i])
            ref = 0.5*(self._refcalon[i] + self._refcaloff[i])
            self._calibrated[i] = tsys * (sig-ref) / ref
            #self.calibrated[i].gbt['row'] = kr
            self._tsys[i] = tsys
