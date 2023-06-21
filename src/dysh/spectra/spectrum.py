from specutils import Spectrum1D, SpectrumList,SpectralRegion
from astropy.modeling.polynomial import Polynomial1D,Chebyshev1D
from astropy.modeling.fitting import LevMarLSQFitter,LinearLSQFitter
from astropy.wcs import WCS
import astropy.units as u
from . import baseline

class Spectrum(Spectrum1D):
    """This class contains a spectrum and its attributes. It is built on
    `~specutils.Spectrum1D` with added attributes like baseline model.
    Note that `~specutils.Spectrum1D` can contain multiple spectra but
    we probably will not use that because the restriction that it can
    have only one spectral axis conflicts with slight Doppler shifts.
    See `~specutils.Spectrum1D` for the instantiation arguments.

    *Note:* `velocity_convention` should be one of {'radio', 'optical', 'relativistic'}; the  `~specutils.Spectrum1D` is wrong (there should not be a 'doppler_' prefix).
    """
    def __init__(self, *args,**kwargs):
        Spectrum1D.__init__(self,*args,**kwargs)
        self._baseline_model = None

    @property
    def baseline_model(self):
        """Returns the computed baseline model or None if it has not yet been computed."""
        return self._baseline_model

    def baseline(self,degree,exclude=None,**kwargs):
        """Compute and optionally remove a baseline.  The model for the
        baseline can be either a 
        `1D polynomial model <https://docs.astropy.org/en/latest/api/astropy.modeling.polynomial.Polynomial1D.html>`_ or a 
        `1D Chebyshev polynomial of the first kind <https://docs.astropy.org/en/latest/api/astropy.modeling.polynomial.Chebyshev1D.html>`_.  The code uses `astropy.modeling`
        and `astropy.fitter` to compute the baseline.  See the documentation for those modules.  This method will set the `baseline_model` attribute to the fitted model function which can be evaluated over a domain.
        
        Parameters
        ----------
            degree : int
                The degree of the polynomial series, a.k.a. baseline order
            exclude : list of 2-tuples of int or ~astropy.units.Quantity, or ~specutils.SpectralRegion
                List of region(s) to exclude from the fit.  The tuple(s) represent a range in the form [lower,upper], inclusive.  
in channel units.  

                Examples: One channel-based region: [11,51], Two channel-based regions: [(11,51),(99,123)]. One ~astropy.units.Quantity region: [110.198*u.GHz,110.204*u.GHz]. One compound ~specutils.SpectralRegion: SpectralRegion([(110.198*u.GHz,110.204*u.GHz),(110.196*u.GHz,110.197*u.GHz)]).

                Default: no exclude region

                TODO: Are these OR'd with the existing mask? make that an option
                TODO: Allow these to be Quantities (spectral axis units or equivalent). See list_to_spectral_region()
            model : str
                One of 'polynomial' or 'chebyshev', Default: 'polynomial'
            fitter  :  `~astropy.fitting._FitterMeta`
                The fitter to use. Default: `~astropy.fitter.LinearLSQFitter` (with `calc_uncertaintes=True).  Be care when choosing a different fitter to be sure it is optimized for this problem.
            remove : bool
                If True, the baseline is removed from the spectrum. Default: False
        """
        kwargs_opts = {
            'remove':False,
            'model':'polynomial',
            'fitter':  LinearLSQFitter(calc_uncertainties=True),
        }
        kwargs_opts.update(kwargs)

        self._baseline_model = baseline(self,degree,exclude,**kwargs)
        if kwargs_opts['remove']:
            s=self.subtract(self._baseline_model(self.spectral_axis))
            self._data = s._data

    def _undo_baseline(self):
        """Undo the most recently computed baseline.  If the baseline
           has been subtracted, it will be added back.  The `baseline_model`
           attribute is set to None.  
        """
        if self._baseline_model is None:
            return
        s=self.add(self._baseline_model(self.spectral_axis))
        self._data = s._data
        self._baseline_model = None

    def _set_exclude_regions(self, exclude):
        """Set the mask for the regions to exclude.
        
        Parameters
        ----------
            exclude : `~numpy.ndarray`-like
                Array where values in the flux to be masked are those that
                astype(bool) converts to True.  
        """
        pass

    def list_to_spectral_region(self,inlist):
        #todo utility code to convert a input list of channels or quantities to a spectral region with units of self.spectral_axis.unit. This could go in core.py
        # combine this with _set_exclude_regions
        pass

    def bshow(self):
        """Show the baseline model""" 
        print(f"baseline model {self._baseline_model}")

    def stats(self):
        """ Compute some statistics of this `Spectrum`.  The mean, rms,
        data minimum and data maximum are calculated.  Note this works
        with slicing, so, e.g.,  `myspectrum[45:153].stats()` will return
        the statistics of the slice.  

        Returns
        -------
        stats : tuple 
            Tuple consisting of (mean,rms,datamin,datamax)
            TODO: maybe make this a dict
        """
        mean = self.mean()
        rms = self.data.std()
        dmin = self.min()
        dmax = self.max()
        return (mean,rms,dmin,dmax)

    def smooth(self):
        #todo use specutils.manipulation.smoothing
        #option to smooth baseline too?
        pass


    @property
    def equivalencies(self):
        # @todo encapsulate take into account velframe and veldef
        """Get the spectral axis equivalencies that can be used in converting the axis between km/s and frequency or wavelength"""
        equiv = u.spectral()
        sa  = self.spectral_axis
        convention=self.velocity_convention
        if sa.doppler_rest is not None:
            rfq = sa.doppler_rest
        elif "RESTFREQ" in self.meta:
            cunit1 = self.meta.get("CUNIT1",self.wcs.wcs.cunit[0])
            #@todo this could be done with a dict str->function
            rfq = self.meta["RESTFREQ"]*cunit1
        else:
            rfq = None
        #print("RESTFREQ is ",rfq)
        if rfq is not None:
            if "radio" in self.velocity_convention:
                # Yeesh, the doppler_convention parameter for SpectralAxis.to does not match the doppler_convention list for Spectrum1D!
# This is actually bug in Spectrum1D documentation https://github.com/astropy/specutils/issues/1067
                convention="radio"
                equiv.extend(u.doppler_radio(rfq))
            elif "optical" in self.velocity_convention:
                convention="optical"
                equiv.extend(u.doppler_optical(rfq))
            elif "relativistic" in self.velocity_convention:
                convention="relativistic"
                equiv.extend(u.doppler_relativistic(rfq))
            elif "redshift" in self.velocity_convention:
                convention="redshift"
                equiv.extend(u.doppler_redshift())
        return equiv
    
