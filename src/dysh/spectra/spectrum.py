from specutils import Spectrum1D, SpectrumList,SpectralRegion
from astropy.modeling.polynomial import Polynomial1D,Chebyshev1D
from astropy.modeling.fitting import LevMarLSQFitter,LinearLSQFitter
from astropy.wcs import WCS
from . import baseline

class Spectrum(Spectrum1D):
    """This class contains a spectrum and its attributes. It is built on
    `~specutils.Spectrum1D` with added attributes like baseline model.
    Note that `~specutils.Spectrum1D` can contain multiple spectra but
    we probably will not use that because the restriction that it can
    have only one spectral axis conflicts with slight Doppler shifts.
    See `~specutils.Spectrum1D` for the instantiation arguments.
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
            exclude: list of 2-tuples
                List of channel-based regions to exclude in the fitting in form [lower,upper]. Default: None
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

