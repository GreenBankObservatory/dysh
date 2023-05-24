
from specutils import Spectrum1D, SpectrumList,SpectralRegion
from astropy.modeling.polynomial import Polynomial1D,Chebyshev1D
from astropy.modeling.fitting import LevMarLSQFitter,LinearLSQFitter
from astropy.wcs import WCS
from . import baseline

class Spectrum(Spectrum1D):
    """basic spectrum builts on Spectrum1D added stuff like baseline model"""
    def __init__(self, *args,**kwargs):
        Spectrum1D.__init__(self,*args,**kwargs)
        self._baseline_model = None

    def baseline(self,order,exclude=None,**kwargs):
        """compute and optionally remove a baseline"""
        kwargs_opts = {
            'remove':False,
            'model':'polynomial',
            'fitter':  LinearLSQFitter(calc_uncertainties=True),
        }
        kwargs_opts.update(kwargs)

        self._baseline_model = baseline(self,order,exclude,**kwargs)
        if kwargs_opts['remove']:
            s=self.subtract(self._baseline_model(self.spectral_axis))
            self._data = s._data

    def _undo_baseline(self):
        s=self.add(self._baseline_model(self.spectral_axis))
        self._data = s._data
        self._baseline_model = None
    
    def bshow(self):
        print(f"baseline model {self._baseline_model}")


    def stats(self):
        """ Return some stats. Note this works with slicing!"""
        mean = self.mean()
        rms = self.data.std()
        dmin = self.min()
        dmax = self.max()
        return (mean,rms,dmin,dmax)

    def smooth(self):
        #todo use specutils.manipulation.smoothing
        #option to smooth baseline too?
        pass

