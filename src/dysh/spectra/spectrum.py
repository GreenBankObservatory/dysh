
from specutils import Spectrum1D, SpectrumList,SpectralRegion

class Spectrum(Spectrum1D):
    """basic spectrum builts on Spectrum1D added stuff like baseline model"""
    def __init__(self, *args,**kwargs):
        Spectrum1D.__init__(self,*args,**kwargs)
        self._baseline_model = None

    def baseline(self,order,exclude=None,**kwargs):
        """compute and optionally remove a baseline"""
        kwargs_opts = {
            'remove':True,
            'model':'polynomial',
            'fitter':  LinearLSQFitter(calc_uncertainties=True),
        }
        kwargs_opts.update(kwargs)

        self._baseline_model = baseline(self,order,exclude,**kwargs)
        if kwargs_opts['remove']:
            self.subtract(self._baseline_model.evaluate(self.spectral_axis))

    def __undo_baseline(self):
        self.add(self._baseline_model.evaluate(self.spectral_axis))
        self._baseline_model = None
    
    def bshow(self);
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

