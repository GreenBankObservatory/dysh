import astropy.units as u
from .spectrum import Spectrum

class PSScan(object):
    """
    Holds a position scan pair
    """
    def __init__(self, sdfits, scan_on, scan_off):
        """
        """
        self.sdfits = sdfits # parent class
        self.status = 0
        #                           # ex1:
        self.nint = 0               # 11
        self.npol = 0               #  2
        self.on = None              # 44
        self.off = None             # 44
        self.calibrated = None      # 22
        self.timeaveraged = None    #  2
        self.polaveraged = None     #  1
        #
        self.nrows = len(scan_on)
        
    def __len__(self):
        return self.nrows
    
class GBTPSScan(PSScan):
    def __init__(self, gbtfits, scans,scanrows,bintable=0):
        PSScan.__init__(self,gbtfits,scans['ON'],scans['OFF'])
        # The rows of the original bintable corresponding to ON and OFF
        self._scanrows = scanrows
        self.npol =  gbtfits.npol(bintable) #TODO deal with bintable
        self.nint = gbtfits.nintegrations(bintable)
        self.on = np.empty(self.nrows, dtype=Spectrum)
        self.off = np.empty(self.nrows, dtype=Spectrum)
        # todo use gbtfits.velocity_convention(veldef,velframe)
        vc = "doppler_radio"
        for i in self._scanrows["ON"]:
            data = np.array([[gbtfits.rawspectrum(bintable,i)]])
            meta = dict(gbtfits.index(bintable).iloc[i])
            #@TODO allow WCS? big performance hit
            self.on[i] = Spectrum(data*u.K,wcs=None,meta=meta,velocity_convention=vc)
        for i in self._scanrows["OFF"]:
            data = np.array([[gbtfits.rawspectrum(bintable,i)]])
            meta = dict(gbtfits.index(bintable).iloc[i])
            #@TODO allow WCS? big performance hit
            self.off[i] = Spectrum(data*u.K,wcs=None,meta=meta,velocity_convention=vc)
