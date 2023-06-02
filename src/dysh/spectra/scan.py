import astropy.units as u
from copy import deepcopy
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
    def __init__(self, gbtfits, scans, scanrows, bintable=0):
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

    def calibrate(self, **kwargs):
        """
        special PS calibration
        There are some arguments how *exactly* this is done
        """
        kwargs_opts = {
            'verbose': False
        }
        kwargs_opts.update(kwargs)

        self.status = 1
        npolint = self.npol * self.nint
        self.calibrated = np.empty(npolint, dtype=Spectrum)
        for i in range(npolint):
            tcal = self.off[2*i].meta['TCAL']
            tcal2= self.on[2*i].meta['TCAL']
            tsys = dcmeantsys(self.off[2*i].data,  self.off[2*i+1].data,tcal)
            tsys2= dcmeantsys(self.on[2*i].data,  self.on[2*i+1].data,tcal2)
            if kwargs_opts['verbose']: 
                print(i,tcal,tsys,tsys2)
            #                 2*i is the CalON     2*i+1 the CalOFF
            #
            sig = 0.5*(self.on[2*i].data + self.on[2*i+1].data)
            ref = 0.5*(self.off[2*i].data + self.off[2*i+1].data)
            #kr = self.on[2*i].gbt['row'] 
            #if True:
            self.calibrated[i] = deepcopy(self.on[2*i])
            #else:
            #    self.calibrated[i] = Spectrum(self.sdfits.data2[kr])
            self.calibrated[i].data = tsys * (sig-ref) / ref
            #self.calibrated[i].gbt['row'] = kr
            self.calibrated[i].meta['tsys'] = tsys
            # fix the meta data ; most of it is ok

