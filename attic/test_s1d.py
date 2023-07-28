import numpy as np
from specutils import Spectrum1D, SpectrumList, SpectralRegion
from specutils.fitting import fit_continuum
import astropy.units as u
from astropy.modeling.models import Gaussian1D
from astropy.modeling.polynomial import Polynomial1D
from memory_profiler import profile
import matplotlib.pyplot as plt
import sys
import resource
import gc
import objgraph

#@profile
def testit(nchan,num):
    # Make no difference if use spectrumlist or []
    if False:
        sl = SpectrumList()
    else:
        sl = []
    #default size 
    raw = np.random.rand(nchan)
    rawsize=sys.getsizeof(raw)
    blocksize=rawsize*num
    x = np.linspace(0,10,nchan)
    g1 = Gaussian1D(amplitude=10,mean=5,stddev=0.3)
    print(f"Size of single {nchan} channel spectrum {rawsize/(1024*1024):.2} Mbytes")
    print(f"X {num} = {blocksize} bytes = {blocksize/(1048576*1024):.2} GB")
    for j in range(num):
        if j%500 == 0:
            print(f"doing row {j}")
        sp = (20.+g1(x)+np.random.normal(size=nchan))*u.K
        p=Spectrum1D(flux=sp, spectral_axis=x*u.um,wcs=None,meta={},velocity_convention="doppler_radio")
        sl.append(p)
    return sl

def baseline(speclist,order,exclude=None,plot=False):
    fc = fit_continuum(p,Polynomial1D(degree=order),exclude_regions=[exclude])
    if plot:
        fig,ax = plt.subplots()
        ax.plot(x,p.flux)
        ax.plot(x,fc(x))
        plt.show()

if __name__ == "__main__":
    sl=testit(nchan=2048,num=2000)
    print(f'Memory usage: %s (Gb) {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(1024*1024)}')
    gc.collect()
    print(f'Post-GC Memory usage: %s (Gb) {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(1024*1024)}')
    objgraph.show_most_common_types()
    p = sl[0]
    x = sl[0].spectral_axis
    nchan = len(x)
    cdelt = x[1]-x[0]
    center = x[nchan // 2]
    width = 0.25*nchan*cdelt
    exclude = SpectralRegion.from_center(center,width)
    i=0
    for sp in sl:
        if (i%500) == 0:
            print("baseline ",i)
        baseline(sp,1,exclude)
        i=1+i
