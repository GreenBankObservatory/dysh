import numpy as np
from specutils import Spectrum1D, SpectrumList
import astropy.units as u
from memory_profiler import profile
import sys

@profile
def testit(nchan,num):
    # Make no difference if use spectrumlist or []
    if False:
        sl = SpectrumList()
    else:
        sl = []
    #default size 
    rawsize=sys.getsizeof(np.random.rand(nchan))
    blocksize=rawsize*num
    print(f"Size of single {nchan} channel spectrum {rawsize/(1024*1024):.2} Mbytes")
    print(f"X {num} = {blocksize} bytes = {blocksize/(1048576*1024):.2} GB")
    for j in range(num):
        if j%500 == 0:
            print(f"doing row {j}")
        sp = np.random.rand(nchan)*u.K
        p=Spectrum1D(flux=sp, wcs=None,meta={},velocity_convention="doppler_radio")
        sl.append(p)

if __name__ == "__main__":
    testit(nchan=32768,num=6040)
