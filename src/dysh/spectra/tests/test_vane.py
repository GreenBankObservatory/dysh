
import numpy as np
import pytest
from astropy import units as u
from astropy.io import fits

import dysh.util as util
from dysh.fits import gbtfitsload, sdfitsload
from dysh.util.files import dysh_data

from dysh.fits.core import mean_data
from dysh.fits.core import getbeam
from dysh.fits.core import calseq
from dysh.fits.core import vanecal
from dysh.fits.core import plot_vegas      # different location ?
from dysh.fits.core import getnod          # should be done by official getnod

class TestVaneScan:
    def test_vane1(self, data_dir):
        """
        test1
        """
        rtol = 1e-4      # do tests on this relative tolerance
        
        sdf_file = dysh_data(test='AGBT21B_024_14/AGBT21B_024_14_test')
        sdf2 = gbtfitsload.GBTFITSLoad(sdf_file)

        beam2 = getbeam(sdf2)
        assert len(beam2) == 2
        assert beam2[0] == 1 and beam2[1] == 9

        #  [1, 9]
        
        tcal = 272
        tsys2 = vanecal(sdf2, [329, 330], feeds=beam2, tcal=tcal)
        #  221.79946241 201.2739606

        sp1,sp2 = getnod(sdf2, [331, 332], beam2, tsys=tsys2)
        sp3,sp4 = getnod(sdf2, [333, 334], beam2, tsys=tsys2)
        sp5 = sp1.average([sp2,sp3,sp4])
        tint = sp5.meta["EXPOSURE"]
        assert np.isclose(tint, 3.942, rtol=rtol)    # 3.941774845123291
        
        tsys = sp5.meta["TSYS"] 
        assert np.isclose(tsys, 210.54, rtol=rtol)  # 210.54325052318202

        #s4 = [float(s) for s in sp5.stats(qac=True).split()]
        rms = sp5.stats()['rms'].value           #  0.2357084722315336
        print(rms)
        assert np.isclose(rms, 0.47509, rtol=rtol)
        
