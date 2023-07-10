
import pytest
import numpy as np

from astropy.io import fits

from dysh.fits import gbtfitsload


class TestGBTPSScan():

    def test_compare_with_GBTIDL(self):

        sdf_file = "../../fits/tests/data/TGBT21A_501_11_ifnum_0_int_0-2.fits"
        gbtidl_file = "../../fits/tests/data/TGBT21A_501_11_ifnum_0_int_0-2_getps_152_plnum_0.fits"

        # Generate the dysh result.
        sdf = gbtfitsload.GBTFITSLoad(sdf_file)
        psscan = sdf.getps(152, plnum=0)
        psscan.calibrate()
        psscan_tavg = psscan.timeaverage(weights="tsys")

        # Load the GBTIDL result.
        hdu = fits.open(gbtidl_file)
        psscan_gbtidl = hdu[1].data["DATA"][0]

        # Compare.
        diff = psscan_tavg.flux.value - psscan_gbtidl
        assert np.nanmedian(diff) == 0.0
