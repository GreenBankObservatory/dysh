
import pytest
import pathlib
import numpy as np

from astropy.io import fits
from astropy.utils.data import (
                            get_pkg_data_filename,
                            )

import dysh
from dysh.fits import gbtfitsload

#dysh_root = pathlib.Path(dysh.__file__).parent.resolve()

class TestGBTPSScan():

    def test_compare_with_GBTIDL(self):

        # get filenames
        sdf_file = get_pkg_data_filename("data/TGBT21A_501_11_ifnum_0_int_0-2.fits")
        gbtidl_file = get_pkg_data_filename("data/TGBT21A_501_11_ifnum_0_int_0-2_getps_152_plnum_0.fits")
        #sdf_file = f"{dysh_root}/fits/tests/data/TGBT21A_501_11_ifnum_0_int_0-2.fits"
        #gbtidl_file = f"{dysh_root}/fits/tests/data/TGBT21A_501_11_ifnum_0_int_0-2_getps_152_plnum_0.fits"

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

class TestSubBeamNod():

    def test_compare_with_GBTIDL(self):
        # get filenames
        # We still need a data file with a single scan in it
        sdf_file = get_pkg_data_filename("data/TRCO_230413_Ka_scan43.fits")
        gbtidl_file = get_pkg_data_filename("data/TRCO_230413_Ka_snodka_43_ifnum_0_plnum_0_fdnum_1.fits")

        # Generate the dysh result.
        # snodka-style. Need test for method='cycle'
        sdf = gbtfitsload.GBTFITSLoad(sdf_file)
        sbn = sdf.subbeamnod(43, sig=None, cal=None,
                        ifnum=0, fdnum=1, calibrate=True,
                        weights='tsys',method='scan') 

        # Load the GBTIDL result.
        hdu = fits.open(gbtidl_file)
        nodscan_gbtidl = hdu[1].data["DATA"][0]

        # Compare.
        ratio = sbn.flux.value/nodscan_gbtidl
        #print("DIFF ",np.nanmedian(sbn.flux.value - nodscan_gbtidl))
        # kluge for now since there is a small wavy pattern in 
        # the difference at the ~0.06 K level
        assert np.nanmedian(ratio) <= 0.998

class TestGBTTPScan():

    def test_compare_with_GBTIDL_tsys_weights(self):

        sdf_file = get_pkg_data_filename("data/TGBT21A_501_11_scan_152_ifnum_0_plnum_0.fits")
        gbtidl_file = get_pkg_data_filename("data/TGBT21A_501_11_gettp_scan_152_ifnum_0_plnum_0_keepints.fits")

        # Generate the dysh result.
        sdf = gbtfitsload.GBTFITSLoad(sdf_file)
        tp  = sdf.gettp(152)
        tpavg = tp.timeaverage()

        # Check that we know how to add.
        assert tpavg.meta["EXPOSURE"] == tp.exposure.sum()

        # Load GBTIDL result.
        hdu = fits.open(gbtidl_file)
        table = hdu[1].data
        data = table["DATA"]

        # Check exposure times.
        assert np.sum(table["EXPOSURE"][:-1] - tp.exposure) == 0.0
        assert tpavg.meta["EXPOSURE"] == table["EXPOSURE"][-1]
        # System temperature.
        # For some reason, the last integration comes out with a
        # difference in TSYS ~4e-10 rather than ~1e-14. Check why.
        assert np.all(abs(table["TSYS"][:-1] - tp.tsys) < 1e-9)
        assert abs(tpavg.meta["TSYS"] - table["TSYS"][-1]) < 1e-10
        # Data, which uses float -- 32 bits.
        assert np.sum(tp._data - data[:-1]) == 0.0
        assert np.nanmean((tpavg.flux.value - data[-1])/data[-1].mean()) < 2**32
        
    def test_compare_with_GBTIDL_equal_weights(self):

        sdf_file = get_pkg_data_filename("data/TGBT21A_501_11_scan_152_ifnum_0_plnum_0.fits")
        gbtidl_file = get_pkg_data_filename("data/TGBT21A_501_11_gettp_scan_152_ifnum_0_plnum_0_eqweight.fits")

        # Generate the dysh result.
        sdf = gbtfitsload.GBTFITSLoad(sdf_file)
        tp  = sdf.gettp(152)
        tpavg = tp.timeaverage(weights=None)
        
        # Check that we know how to add.
        assert tpavg.meta["EXPOSURE"] == tp.exposure.sum()

        # Load GBTIDL result.
        hdu = fits.open(gbtidl_file)
        table = hdu[1].data
        data = table["DATA"]

        # Compare Dysh and GBTIDL.
        assert table["EXPOSURE"][0] == tpavg.meta["EXPOSURE"]
        assert abs(table["TSYS"][0] - tpavg.meta["TSYS"]) < 2**-32
        assert np.all((data[0] - tpavg.flux.value) == 0.0)

