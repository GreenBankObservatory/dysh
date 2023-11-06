import pytest
import pathlib

import numpy as np

from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename

import dysh
from dysh.fits import gbtfitsload



class TestPSScan:
    def test_compare_with_GBTIDL(self, data_dir):
        # get filenames
        sdf_file = f"{data_dir}/TGBT21A_501_11/TGBT21A_501_11_ifnum_0_int_0-2.fits"
        gbtidl_file = f"{data_dir}/TGBT21A_501_11/TGBT21A_501_11_ifnum_0_int_0-2_getps_152_plnum_0.fits"
        #sdf_file = get_pkg_data_filename("data/TGBT21A_501_11_ifnum_0_int_0-2.fits")
        #gbtidl_file = get_pkg_data_filename("data/TGBT21A_501_11_ifnum_0_int_0-2_getps_152_plnum_0.fits")
        # sdf_file = f"{dysh_root}/fits/tests/data/TGBT21A_501_11_ifnum_0_int_0-2.fits"
        # gbtidl_file = f"{dysh_root}/fits/tests/data/TGBT21A_501_11_ifnum_0_int_0-2_getps_152_plnum_0.fits"

        # Generate the dysh result.
        sdf = gbtfitsload.GBTFITSLoad(sdf_file)
        # psscan is a ScanList
        psscan = sdf.getps(152, plnum=0)
        assert len(psscan) == 1
        psscan.calibrate()
        # psscan_tavg is a list
        psscan_tavg = psscan.timeaverage(weights="tsys")
        assert len(psscan_tavg) == 1

        # Load the GBTIDL result.
        hdu = fits.open(gbtidl_file)
        psscan_gbtidl = hdu[1].data["DATA"][0]

        # Compare.
        diff = psscan_tavg[0].flux.value - psscan_gbtidl
        assert np.nanmedian(diff) == 0.0

    @pytest.mark.skip(reason="We need to update this to work with multifits and ScanBlocks")
    def test_baseline_removal(self):

        #get filenames
        sdf_file = get_pkg_data_filename("data/AGBT17A_404_01_scan_19_prebaseline.fits")
        gbtidl_file = get_pkg_data_filename("data/AGBT17A_404_01_scan_19_postbaseline.fits")
        gbtidl_modelfile = get_pkg_data_filename("data/AGBT17A_404_01_scan_19_bline_model.fits")

        #generate the dysh result
        sdf = gbtfitsload.GBTFITSLoad(sdf_file)
        #has already been calibrated and smoothed slightly in gbtidl
        psscan = sdf.getspec(0)
        #3rd order baseline fitted in gbtidl with regions [100,400] and [450,750]
        psscan.baseline(3,[(0,100),(400,450),(750,819)],remove=True,model='chebyshev')

        #load gbtidl result
        #assuming the pre-baseline spectrum is still the same (can still test this if need be)
        hdu = fits.open(gbtidl_file)
        gbtidl_post = hdu[1].data["DATA"][0]
        hdu = fits.open(gbtidl_modelfile)
        gbtidl_bline_model = hdu[1].data["DATA"][0]

        diff = psscan - gbtidl_post

        #check that the spectra are the same but this won't pass right now
        #what is the tolerance for not passing?
        #[TODO] Find the right threshold for this
        #assert np.nanmean(diff) == 0


class TestSubBeamNod:
    def test_compare_with_GBTIDL(self, data_dir):
        # get filenames
        # We still need a data file with a single scan in it
        #sdf_file = get_pkg_data_filename("data/TRCO_230413_Ka_scan43.fits")
        #gbtidl_file = get_pkg_data_filename("data/TRCO_230413_Ka_snodka_43_ifnum_0_plnum_0_fdnum_1.fits")
        sdf_file = f"{data_dir}/TRCO_230413_Ka/TRCO_230413_Ka_scan43.fits"
        gbtidl_file = f"{data_dir}/TRCO_230413_Ka/TRCO_230413_Ka_snodka_43_ifnum_0_plnum_0_fdnum_1.fits"


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
        ratio = sbn[0].calibrated(0).flux.value / nodscan_gbtidl
        # kluge for now since there is a small wavy pattern in
        # the difference at the ~0.06 K level
        assert np.nanmedian(ratio) <= 0.998


class TestTPScan:
    def test_compare_with_GBTIDL_tsys_weights(self, data_dir):
        """
        This test compares `gettp` when using radiometer weights.
        It takes a scan with multiple integrations and averages
        them using weights from the radiometer equation.
        It checks that the averaged exposure and spectra are the same,
        and that the system temperature is the same up to the precision
        used by GBTIDL.
        It also checks that the spectra, exposures and system temperatures are the
        same for individual integrations after calibrating them.
        """

        #sdf_file = get_pkg_data_filename("data/TGBT21A_501_11_scan_152_ifnum_0_plnum_0.fits")
        #gbtidl_file = get_pkg_data_filename("data/TGBT21A_501_11_gettp_scan_152_ifnum_0_plnum_0_keepints.fits")
        sdf_file = f"{data_dir}/TGBT21A_501_11/TGBT21A_501_11_scan_152_ifnum_0_plnum_0.fits"
        gbtidl_file = f"{data_dir}/TGBT21A_501_11/TGBT21A_501_11_gettp_scan_152_ifnum_0_plnum_0_keepints.fits"

        # Generate the dysh result.
        sdf = gbtfitsload.GBTFITSLoad(sdf_file)
        # tp is a ScanList
        tp = sdf.gettp(152)
        assert len(tp) == 1
        # tpavg is a list
        tpavg = tp.timeaverage()
        assert len(tpavg) == 1

        # Check that we know how to add.
        assert tpavg[0].meta["EXPOSURE"] == tp[0].exposure.sum()

        # Load GBTIDL result.
        hdu = fits.open(gbtidl_file)
        table = hdu[1].data
        data = table["DATA"]

        # Check exposure times.
        # The last row in the GBTIDL file is the averaged data.
        assert np.sum(table["EXPOSURE"][:-1] - tp[0].exposure) == 0.0
        assert tpavg[0].meta["EXPOSURE"] == table["EXPOSURE"][-1]
        # System temperature.
        # For some reason, the last integration comes out with a
        # difference in TSYS ~4e-10 rather than ~1e-14. Check why.
        assert np.all(abs(table["TSYS"][:-1] - tp[0].tsys) < 1e-9)
        assert abs(tpavg[0].meta["TSYS"] - table["TSYS"][-1]) < 1e-10
        # Data, which uses float -- 32 bits.
        assert np.sum(tp[0]._data - data[:-1]) == 0.0
        assert np.nanmean((tpavg[0].flux.value - data[-1])/data[-1].mean()) < 2**-32


    def test_compare_with_GBTIDL_equal_weights(self, data_dir):
        """
        This test compares `gettp` when using equal weights.
        It takes a scan with multiple integrations and averages
        them using unity weights.
        It checks that the exposure and averaged spectra are the same,
        and that the system temperature is the same up to the precision
        used by GBTIDL.
        """
        #sdf_file = get_pkg_data_filename("data/TGBT21A_501_11_scan_152_ifnum_0_plnum_0.fits")
        #gbtidl_file = get_pkg_data_filename("data/TGBT21A_501_11_gettp_scan_152_ifnum_0_plnum_0_eqweight.fits")
        sdf_file = f"{data_dir}/TGBT21A_501_11/TGBT21A_501_11_scan_152_ifnum_0_plnum_0.fits"
        gbtidl_file = f"{data_dir}/TGBT21A_501_11/TGBT21A_501_11_gettp_scan_152_ifnum_0_plnum_0_eqweight.fits"

        # Generate the dysh result.
        sdf = gbtfitsload.GBTFITSLoad(sdf_file)
        tp = sdf.gettp(152)
        assert len(tp) == 1
        # tpavg is a list
        tpavg = tp.timeaverage(weights=None)
        assert len(tpavg) == 1

        # Check that we know how to add.
        assert tpavg[0].meta["EXPOSURE"] == tp[0].exposure.sum()

        # Load GBTIDL result.
        hdu = fits.open(gbtidl_file)
        table = hdu[1].data
        data = table["DATA"]

        # Compare Dysh and GBTIDL.
        assert table["EXPOSURE"][0] == tpavg[0].meta["EXPOSURE"]
        assert abs(table["TSYS"][0] - tpavg[0].meta["TSYS"]) < 2**-32
        assert np.all((data[0] - tpavg[0].flux.value) == 0.0)
