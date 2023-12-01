import pathlib

import numpy as np
import pytest
from astropy.io import fits

import dysh
from dysh.fits import gbtfitsload


class TestPSScan:
    def test_tsys(self, data_dir):
        """
        Test that `getps` results in the same system temperature as GBTIDL.
        """
        sdf_file = f"{data_dir}/TGBT21A_501_11/TGBT21A_501_11.raw.vegas.fits"
        gbtidl_file = f"{data_dir}/TGBT21A_501_11/TGBT21A_501_11_getps_scan_152_intnum_0_ifnum_0_plnum_0.fits"

        sdf = gbtfitsload.GBTFITSLoad(sdf_file)
        tps = sdf.getps(152)
        tsys = tps[0].tsys

        hdu = fits.open(gbtidl_file)
        gbtidl_table = hdu[1].data
        gbtidl_tsys = gbtidl_table["TSYS"]

        assert (tsys - gbtidl_tsys)[0] == 0.0

    def test_compare_with_GBTIDL(self, data_dir):
        """ """
        # get filenames
        sdf_file = f"{data_dir}/TGBT21A_501_11/TGBT21A_501_11_ifnum_0_int_0-2.fits"
        gbtidl_file = f"{data_dir}/TGBT21A_501_11/TGBT21A_501_11_ifnum_0_int_0-2_getps_152_plnum_0.fits"

        # Generate the dysh result.
        sdf = gbtfitsload.GBTFITSLoad(sdf_file)
        # psscan is a ScanList.
        psscan = sdf.getps(152, plnum=0)
        assert len(psscan) == 1
        psscan.calibrate()
        # psscan_tavg is a list.
        psscan_tavg = psscan.timeaverage(weights="tsys")
        assert len(psscan_tavg) == 1

        # Load the GBTIDL result.
        hdu = fits.open(gbtidl_file)
        table = hdu[1].data
        psscan_gbtidl = table["DATA"][0]

        # Compare.
        diff = psscan_tavg[0].flux.value - psscan_gbtidl
        assert abs(np.nanmedian(diff)) < 1e-9
        assert psscan_tavg[0].meta["EXPOSURE"] == table["EXPOSURE"][0]
        assert psscan_tavg[0].meta["TSYS"] == pytest.approx(table["TSYS"][0], 1e-14)

    def test_compare_with_GBTIDL_2(self, data_dir):
        """
        Test `getps` when working with multiple scans.
        """

        data_path = f"{data_dir}/TGBT21A_501_11/NGC2782"
        sdf_file = f"{data_path}/TGBT21A_501_11_NGC2782.raw.vegas.A.fits"
        gbtidl_file = f"{data_path}/TGBT21A_501_11_getps_scans_156-158_ifnum_0_plnum_0_timeaverage.fits"

        sdf = gbtfitsload.GBTFITSLoad(sdf_file)
        ps_scans = sdf.getps([156, 158], plnum=0)
        ta = ps_scans.timeaverage()

        hdu = fits.open(gbtidl_file)
        table = hdu[1].data

        assert ta[0].meta["TSYS"] == pytest.approx(table["TSYS"], 1e-14)
        assert ta[0].meta["EXPOSURE"] == table["EXPOSURE"]
        assert np.all(np.abs(table["DATA"][0] - ta[0].flux.value) < 3e-7)

    @pytest.mark.skip(reason="We need to update this to work with multifits and ScanBlocks")
    def test_baseline_removal(self, data_dir):
        # get filenames
        data_path = f"{data_dir}/AGBT17A_404_01"
        sdf_file = f"{data_path}/AGBT17A_404_01_scan_19_prebaseline.fits"
        gbtidl_file = f"{data_path}/AGBT17A_404_01_scan_19_postbaseline.fits"
        gbtidl_modelfile = f"{data_path}/AGBT17A_404_01_scan_19_bline_model.fits"

        # generate the dysh result
        sdf = gbtfitsload.GBTFITSLoad(sdf_file)
        # has already been calibrated and smoothed slightly in gbtidl
        psscan = sdf.getspec(0)
        # 3rd order baseline fitted in gbtidl with regions [100,400] and [450,750]
        psscan.baseline(3, [(0, 100), (400, 450), (750, 819)], remove=True, model="chebyshev")

        # load gbtidl result
        # assuming the pre-baseline spectrum is still the same (can still test this if need be)
        hdu = fits.open(gbtidl_file)
        gbtidl_post = hdu[1].data["DATA"][0]
        hdu = fits.open(gbtidl_modelfile)
        gbtidl_bline_model = hdu[1].data["DATA"][0]

        diff = psscan - gbtidl_post

        # check that the spectra are the same but this won't pass right now
        # what is the tolerance for not passing?
        # [TODO] Find the right threshold for this
        # assert np.nanmean(diff) == 0

    def test_blank_integrations(self, data_dir):
        """
        Test `getps` when there are blanked integrations.
        """
        data_path = f"{data_dir}/TGBT21A_501_11/NGC2782_blanks"
        sdf_file = f"{data_path}/NGC2782.raw.vegas.A.fits"

        sdf = gbtfitsload.GBTFITSLoad(sdf_file)
        ps_sb = sdf.getps([156], plnum=0)
        ta1 = ps_sb.timeaverage()
        # This should not raise any errors.
        ta2 = ps_sb.timeaverage(None)
        # Check if the time average is all NaNs.
        all_nan = np.isnan(ta1[0].flux.value).sum() == len(ta1[0].flux)
        assert ~all_nan
        # Check that the metadata is accurate.
        # The system temperature is different because of the squared averaging.
        assert abs(ps_sb[0].calibrated(0).meta["TSYS"] - ta1[0].meta["TSYS"]) < 5e-16
        assert (ps_sb[0].calibrated(0).meta["EXPOSURE"] - ta1[0].meta["EXPOSURE"]) == 0.0
        # Check if the time averaged data matches that from the first integration.
        # assert np.all(abs(ps_sb[0].calibrated(0).flux.value - ta1[0].flux.value) < 2e-19)
        # Set to 5E-16 because Windows OS tests fail below that.  Need to understand why.
        assert np.all(abs(ps_sb[0].calibrated(0).flux.value - ta1[0].flux.value) < 5e-16)


class TestSubBeamNod:
    def test_compare_with_GBTIDL(self, data_dir):
        # get filenames
        # We still need a data file with a single scan in it
        sdf_file = f"{data_dir}/TRCO_230413_Ka/TRCO_230413_Ka_scan43.fits"
        gbtidl_file = f"{data_dir}/TRCO_230413_Ka/TRCO_230413_Ka_snodka_43_ifnum_0_plnum_0_fdnum_1.fits"

        # Generate the dysh result.
        # snodka-style. Need test for method='cycle'
        sdf = gbtfitsload.GBTFITSLoad(sdf_file)
        sbn = sdf.subbeamnod(43, sig=None, cal=None, ifnum=0, fdnum=1, calibrate=True, weights="tsys", method="scan")

        # Load the GBTIDL result.
        hdu = fits.open(gbtidl_file)
        nodscan_gbtidl = hdu[1].data["DATA"][0]

        # Compare.
        ratio = sbn[0].calibrated(0).flux.value / nodscan_gbtidl
        # kluge for now since there is a small wavy pattern in
        # the difference at the ~0.06 K level
        assert np.nanmedian(ratio) <= 0.998


class TestTPScan:
    def test_tsys(self, data_dir):
        """
        Test that `gettp` produces the same system temperature as GBTDIL.
        """

        sdf_file = f"{data_dir}/TGBT21A_501_11/TGBT21A_501_11.raw.vegas.fits"
        gbtidl_file = f"{data_dir}/TGBT21A_501_11/TGBT21A_501_11_getps_scan_152_intnum_0_ifnum_0_plnum_0.fits"

        sdf = gbtfitsload.GBTFITSLoad(sdf_file)
        tps = sdf.gettp(153)
        tsys = tps[0].tsys

        hdu = fits.open(gbtidl_file)
        gbtidl_table = hdu[1].data
        gbtidl_tsys = gbtidl_table["TSYS"]

        assert (tsys - gbtidl_tsys)[0] == 0.0

    def test_exposure(self, data_dir):
        """
        Test that `gettp` holds the same exposure times as GBTIDL.
        """

        sdf_file = f"{data_dir}/TGBT21A_501_11/TGBT21A_501_11_scan_152_ifnum_0_plnum_0.fits"
        gbtidl_file = f"{data_dir}/TGBT21A_501_11/TGBT21A_501_11_gettp_scan_152_ifnum_0_plnum_0_keepints.fits"

        sdf = gbtfitsload.GBTFITSLoad(sdf_file)
        tp = sdf.gettp(152)

        hdu = fits.open(gbtidl_file)
        table = hdu[1].data

        for i in range(len(tp[0]._scanrows) // 2):
            assert tp[0].total_power(i).meta["EXPOSURE"] == table["EXPOSURE"][i]

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
        assert np.nanmean((tpavg[0].flux.value - data[-1]) / data[-1].mean()) < 2**-32

    def test_compare_with_GBTIDL_equal_weights(self, data_dir):
        """
        This test compares `gettp` when using equal weights.
        It takes a scan with multiple integrations and averages
        them using unity weights.
        It checks that the exposure and averaged spectra are the same,
        and that the system temperature is the same up to the precision
        used by GBTIDL.
        """
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
        assert np.all((data[0] - tpavg[0].flux.value.astype(np.float32)) == 0.0)
