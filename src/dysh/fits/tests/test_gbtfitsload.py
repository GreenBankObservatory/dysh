import glob
import os
import pathlib
import platform
import shutil
from copy import deepcopy
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from astropy.io import fits
from pandas.testing import assert_frame_equal, assert_series_equal

from dysh import util
from dysh.fits import gbtfitsload


class TestGBTFITSLoad:
    """ """

    def setup_method(self):
        self.root_dir = util.get_project_root()
        self.data_dir = f"{self.root_dir}/testdata"
        self._file_list = glob.glob(f"{self.data_dir}/TGBT21A_501_11/*.fits")

    def test_load(self):
        """
        Test loading 8 different sdfits files.
        Check: number of pandas rows loaded is equal to the expected number.
        """
        expected = {
            "TGBT21A_501_11.raw.vegas.fits": 4,
            "TGBT21A_501_11_2.raw.vegas.fits": 8,
            "TGBT21A_501_11_getps_scan_152_intnum_0_ifnum_0_plnum_0.fits": 1,
            "TGBT21A_501_11_gettp_scan_152_ifnum_0_plnum_0_cal_state_0.fits": 1,
            "TGBT21A_501_11_gettp_scan_152_ifnum_0_plnum_0_cal_state_1.fits": 1,
            "TGBT21A_501_11_ifnum_0_int_0-2.fits": 24,
            "TGBT21A_501_11_ifnum_0_int_0-2_getps_152_plnum_0.fits": 1,
            "TGBT21A_501_11_ifnum_0_int_0-2_getps_152_plnum_1.fits": 1,
            "TGBT21A_501_11_gettp_scan_152_ifnum_0_plnum_0.fits": 1,
            "TGBT21A_501_11_gettp_scan_152_ifnum_0_plnum_0_eqweight.fits": 1,
            "TGBT21A_501_11_gettp_scan_152_ifnum_0_plnum_0_keepints.fits": 152,
            "TGBT21A_501_11_scan_152_ifnum_0_plnum_0.fits": 302,
            "getps_154_ifnum_0_plnum_0_intnum_0.fits": 1,
            "TGBT21A_501_11.raw.156.fits": 7,
            "testselection.fits": 50,
            "TGBT21A_501_11_getps_scan_152_ifnum_0_plnum_0_smthoff_15.fits": 1,
            "TGBT21A_504_01/TGBT21A_504_01_gettp_scan_20_ifnum_0_plnum_1_sig_state_0_cal_state_1.fits": 1,
            "TGBT21A_504_01/TGBT21A_504_01_gettp_scan_20_ifnum_0_plnum_1_sig_state_0_cal_state_0.fits": 1,
            "TGBT21A_504_01/TGBT21A_504_01_gettp_scan_20_ifnum_0_plnum_1_sig_state_0_cal_all.fits": 1,
        }

        for fnm in self._file_list:
            print(fnm)

            filename = os.path.basename(fnm)
            sdf = gbtfitsload.GBTFITSLoad(fnm)
            assert len(sdf.index(bintable=0)) == expected[filename]

    def test_getspec(self):
        """
        Test that a GBTFITSLoad object can use the `getspec` function.
        """

        sdf_file = f"{self.data_dir}/TGBT21A_501_11/TGBT21A_501_11.raw.vegas.fits"
        sdf = gbtfitsload.GBTFITSLoad(sdf_file)
        sdf.flag_channel([[0, 100]])
        sdf.apply_flags()
        spec = sdf.getspec(0, setmask=False)
        spec2 = sdf.getspec(0, setmask=True)
        assert any(spec.mask != spec2.mask)
        assert all(spec2.mask[0:101])
        assert all(spec2.mask[102:] == False)

    def test_getps_single_int(self):
        """
        Compare gbtidl result to dysh for a getps spectrum from a single integration/pol/feed.
        For the differenced spectrum (gbtidl - dysh) we check:
         - median is 0.0
         - all non-nan values are less than 5e-7
         - index 3072 is nan (why?)
        """
        gbtidl_file = f"{self.data_dir}/TGBT21A_501_11/TGBT21A_501_11_getps_scan_152_intnum_0_ifnum_0_plnum_0.fits"
        # We should probably use dysh to open the file...
        hdu = fits.open(gbtidl_file)
        gbtidl_getps = hdu[1].data["DATA"][0]

        sdf_file = f"{self.data_dir}/TGBT21A_501_11/TGBT21A_501_11.raw.vegas.fits"
        sdf = gbtfitsload.GBTFITSLoad(sdf_file)
        psscan = sdf.getps(scan=152)
        assert len(psscan) == 1
        psscan.calibrate()
        dysh_getps = psscan[0].calibrated(0).flux.to("K").value

        diff = gbtidl_getps - dysh_getps
        hdu.close()
        assert np.nanmedian(diff) == 0.0
        assert np.all(abs(diff[~np.isnan(diff)]) < 5e-7)
        assert np.isnan(diff[3072])

        # add a test to ensure that a bad scan number raises an ValueError (issue 462)
        with pytest.raises(ValueError):
            sdf.getps(scan=99999)

    def test_getps_acs(self):
        """
        Compare `GBTIDL` result to `dysh` with ACS data.
        """

        data_dir = util.get_project_testdata() / "AGBT05B_047_01"
        sdf_file = data_dir / "AGBT05B_047_01.raw.acs"
        idl_file = data_dir / "gbtidl" / "AGBT05B_047_01.getps.acs.fits"

        sdf = gbtfitsload.GBTFITSLoad(sdf_file)
        getps = sdf.getps(scan=51, ifnum=0, plnum=0)
        ps = getps.timeaverage()
        ps_vals = ps.flux.value

        hdu = fits.open(idl_file)
        table = hdu[1].data
        gbtidl_spec = table["DATA"][0]
        hdu.close()

        # Do not compare NaN values.
        mask = np.isnan(ps_vals) | np.isnan(gbtidl_spec)

        # Compare data.
        diff = ps_vals[~mask] - gbtidl_spec[~mask]
        # try:
        assert np.all(diff < 1e-3)
        # except AssertionError:
        #    print(f"Comparison with GBTIDL ACS Spectrum failed, mean difference is {np.nanmean(diff)}")

        for col in table.names:
            if col not in ["DATA"]:
                try:
                    ps.meta[col]
                except KeyError:
                    continue
                try:
                    assert ps.meta[col] == pytest.approx(table[col][0], 1e-3)
                except AssertionError:
                    print(f"{col} fails: {ps.meta[col]}, {table[col][0]}")

    def test_gettp(self):
        """
        Compare gbtidl result to dysh for a gettp spectrum from a single polarization and feed and
        different cal states.
        For the differenced spectrum (gbtidl - dysh) we check:
        For the noise calibration diode on, off, and both:
         - mean value is 0.0
        """
        # @todo refactor the repeated gbtidl/tp0 sections here.
        # Get the answer from GBTIDL.
        gbtidl_file = f"{self.data_dir}/TGBT21A_501_11/TGBT21A_501_11_gettp_scan_152_ifnum_0_plnum_0_cal_state_1.fits"
        hdu = fits.open(gbtidl_file)
        gbtidl_gettp = hdu[1].data["DATA"][0]
        gbtidl_exp = hdu[1].data["EXPOSURE"][0]
        gbtidl_tsys = hdu[1].data["TSYS"][0]
        hdu.close()

        # Get the answer from dysh.
        sdf_file = f"{self.data_dir}/TGBT21A_501_11/TGBT21A_501_11.raw.vegas.fits"
        sdf = gbtfitsload.GBTFITSLoad(sdf_file)
        tps_on = sdf.gettp(scan=152, sig=True, cal=True, calibrate=True, ifnum=0, plnum=0)
        assert len(tps_on) == 1

        # Compare.
        tp0 = tps_on[0].total_power(0)
        diff = tp0.flux.value - gbtidl_gettp
        assert np.nanmean(diff) == 0.0
        assert tp0.meta["TSYS"] == pytest.approx(gbtidl_tsys)
        assert tp0.meta["EXPOSURE"] == pytest.approx(gbtidl_exp)

        # Now with the noise diode Off.
        tps_off = sdf.gettp(scan=152, sig=True, cal=False, calibrate=True, ifnum=0, plnum=0)
        assert len(tps_off) == 1
        gbtidl_file = f"{self.data_dir}/TGBT21A_501_11/TGBT21A_501_11_gettp_scan_152_ifnum_0_plnum_0_cal_state_0.fits"
        hdu = fits.open(gbtidl_file)
        gbtidl_gettp = hdu[1].data["DATA"][0]
        gbtidl_exp = hdu[1].data["EXPOSURE"][0]
        gbtidl_tsys = hdu[1].data["TSYS"][0]
        tp0 = tps_off[0].total_power(0)
        diff = tp0.flux.value - gbtidl_gettp
        hdu.close()
        assert np.nanmean(diff) == 0.0
        assert tp0.meta["TSYS"] == pytest.approx(gbtidl_tsys)
        assert tp0.meta["EXPOSURE"] == pytest.approx(gbtidl_exp)

        # Now, both on and off.
        tps = sdf.gettp(scan=152, sig=None, cal=None, calibrate=True, ifnum=0, plnum=0)
        assert len(tps) == 1
        gbtidl_file = f"{self.data_dir}/TGBT21A_501_11/TGBT21A_501_11_gettp_scan_152_ifnum_0_plnum_0.fits"
        hdu = fits.open(gbtidl_file)
        gbtidl_gettp = hdu[1].data["DATA"][0]
        gbtidl_exp = hdu[1].data["EXPOSURE"][0]
        gbtidl_tsys = hdu[1].data["TSYS"][0]
        tp0 = tps[0].total_power(0)
        diff = tp0.flux.value - gbtidl_gettp
        hdu.close()
        assert np.nanmean(diff) == 0.0
        assert tp0.meta["TSYS"] == pytest.approx(gbtidl_tsys)
        assert tp0.meta["EXPOSURE"] == pytest.approx(gbtidl_exp)

        # Now do some sig=F data.
        sdf_file = f"{self.data_dir}/TGBT21A_504_01/TGBT21A_504_01.raw.vegas/TGBT21A_504_01.raw.vegas.A.fits"
        sdf = gbtfitsload.GBTFITSLoad(sdf_file)
        tps = sdf.gettp(scan=20, ifnum=0, plnum=1, sig=False, cal=True)
        gbtidl_file = (
            f"{self.data_dir}/TGBT21A_504_01/TGBT21A_504_01_gettp_scan_20_ifnum_0_plnum_1_sig_state_0_cal_state_1.fits"
        )
        hdu = fits.open(gbtidl_file)
        gbtidl_gettp = hdu[1].data["DATA"][0]
        gbtidl_exp = hdu[1].data["EXPOSURE"][0]
        gbtidl_tsys = hdu[1].data["TSYS"][0]
        tp0 = tps[0].timeaverage()
        diff = (tp0.flux.value - gbtidl_gettp) / gbtidl_gettp
        hdu.close()
        assert np.nanmean(diff) < 3e-8
        assert tp0.meta["TSYS"] == pytest.approx(gbtidl_tsys)
        assert tp0.meta["EXPOSURE"] == pytest.approx(gbtidl_exp)

        tps = sdf.gettp(scan=20, ifnum=0, plnum=1, sig=False, cal=False)
        gbtidl_file = (
            f"{self.data_dir}/TGBT21A_504_01/TGBT21A_504_01_gettp_scan_20_ifnum_0_plnum_1_sig_state_0_cal_state_0.fits"
        )
        hdu = fits.open(gbtidl_file)
        gbtidl_gettp = hdu[1].data["DATA"][0]
        gbtidl_exp = hdu[1].data["EXPOSURE"][0]
        gbtidl_tsys = hdu[1].data["TSYS"][0]
        tp0 = tps[0].timeaverage()
        diff = (tp0.flux.value - gbtidl_gettp) / gbtidl_gettp
        hdu.close()
        assert np.nanmean(diff) < 1e-7
        assert tp0.meta["TSYS"] == pytest.approx(gbtidl_tsys)
        assert tp0.meta["EXPOSURE"] == pytest.approx(gbtidl_exp)

        tps = sdf.gettp(scan=20, ifnum=0, plnum=1, sig=False, cal=None)
        gbtidl_file = (
            f"{self.data_dir}/TGBT21A_504_01/TGBT21A_504_01_gettp_scan_20_ifnum_0_plnum_1_sig_state_0_cal_all.fits"
        )
        hdu = fits.open(gbtidl_file)
        gbtidl_gettp = hdu[1].data["DATA"][0]
        gbtidl_exp = hdu[1].data["EXPOSURE"][0]
        gbtidl_tsys = hdu[1].data["TSYS"][0]
        tp0 = tps[0].timeaverage()
        diff = (tp0.flux.value - gbtidl_gettp) / gbtidl_gettp
        hdu.close()
        assert np.nanmean(diff) < 1e-7
        assert tp0.meta["TSYS"] == pytest.approx(gbtidl_tsys)
        assert tp0.meta["EXPOSURE"] == pytest.approx(gbtidl_exp)

        # Self consistency check.
        # This only makes sure that the output matches what is expected given the data selection.
        data_path = f"{self.data_dir}/AGBT18B_354_03/AGBT18B_354_03.raw.vegas/AGBT18B_354_03.raw.vegas.A.fits"
        sdf = gbtfitsload.GBTFITSLoad(data_path, verbose=False)
        tests = {
            0: {"SCAN": 6, "IFNUM": 2, "PLNUM": 0, "CAL": None, "SIG": None},
            1: {"SCAN": 6, "IFNUM": 2, "PLNUM": 0, "CAL": True, "SIG": None},
            2: {"SCAN": 6, "IFNUM": 2, "PLNUM": 0, "CAL": False, "SIG": None},
            3: {"SCAN": 6, "IFNUM": 2, "PLNUM": 0, "CAL": None, "SIG": True},
            4: {"SCAN": 6, "IFNUM": 2, "PLNUM": 0, "CAL": None, "SIG": False},
            5: {"SCAN": 6, "IFNUM": 2, "PLNUM": 0, "CAL": True, "SIG": True},
            6: {"SCAN": 6, "IFNUM": 2, "PLNUM": 0, "CAL": True, "SIG": False},
            7: {"SCAN": 6, "IFNUM": 2, "PLNUM": 0, "CAL": False, "SIG": False},
            8: {"SCAN": 6, "IFNUM": 2, "PLNUM": 0, "CAL": False, "SIG": True},
        }
        for k, v in tests.items():
            if v["SIG"] == False:
                with pytest.raises(Exception):
                    tps = sdf.gettp(scan=v["SCAN"], ifnum=v["IFNUM"], plnum=v["PLNUM"], cal=v["CAL"], sig=v["SIG"])
                continue
            tps = sdf.gettp(scan=v["SCAN"], ifnum=v["IFNUM"], plnum=v["PLNUM"], cal=v["CAL"], sig=v["SIG"])
            if v["CAL"]:
                assert np.all(tps[0]._refcalon[0] == tps[0].total_power(0).flux.value)
            tp = tps.timeaverage(weights=None)
            if v["CAL"] is None:
                cal = (0.5 * (tps[0]._refcalon + tps[0]._refcaloff)).astype(np.float64)
            elif not v["CAL"]:
                # CAL=False
                cal = tps[0]._refcaloff.astype(np.float64)
            else:
                # CAL=True
                cal = tps[0]._refcalon.astype(np.float64)
            # diff = tp.flux.value - np.nanmean(cal, axis=0)
            assert np.all(tp.flux.value - np.nanmean(cal, axis=0) == 0)
        # Check that selection is being applied properly.
        tp_scans = sdf.gettp(scan=[6, 7], plnum=0)
        # Weird that the results are different for a bunch of channels.
        # This has to do with slight differences in Tsys weighting in ScanBlock.timeaverage() vs. Scan.timeaverage()
        assert np.all((sdf.gettp(scan=6, plnum=0).timeaverage().flux - tp_scans[0].timeaverage().flux).value < 2e-6)
        assert np.all((sdf.gettp(scan=7, plnum=0).timeaverage().flux - tp_scans[1].timeaverage().flux).value < 2e-6)
        assert np.all(
            (
                sdf.gettp(scan=6, plnum=0).timeaverage(weights=None).flux - tp_scans[0].timeaverage(weights=None).flux
            ).value
            == 0
        )
        assert np.all(
            (
                sdf.gettp(scan=7, plnum=0).timeaverage(weights=None).flux - tp_scans[1].timeaverage(weights=None).flux
            ).value
            == 0
        )

    def test_load_multifits(self):
        """
        Loading multiple SDFITS files under a directory.
        It checks that

        * the GBTFITSLoad object has the correct number of IFs
        * the GBTFITSLoad object has the correct number of polarizations
        * the GBTFITSLoad object has the correct number of beams
        """

        fits_path = f"{self.data_dir}/AGBT18B_354_03/AGBT18B_354_03.raw.vegas"

        sdf = gbtfitsload.GBTFITSLoad(fits_path)

        # Check IFNUMs.
        ifnums = np.sort(sdf.udata("IFNUM"))
        assert np.sum(np.subtract(ifnums, [0, 1, 2, 3])) == 0

        # Check PLNUMs.
        plnums = np.sort(sdf.udata("PLNUM"))
        assert np.sum(np.subtract(plnums, [0, 1])) == 0

        # Check FDNUMs.
        fdnums = np.sort(sdf.udata("FDNUM"))
        assert len(fdnums) == 1
        assert np.sum(np.subtract(fdnums, [0])) == 0

    # @pytest.mark.skip(reason="We need to update this to work with multifits and ScanBlocks")
    def test_multifits_getps_offon(self):
        """
        Loading multiple SDFITS files under a directory.
        It checks that

        * getps works on all IFs
        *
        """
        proj_dir = f"{self.data_dir}/AGBT18B_354_03"
        # Load the data.
        fits_dir = f"{proj_dir}/AGBT18B_354_03.raw.vegas"
        sdf = gbtfitsload.GBTFITSLoad(fits_dir)

        # Check one IF.
        ps_scans = sdf.getps(6, ifnum=2)
        ps_spec = ps_scans[0].calibrated(0).flux.to("K").value

        gbtidl_dir = f"{proj_dir}/gbtidl_outputs"
        hdu = fits.open(f"{gbtidl_dir}/getps_scan_6_ifnum_2_plnum_0_intnum_0.fits")
        table = hdu[1].data
        gbtidl_spec = table["DATA"]

        diff = ps_spec.astype(np.float32) - gbtidl_spec[0]
        hdu.close()
        # assert np.all((ps_spec.astype(np.float32) - gbtidl_spec) == 0)
        assert np.all(abs(diff[~np.isnan(diff)]) < 7e-5)
        assert table["EXPOSURE"] == ps_scans[0].calibrated(0).meta["EXPOSURE"]
        assert np.all(abs(diff[~np.isnan(diff)]) < 7e-5)

    def test_summary(self):
        """Test that some of the columns in the summary
        match the ones produced by `GBTIDL`."""

        def read_gbtidl_summary(filename, idx=1):
            """ """
            with open(filename, "r") as log:
                lines = log.readlines()
            lines.pop(idx)
            tmp = Path(f"{filename}.tmp")
            with open(tmp, "w") as f:
                for line in lines:
                    f.write(line)

            df = pd.read_fwf(f"{filename}.tmp")
            # Clean up.
            tmp.unlink()

            return df

        cols = {
            "SCAN": "Scan",
            "OBJECT": "Source",
            "VELOCITY": "Vel",
            # "PROC": "Proc", # GBTIDL trims the names.
            "PROCSEQN": "Seq",
            "# IF": "nIF",
            "# INT": "nInt",
            "# FEED": "nFd",
        }

        path = util.get_project_testdata() / "AGBT05B_047_01"
        sdf_file = path / "AGBT05B_047_01.raw.acs"
        sdf = gbtfitsload.GBTFITSLoad(sdf_file)
        dysh_df = sdf.summary()

        gbtidl_summary = read_gbtidl_summary(path / "gbtidl" / "AGBT05B_047_01.summary")

        for col in cols.items():
            assert_series_equal(
                dysh_df[col[0]],  # .sort_values(),
                gbtidl_summary[col[1]],  # .sort_values(),
                check_dtype=False,
                check_names=False,
                check_index=False,
            )

    def test_contruct_integration_number(self):
        """Test that construction of integration number (intnum) during FITS load matches
        that in the GBTIDL index file
        """
        p = util.get_project_testdata() / "AGBT20B_014_03.raw.vegas"
        index_file = p / "AGBT20B_014_03.raw.vegas.A6.index"
        data_file = p / "AGBT20B_014_03.raw.vegas.A6.fits"
        g = gbtfitsload.GBTFITSLoad(data_file)
        gbtidl_index = pd.read_csv(index_file, skiprows=10, sep="\s+")
        assert np.all(g._index["INTNUM"] == gbtidl_index["INT"])

    def test_getps_smoothref(self):
        """ """

        path = util.get_project_testdata() / "TGBT21A_501_11"
        data_file = path / "TGBT21A_501_11.raw.vegas.fits"
        gbtidl_file = path / "TGBT21A_501_11_getps_scan_152_ifnum_0_plnum_0_smthoff_15.fits"

        sdf = gbtfitsload.GBTFITSLoad(data_file)
        sb = sdf.getps(scan=152, ifnum=0, plnum=0, smoothref=15)
        ta = sb.timeaverage()

        hdu = fits.open(gbtidl_file)
        table = hdu[1].data
        hdu.close()
        gbtidl_spec = table["DATA"][0]

        diff = ta.flux.to("K").value.astype(np.float32) - gbtidl_spec

        assert np.all(abs(diff[~np.isnan(diff)]) < 7e-5)
        assert ta.meta["EXPOSURE"] == table["EXPOSURE"][0]
        for k, v in ta.meta.items():
            if k in ["DURATION", "TUNIT7", "VSPRPIX", "CAL"]:
                continue
            try:
                assert v == table[k][0]
            except KeyError:
                continue

    def test_getps_flagging(self):
        path = util.get_project_testdata() / "TGBT21A_501_11"
        data_file = path / "TGBT21A_501_11.raw.vegas.fits"
        sdf = gbtfitsload.GBTFITSLoad(data_file)
        sdf.flag_channel([[10, 20], [30, 41]])
        sb = sdf.getps(scan=152, ifnum=0, plnum=0, apply_flags=True)
        ta = sb.timeaverage()
        # average_spectra masks out the NaN in channel 3072
        expected_mask = np.hstack([np.arange(10, 21), np.arange(30, 42), np.array([3072])])
        assert np.all(np.where(ta.mask) == expected_mask)

    def test_write_single_file(self, tmp_path):
        "Test that writing an SDFITS file works when subselecting data"
        p = util.get_project_testdata() / "AGBT20B_014_03.raw.vegas"
        data_file = p / "AGBT20B_014_03.raw.vegas.A6.fits"
        g = gbtfitsload.GBTFITSLoad(data_file)
        o = tmp_path / "sub"
        o.mkdir()
        out = o / "test_write_single.fits"
        g.write(out, plnum=1, intnum=2, overwrite=True)
        t = gbtfitsload.GBTFITSLoad(out)
        assert set(t._index["PLNUM"]) == set([1])
        # assert set(t._index["INT"]) == set([2])  # this exists because GBTIDL wrote it
        assert set(t._index["INTNUM"]) == set([2])

    def test_write_multi_file(self, tmp_path):
        "Test that writing multiple SDFITS files works, including subselection of data"
        f = util.get_project_testdata() / "AGBT18B_354_03/AGBT18B_354_03.raw.vegas/"
        g = gbtfitsload.GBTFITSLoad(f)
        # writes testmultia0,1,2,3.fits
        o = tmp_path / "sub"
        o.mkdir()
        output = o / "testmulti.fits"
        g.write(output, multifile=True, scan=6, overwrite=True)
        # @todo remove test output files in a teardown method
        sdf = gbtfitsload.GBTFITSLoad(o)
        assert set(sdf["SCAN"]) == set([6])

    def test_write_all(self, tmp_path):
        """Test that we can write a loaded SDFITS file without any changes"""
        p = util.get_project_testdata() / "AGBT20B_014_03.raw.vegas"
        data_file = p / "AGBT20B_014_03.raw.vegas.A6.fits"
        org_sdf = gbtfitsload.GBTFITSLoad(data_file)
        d = tmp_path / "sub"
        d.mkdir()
        output = d / "test_write_all.fits"
        # don't write flags to avoid TDIM84 new column
        org_sdf.write(output, overwrite=True, flags=False)
        new_sdf = gbtfitsload.GBTFITSLoad(output)
        # Compare the index for both SDFITS.
        # Note we now auto-add a HISTORY card at instantiation, so drop that
        # from the comparison
        assert_frame_equal(org_sdf._index, new_sdf._index.drop(columns="HISTORY"))

    def test_get_item(self):
        f = util.get_project_testdata() / "AGBT18B_354_03/AGBT18B_354_03.raw.vegas/"
        g = gbtfitsload.GBTFITSLoad(f)
        assert list(set(g["PLNUM"])) == [0, 1]
        assert list(set(g[["SCAN", "IFNUM"]].loc[0])) == [2, 6]
        # test case insensitivity
        assert list(set(g["plnum"])) == [0, 1]
        with pytest.raises(KeyError):
            g["FOOBAR"]

    def test_set_item(self, tmp_path):
        # First test with a single number or string
        keyval = {
            "IFNUM": 12,
            "FRONTEND": "Rx999",
            "sitelong": 42.21,
        }
        d = util.get_project_testdata()
        files = [
            d / "AGBT20B_014_03.raw.vegas/AGBT20B_014_03.raw.vegas.A6.fits",  # single file
            d / "AGBT18B_354_03/AGBT18B_354_03.raw.vegas/",  # multiple files
        ]
        o = tmp_path / "gsetitem"
        o.mkdir()

        i = 0
        for f in files:
            g = gbtfitsload.GBTFITSLoad(f)
            for key, val in keyval.items():
                _set = set([val])
                g[key] = val
                assert set(g[key]) == _set
                for sdf in g._sdf:
                    assert set(sdf[key]) == _set
                    for b in sdf._bintable:
                        assert set(b.data[key]) == _set

            # check that thing were written correctly.
            out = o / f"test_write_gsetitem{i}.fits"
            print(f"trying to writing file #{i} {out}")
            g.write(out, overwrite=True, flags=False)
            i += 1
            if "A6" in f.name:
                g = gbtfitsload.GBTFITSLoad(out)
            else:
                g = gbtfitsload.GBTFITSLoad(o)
            for key, val in keyval.items():
                _set = set([val])
                g[key] = val
                assert set(g[key]) == _set
                for sdf in g._sdf:
                    assert set(sdf[key]) == _set
                    for b in sdf._bintable:
                        assert set(b.data[key]) == _set

        # now test array of numbers or strings
        for f in files:
            g = gbtfitsload.GBTFITSLoad(f)
            for key, val in keyval.items():
                array = [val] * g.total_rows
                _set = set([val])
                g[key] = array
                assert set(g[key]) == _set
                for sdf in g._sdf:
                    for b in sdf._bintable:
                        assert set(b.data[key]) == _set

            # check that thing were written correctly.
            out = o / f"test_write_gsetitem{i}.fits"
            print(f"trying to writing file #{i} {out}")
            g.write(out, overwrite=True, flags=False)
            i += 1
            if "A6" in f.name:
                g = gbtfitsload.GBTFITSLoad(out)
            else:
                g = gbtfitsload.GBTFITSLoad(o)
            for key, val in keyval.items():
                _set = set([val])
                g[key] = val
                assert set(g[key]) == _set
                for sdf in g._sdf:
                    assert set(sdf[key]) == _set
                    for b in sdf._bintable:
                        assert set(b.data[key]) == _set

        # check that exception is handled for incorrect length
        for f in files:
            g = gbtfitsload.GBTFITSLoad(f)
            for key, val in keyval.items():
                array = [val] * 2 * g.total_rows
                with pytest.raises(ValueError):
                    g[key] = array

        # test that changed a previously selection column results in a warning
        g = gbtfitsload.GBTFITSLoad(files[0])
        g.select(ifnum=2)
        with pytest.warns(UserWarning):
            g["ifnum"] = 3

        def test_data_access(self):
            """test getting and setting the DATA column of SDFITS"""
            # File with a single BinTableHDU
            d = util.get_project_testdata()
            f = d / "AGBT18B_354_03/AGBT18B_354_03.raw.vegas/AGBT18B_354_03.raw.vegas.A.fits"
            g = gbtfitsload.GBTFITSLoad(f)
            data = g["DATA"]
            assert data.shape == (32, 131072)
            g["DATA"] = np.zeros([32, 131072])
            assert np.all(g["DATA"] == 0)

            # File with multiple BinTableHDUs
            f = d / "TGBT17A_506_11/TGBT17A_506_11.raw.vegas.A_truncated_rows.fits"
            g = gbtfitsload.GBTFITSLoad(f)
            # The binary tables have different shapes, so setting and getting is not allowed.
            with pytest.raises(Exception):
                g["DATA"]
            with pytest.raises(Exception):
                g["DATA"] = np.random.rand(1024)
            # Multiple files
            f = util.get_project_testdata() / "AGBT18B_354_03/AGBT18B_354_03.raw.vegas/"
            g = gbtfitsload.GBTFITSLoad(f)

    def test_azel_coords(self, tmp_path):
        """
        Test that observations using AzEl coordinates can produce a valid `Spectrum`.
        """

        # Reuse an existing file in testdata.
        fits_path = util.get_project_testdata() / "AGBT18B_354_03/AGBT18B_354_03.raw.vegas"
        new_path = tmp_path / "o"
        new_path.mkdir(parents=True)
        sdf_org = gbtfitsload.GBTFITSLoad(fits_path)
        sdf_org["RADESYS"] = ""
        sdf_org["CTYPE2"] = "AZ"
        sdf_org["CTYPE3"] = "EL"

        # Create a temporary directory and write the modified SDFITS.
        new_path = tmp_path / "o"
        new_path.mkdir(parents=True, exist_ok=True)
        sdf_org.write(new_path / "test_azel.fits", overwrite=True, flags=True)

        # Now the actual test.
        sdf = gbtfitsload.GBTFITSLoad(new_path)
        # Not this part of the test, but just to make sure.
        assert np.all(sdf["RADESYS"] == "AltAz")
        # Test that we can create a `Spectrum` object.
        tp = sdf.gettp(scan=6, plnum=0)[0].total_power(0)

    def test_hadec_coords(self, tmp_path):
        """
        Test that observations using HADec coordinates can produce a valid `Spectrum`.
        """

        # Reuse an existing file in testdata.
        fits_path = util.get_project_testdata() / "AGBT18B_354_03/AGBT18B_354_03.raw.vegas"
        new_path = tmp_path / "o"
        new_path.mkdir(parents=True)
        sdf_org = gbtfitsload.GBTFITSLoad(fits_path)
        sdf_org["RADESYS"] = ""
        sdf_org["CTYPE2"] = "HA"

        # Create a temporary directory and write the modified SDFITS.
        new_path = tmp_path / "o"
        new_path.mkdir(parents=True, exist_ok=True)
        sdf_org.write(new_path / "test_hadec.fits", overwrite=True)

        # Now the actual test.
        sdf = gbtfitsload.GBTFITSLoad(new_path)
        # Not this part of the test, but just to make sure.
        assert np.all(sdf["RADESYS"] == "hadec")
        # Test that we can create a `Spectrum` object.
        tp = sdf.gettp(scan=6, plnum=0)[0].total_power(0)

    def test_galactic_coords(self, tmp_path):
        """
        Test that observations using Galactic coordinates can produce a valid `Spectrum`.
        """

        # Reuse an existing file in testdata.
        fits_path = util.get_project_testdata() / "AGBT18B_354_03/AGBT18B_354_03.raw.vegas"
        new_path = tmp_path / "o"
        new_path.mkdir(parents=True)
        sdf_org = gbtfitsload.GBTFITSLoad(fits_path)
        sdf_org["RADESYS"] = ""
        sdf_org["CTYPE2"] = "GLON"
        sdf_org["CTYPE3"] = "GLAT"

        # Create a temporary directory and write the modified SDFITS.
        new_path = tmp_path / "o"
        new_path.mkdir(parents=True, exist_ok=True)
        sdf_org.write(new_path / "test_galactic.fits", overwrite=True)

        # Now the actual test.
        sdf = gbtfitsload.GBTFITSLoad(new_path)
        # Not this part of the test, but just to make sure.
        assert np.all(sdf["RADESYS"] == "galactic")
        # Test that we can create a `Spectrum` object.
        tp = sdf.gettp(scan=6, plnum=0)[0].total_power(0)

    def test_add_history_comments(self):
        fits_path = util.get_project_testdata() / "AGBT18B_354_03/AGBT18B_354_03.raw.vegas"
        sdf = gbtfitsload.GBTFITSLoad(fits_path)
        sb = sdf.gettp(scan=6, plnum=0)
        sdf.add_comment("My dear Aunt Sally")
        sdf.add_history("ran the test for history and comments")
        assert "My dear Aunt Sally" in sdf.comments
        assert "ran the test for history and comments" in sdf.history
        assert any("Project ID: AGBT18B_354_03" in substr for substr in sb.history)

    def test_online(self, tmp_path):
        f1 = util.get_project_testdata() / "TGBT21A_501_11/TGBT21A_501_11.raw.vegas.fits"
        f2 = util.get_project_testdata() / "TGBT21A_501_11/TGBT21A_501_11_2.raw.vegas.fits"
        #
        sdfits = tmp_path / "sdfits"
        sdfits.mkdir()
        os.environ["SDFITS_DATA"] = str(sdfits)
        print("PJT1", sdfits)
        o1 = sdfits / "online.fits"
        print("PJT2", o1)
        #
        shutil.copyfile(f1, o1)
        sdf = gbtfitsload.GBTOnline()
        s = sdf.summary()
        n = len(sdf._index)
        assert n == 4
        if sdf._platform == "Windows":
            # os.remove(o1)
            # pathlib.Path.unlink(o1)
            print("Windows seems to lock the file, can't remover or overwite")
        else:
            shutil.copyfile(f2, o1)
            s = sdf.summary()
            n = len(sdf._index)
            assert n == 8

    def test_write_read_flags(self, tmp_path):
        fits_path = util.get_project_testdata() / "AGBT18B_354_03/AGBT18B_354_03.raw.vegas"
        sdf = gbtfitsload.GBTFITSLoad(fits_path)
        sdf.flag(ifnum=3)
        sdf.apply_flags()

        new_path = tmp_path / "flagtest_multi"
        new_path.mkdir(parents=True, exist_ok=True)
        # first test with multifile=True. Default is flags=True
        sdf.write(new_path / "test.fits", multifile=True, overwrite=True)
        mdf = gbtfitsload.GBTFITSLoad(new_path)
        for i in range(len(sdf._sdf)):
            assert np.all(sdf._sdf[i]._flagmask[0] == mdf._sdf[i]._flagmask[0])
        # now with multifile = False
        new_path = tmp_path / "flagtest_single"
        new_path.mkdir(parents=True, exist_ok=True)
        sdf.write(new_path / "foobar.fits", multifile=False, overwrite=True)
        g2 = gbtfitsload.GBTFITSLoad(new_path / "foobar.fits", verbose=True)
        flags2 = g2._sdf[0]._flagmask
        for i in range(len(sdf._sdf)):
            assert np.all(sdf._sdf[i]._flagmask[0] == flags2[i])

    def test_additive_flags(self, tmp_path):
        # This test is a regression for issue #429
        # https://github.com/GreenBankObservatory/dysh/issues/429
        # We copy the flag mask for each rule
        # then check that the final flag mask is the logical OR of them
        fits_path = (
            util.get_project_testdata() / "AGBT18B_354_03/AGBT18B_354_03.raw.vegas/AGBT18B_354_03.raw.vegas.A.fits"
        )
        sdf = gbtfitsload.GBTFITSLoad(fits_path, skipflags=True)
        sdf.flag(scan=6, channel=[[2500, 3000]], intnum=range(0, 3))
        sdf.apply_flags()
        flag1 = sdf._sdf[0]._flagmask[0].copy()
        sdf.clear_flags()
        sdf.flag_channel(channel=[[80000, 100000]])
        sdf.apply_flags()
        flag2 = sdf._sdf[0]._flagmask[0].copy()
        sdf.clear_flags()
        sdf.flag(scan=6, channel=[[2500, 3000]], intnum=range(0, 3))
        sdf.flag_channel(channel=[[80000, 100000]])
        sdf.apply_flags()
        assert np.all(sdf._sdf[0]._flagmask[0] == flag1 | flag2)

    def test_read_gbtidl_flags(self):
        """
        Test reading a flag file generated by GBTIDL.
        """
        fits_path = util.get_project_testdata() / "AGBT17A_404_01/AGBT17A_404_01.raw.vegas"
        sdf = gbtfitsload.GBTFITSLoad(fits_path)
        # The data should not be flagged, as the flags are only for intnum=arange(43,52)
        ps = sdf.getps(scan=19, plnum=0, apply_flags=True).timeaverage()
        assert np.all(ps.mask == False)
        # The data should be flagged for these integrations.
        ps = sdf.getps(scan=19, plnum=0, apply_flags=True, intnum=[i for i in range(43, 52)]).timeaverage()
        assert np.all(ps.mask[2299:] == True)
        assert np.all(ps.mask[:2299] == False)

    def test_rawspectrum(self):
        """regression test for issue 442"""
        fits_path = util.get_project_testdata() / "AGBT05B_047_01/AGBT05B_047_01.raw.acs/AGBT05B_047_01.raw.acs.fits"
        sdf = gbtfitsload.GBTFITSLoad(fits_path)
        psscan = sdf.getps(plnum=0)
        spec = psscan.timeaverage()
        exp1 = spec.meta["EXPOSURE"]  # 214.86978780485427
        sdf.flag_range(elevation=((None, 18.4)))
        psscan2 = sdf.getps(plnum=0)
        spec2 = psscan2.timeaverage()
        exp2 = spec2.meta["EXPOSURE"]  # 58.59014643665782
        assert exp2 < exp1

    def test_getnod_wcal(self):
        """
        Test for getnod using data with noise diode.
        """

        # Reduce with dysh.
        fits_path = util.get_project_testdata() / "TGBT22A_503_02/TGBT22A_503_02.raw.vegas"
        sdf = gbtfitsload.GBTFITSLoad(fits_path)
        nodsb = sdf.getnod(scan=62, ifnum=0, plnum=0)
        nodsp0 = nodsb[0].timeaverage()
        nodsp1 = nodsb[1].timeaverage()

        # Load GBTIDL reduction.
        # row 0 is `fdnum=2`.
        # row 1 is `fdnum=6`.
        hdu = fits.open(util.get_project_testdata() / "TGBT22A_503_02/TGBT22A_503_02.cal.vegas.fits")
        table = hdu[1].data

        # Compare.
        assert nodsp0.meta["EXPOSURE"] == pytest.approx(table["EXPOSURE"][0])
        assert nodsp1.meta["EXPOSURE"] == pytest.approx(table["EXPOSURE"][1])
        # These assert internally.
        np.testing.assert_allclose(nodsp0.data, table["DATA"][0], rtol=2e-7, equal_nan=False)
        np.testing.assert_allclose(nodsp1.data, table["DATA"][1], rtol=2e-7, equal_nan=False)
        assert table["TSYS"][0] == pytest.approx(nodsp0.meta["TSYS"])
        assert table["TSYS"][1] == pytest.approx(nodsp1.meta["TSYS"])

    def test_getnod_nocal(self):
        """
        Test for getnod using data without noise diode.
        """

        # Reduce with dysh.
        fits_path = util.get_project_testdata() / "TSCAL_220105_W/TSCAL_220105_W.raw.vegas"
        sdf = gbtfitsload.GBTFITSLoad(fits_path)
        nodsb = sdf.getnod(scan=24, ifnum=0, plnum=0)
        nodsp0 = nodsb[0].timeaverage()
        nodsp1 = nodsb[1].timeaverage()

        # Load GBTIDL reduction.
        # row 0 is `fdnum=0`.
        # row 1 is `fdnum=1`.
        hdu = fits.open(util.get_project_testdata() / "TSCAL_220105_W/TSCAL_220105_W.cal.vegas.fits")
        table = hdu[1].data

        # Compare.
        assert nodsp0.meta["EXPOSURE"] == pytest.approx(table["EXPOSURE"][0])
        assert nodsp1.meta["EXPOSURE"] == pytest.approx(table["EXPOSURE"][1])
        # These assert internally.
        np.testing.assert_allclose(nodsp0.data, table["DATA"][0], rtol=2e-7, equal_nan=False)
        np.testing.assert_allclose(nodsp1.data, table["DATA"][1], rtol=2e-7, equal_nan=False)
        assert table["TSYS"][0] == pytest.approx(nodsp0.meta["TSYS"])
        assert table["TSYS"][1] == pytest.approx(nodsp1.meta["TSYS"])

    def test_subbeamnod(self):
        """simple check of subbeamnod for two different cases.  this mimics the notebook example"""
        sdf_file = f"{self.data_dir}/AGBT13A_124_06/AGBT13A_124_06.raw.acs/AGBT13A_124_06.raw.acs.fits"
        sdf = gbtfitsload.GBTFITSLoad(sdf_file)
        # don't scale
        sb = sdf.subbeamnod(scan=44, fdnum=1, ifnum=0, plnum=0, method="cycle")
        sb2 = sdf.subbeamnod(scan=44, fdnum=1, ifnum=0, plnum=0, method="scan")
        s = sb.timeaverage() - sb2.timeaverage()
        assert np.nanmean(s.data) == pytest.approx(0.0022912487, abs=1e-8)

    def test_scale(self):
        # Check that scaling to Ta* or Jy works.
        # The code that computes the scale factors is tested in test_gaincorrection.py, so
        # we don't need to recheck that here.
        # PSScan
        sdf_file = f"{self.data_dir}/TGBT21A_501_11/TGBT21A_501_11.raw.vegas.fits"
        sdf = gbtfitsload.GBTFITSLoad(sdf_file)
        sba = sdf.getps(scan=152, bunit="ta*", zenith_opacity=0.05)
        sbb = sdf.getps(scan=152, bunit="jy", zenith_opacity=0.05)
        # The ratio of these scale factors should be the Jy/K of the telescope
        jyperk = sbb[0].bscale / sba[0].bscale
        gc = util.gaincorrection.GBTGainCorrection()
        assert jyperk == pytest.approx(gc.jyperk.value, 1e-6)
        assert sba[0].bunit == "ta*"
        assert sbb[0].bunit == "jy"
        assert sba[0].is_scaled
        assert sbb[0].is_scaled
        # Now test scaling after the fact
        sbd = sdf.getps(scan=152)
        sbd[0].scale("jy", zenith_opacity=0.1)
        assert sbd[0].bunit == "jy"
        assert sbd[0].is_scaled

        # try a bad scale type
        with pytest.raises(ValueError):
            sba = sdf.getps(scan=152, bunit="foobar", zenith_opacity=0.05)
        # try a bad tau
        with pytest.raises(ValueError):
            sba = sdf.getps(scan=152, bunit="jy", zenith_opacity=-1)

        # test that scaling a ScanBlock works, also case insensitivity
        sb = sdf.getps(scan=152, bunit="Ta*", zenith_opacity=0.05)
        assert sb.bunit == "ta*"
        with pytest.raises(ValueError):
            sb.scale("not a valid bunit", zenith_opacity=0.2)

    def test_getps_nocal(self):
        """
        Test for `getps` without noise diodes.
        """

        sdf_file = f"{self.data_dir}/TGBT24B_613_15/TGBT24B_613_15.raw.vegas"
        sdf = gbtfitsload.GBTFITSLoad(sdf_file)
        ps = sdf.getps(scan=23, fdnum=10, plnum=0).timeaverage()

        # Load GBTIDL reduction.
        hdu = fits.open(util.get_project_testdata() / "TGBT24B_613_15/TGBT24B_613_15.cal.vegas.fits")
        table = hdu[1].data

        # Compare.
        assert ps.meta["EXPOSURE"] == pytest.approx(table["EXPOSURE"][0])
        assert table["TSYS"][0] == pytest.approx(ps.meta["TSYS"])
        np.testing.assert_allclose(
            ps.data[~np.isnan(table["DATA"][0])],
            table["DATA"][0][~np.isnan(table["DATA"][0])],
            rtol=2e-7,
            equal_nan=False,
        )
