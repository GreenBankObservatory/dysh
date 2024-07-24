import glob
import os

import numpy as np
import pytest

from dysh import util
from dysh.fits.sdfitsload import SDFITSLoad


class TestSDFITSLoad:
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
            sdf = SDFITSLoad(fnm)
            assert len(sdf._index) == expected[filename]

    def test_getspec(self):
        """
        Test that a SDFITSLoad object can use the `getspec` function.
        """

        sdf_file = f"{self.data_dir}/TGBT21A_501_11/TGBT21A_501_11.raw.vegas.fits"
        sdf = SDFITSLoad(sdf_file)
        spec = sdf.getspec(0)
        assert np.nanmean(spec.data) == pytest.approx(510458900.0)
        assert spec.meta["BACKEND"] == "VEGAS"

    def test_write_single_file(self, tmp_path):
        "Test that writing an SDFITS file works when subselecting data"
        p = util.get_project_testdata() / "AGBT20B_014_03.raw.vegas"
        data_file = p / "AGBT20B_014_03.raw.vegas.A6.fits"
        g = SDFITSLoad(data_file)
        o = tmp_path / "sub"
        o.mkdir()
        out = o / "test_write_single.fits"
        g.write(out, overwrite=True)
        t = SDFITSLoad(out)
        assert set(t["PLNUM"]) == set([0, 1])
        assert set(t["INT"]) == set([0, 1, 2])  # this exists because GBTIDL wrote it

    def test_get_item(self):
        f = util.get_project_testdata() / "AGBT18B_354_03/AGBT18B_354_03.raw.vegas/AGBT18B_354_03.raw.vegas.A.fits"
        g = SDFITSLoad(f)
        assert list(set(g["PLNUM"])) == [0, 1]
        assert list(set(g[["SCAN", "IFNUM"]].loc[0])) == [2, 6]
        with pytest.raises(KeyError):
            g["FOOBAR"]

    def test_set_item(self):
        # File with a single BinTableHDU
        f = util.get_project_testdata() / "AGBT18B_354_03/AGBT18B_354_03.raw.vegas/AGBT18B_354_03.raw.vegas.A.fits"
        g = SDFITSLoad(f)
        # all rows of a column to a single value
        g["FREQRES"] = 1500.0  # Hz
        # test that the selection (index) was set
        assert list(set(g["FREQRES"]))[0] == 1500
        # test that the BinTableHDU data was set
        assert np.all(g._bintable[0].data["FREQRES"] == 1500.0)
        # rows of a column to different values
        x = 3.1415 * np.arange(32, dtype=np.float64)
        # make sure lower/varying case also works
        g["dopfreq"] = x
        assert np.all(g["DoPFreQ"] == x)
        assert np.all(g._bintable[0].data["DOPFREQ"] == x)
        # Wrong length array (except single value which sets all rows) should raise ValueError.
        # We re-raise with additional context as Exception.
        with pytest.raises(Exception):
            g["TWARM"] = np.arange(99)

        # File with multiple BinTableHDUs

    def test_add_bintable_column(self):
        # File with a single BinTableHDU
        assert True
        # File with multiple BinTableHDUs
        assert True
