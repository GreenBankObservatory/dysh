import glob
import os

import numpy as np
import pytest

from dysh import util
from dysh.fits.sdfitsload import SDFITSLoad


class TestSDFITSLoad:
    """ """

    def setup_method(self):
        self.data_dir = util.get_project_testdata()
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
            sdf = SDFITSLoad(fnm)
            assert len(sdf._index) == expected[filename]

    def test_names(self):
        """
        Test basic filename
        """
        filename = "TGBT21A_501_11.raw.vegas.fits"  # noqa: F841
        fnm = self._file_list[0]  # note fnm is a string
        sdf = SDFITSLoad(fnm)
        assert fnm == sdf.filename

    def test_getspec(self):
        """
        Test that a SDFITSLoad object can use the `getspec` function.
        """
        index = 2
        sdf_file = f"{self.data_dir}/TGBT21A_501_11/TGBT21A_501_11.raw.vegas.fits"
        sdf = SDFITSLoad(sdf_file)
        spec = sdf.getspec(index)
        assert spec.flux.unit == "ct"
        assert np.nanmean(spec.data) == pytest.approx(504480960.0)
        assert spec.meta["BACKEND"] == "VEGAS"
        assert spec.meta["ROW"] == index
        assert spec.meta["CRVAL4"] == -6

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

    def test_set_item(self, tmp_path):
        # File with a single BinTableHDU
        d = util.get_project_testdata()
        f = d / "AGBT18B_354_03/AGBT18B_354_03.raw.vegas/AGBT18B_354_03.raw.vegas.A.fits"
        g = SDFITSLoad(f)
        # all rows of a column to a single value
        with pytest.warns(UserWarning):
            g["FREQRES"] = 1500.0  # Hz
        # test that the selection (index) was set
        assert list(set(g["FREQRES"]))[0] == 1500  # noqa: RUF015
        # test that the BinTableHDU data was set
        assert np.all(g._bintable[0].data["FREQRES"] == 1500.0)
        # rows of a column to different values
        x = 3.1415 * np.arange(g.nrows(0), dtype=np.float64)
        # make sure lower/varying case also works
        with pytest.warns(UserWarning):
            g["dopfreq"] = x
        assert np.all(g["DoPFreQ"] == x)
        assert np.all(g._bintable[0].data["DOPFREQ"] == x)
        # Wrong length array (except single value which sets all rows) should raise ValueError.
        # We re-raise with additional context as Exception.
        with pytest.warns(UserWarning):
            with pytest.raises(Exception):
                g["TWARM"] = np.arange(g.nrows(0) + 99)

        # File with multiple BinTableHDUs
        # This is a rare case and in any event the multiple bintables make likely  have
        # different shapes, in which case we can't do setitem operations.
        # This file has 74 columns in each bintable, 3 rows in first table and 5 rows in 2nd
        # (after adding primary HDU, it will have 91 columns)
        f = d / "TGBT17A_506_11/TGBT17A_506_11.raw.vegas.A_truncated_rows.fits"
        g = SDFITSLoad(f)
        # first test setting all rows to same value
        c = "MBM12"
        with pytest.warns(UserWarning):
            g["OBJECT"] = c
        assert np.all(g["object"] == c)
        assert np.all(g.bintable[0].data["OBJECT"] == c)
        assert np.all(g.bintable[1].data["OBJECT"] == c)
        # now an array
        c = ["NGC123"] * g.nrows(0) + ["3C111"] * g.nrows(1)
        with pytest.warns(UserWarning):
            g["object"] = c
        assert np.all(g["object"] == c)
        assert np.all(g.bintable[0].data["OBJECT"] == c[0 : g.nrows(0)])
        assert np.all(g.bintable[1].data["OBJECT"] == c[g.nrows(0) :])
        c.append("ONETOOMANY")
        with pytest.warns(UserWarning):
            with pytest.raises(Exception):
                g["object"] = c
        # now a number
        num = 1.23e9
        with pytest.warns(UserWarning):
            g["RESTFREQ"] = num
        assert np.all(g["restfreq"] == num)
        assert np.all(g.bintable[0].data["RESTFREQ"] == num)
        assert np.all(g.bintable[1].data["RESTFREQ"] == num)
        # now a sequence of numbers
        num = np.arange(8) * 9.0e9
        with pytest.warns(UserWarning):
            g["RESTFREQ"] = num
        assert np.all(g["restfreq"] == num)
        assert np.all(g.bintable[0].data["restfreq"] == num[0 : g.nrows(0)])  # wow astropy allows lowercase
        assert np.all(g.bintable[1].data["restfreq"] == num[g.nrows(0) :])
        # check that the change gets written
        o = tmp_path / "setitem"
        o.mkdir()
        out = o / "test_write_setitem.fits"
        g.write(out, overwrite=True)
        g = SDFITSLoad(out)
        assert np.all(g["object"] == c[0 : g.total_rows])
        assert np.all(g.bintable[0].data["OBJECT"] == c[0 : g.nrows(0)])
        assert np.all(g.bintable[1].data["OBJECT"] == c[g.nrows(0) : g.total_rows])
        assert np.all(g["restFREQ"] == num)
        assert np.all(g.bintable[0].data["RESTFREQ"] == num[0 : g.nrows(0)])
        assert np.all(g.bintable[1].data["RESTFREQ"] == num[g.nrows(0) :])

    def test_rename_column(self, tmp_path):
        # File with a single BinTableHDU
        d = util.get_project_testdata()
        f = d / "AGBT18B_354_03/AGBT18B_354_03.raw.vegas/AGBT18B_354_03.raw.vegas.A.fits"
        g = SDFITSLoad(f)
        g.rename_column("object", "target")
        assert "TARGET" in g.columns
        assert "TARGET" in g._bintable[0].data.names
        o = tmp_path / "rename"
        o.mkdir()
        out = o / "test_write_setitem.fits"
        g.write(out, overwrite=True)
        g = SDFITSLoad(out)
        assert "TARGET" in g.columns
        assert "TARGET" in g._bintable[0].data.names

        # File with multiple BinTableHDUs
        f = d / "TGBT17A_506_11/TGBT17A_506_11.raw.vegas.A_truncated_rows.fits"
        g = SDFITSLoad(f)
        g.rename_column(
            "AZIMUTH", "THESPINNYDIRECTION"
        )  # believe it or not, FITS allows TTYPEs to be more than 8 chars
        assert "THESPINNYDIRECTION" in g.columns
        for b in g.bintable:
            assert "THESPINNYDIRECTION" in b.data.names
        # check that the change gets written
        out = o / "test_write_setitem.fits"
        g.write(out, overwrite=True)
        g = SDFITSLoad(out)
        assert "THESPINNYDIRECTION" in g.columns
        for b in g.bintable:
            assert "THESPINNYDIRECTION" in b.data.names

    def test_add_column(self, tmp_path):
        # File with a single BinTableHDU
        d = util.get_project_testdata()
        f = d / "AGBT18B_354_03/AGBT18B_354_03.raw.vegas/AGBT18B_354_03.raw.vegas.A.fits"
        g = SDFITSLoad(f)
        # single value
        num = 7.2e-6
        meow = "Purrrrr"
        # multivalue
        e = 2.718
        val = e * np.arange(g.nrows(0))
        c = ["boy wonder"] * g.nrows(0)
        colval = {"superman": num, "catwoman": meow, "batman": val, "robin": c}
        for k, v in colval.items():
            g[k] = v
            assert np.all(g[k] == v)
            assert np.all(g.bintable[0].data[k] == v)
        with pytest.raises(Exception):
            c.append("nope")
            with pytest.warns(UserWarning):
                g["robin"] = c
        c = ["boy wonder"] * g.nrows(0)
        colval["robin"] = c
        # write and check
        o = tmp_path / "addcol"
        o.mkdir()
        out = o / "test_write_setitem.fits"
        g.write(out, overwrite=True)
        g = SDFITSLoad(out)
        for k, v in colval.items():
            assert np.all(g[k] == v)
            assert np.all(g.bintable[0].data[k] == v)

        # File with multiple BinTableHDUs
        f = d / "TGBT17A_506_11/TGBT17A_506_11.raw.vegas.A_truncated_rows.fits"
        g = SDFITSLoad(f)
        val = e * np.arange(g.total_rows)
        c = ["boy wonder"] * g.total_rows
        colval = {"superman": num, "catwoman": meow, "batman": val, "robin": c}
        for k, v in colval.items():
            g[k] = v
            assert np.all(g[k] == v)
        #    assert np.all(g.bintable[0].data[k] == v)
        #    assert np.all(g.bintable[1].data[k] == v)
        with pytest.raises(Exception):
            c.append("nope")
            with pytest.warns(UserWarning):
                g["robin"] = c
        c = ["boy wonder"] * g.total_rows
        colval["robin"] = c

    def test_delete_column(self):
        # File with a single BinTableHDU
        d = util.get_project_testdata()
        f = d / "AGBT18B_354_03/AGBT18B_354_03.raw.vegas/AGBT18B_354_03.raw.vegas.A.fits"
        g = SDFITSLoad(f)
        with pytest.warns(UserWarning):
            g.delete_column("dopfreq")
        assert "DOPFREQ" not in g.columns
        assert "DOPFREQ" not in g._bintable[0].data.names

        # File with multiple BinTableHDUs
        f = d / "TGBT17A_506_11/TGBT17A_506_11.raw.vegas.A_truncated_rows.fits"
        g = SDFITSLoad(f)
        with pytest.warns(UserWarning):
            g.delete_column("TWARM")
        assert "TWARM" not in g.columns
        for b in g._bintable:
            assert "TWARM" not in b.data.names

        # DATA columns can't be deleted
        with pytest.warns(UserWarning):
            with pytest.raises(Exception):
                g.delete_column("dAtA")

    def test_data_access(self):
        """test getting and setting the DATA column of SDFITS"""
        # File with a single BinTableHDU
        d = util.get_project_testdata()
        f = d / "AGBT18B_354_03/AGBT18B_354_03.raw.vegas/AGBT18B_354_03.raw.vegas.A.fits"
        g = SDFITSLoad(f)
        data = g["DATA"]
        assert data.shape == (32, 131072)
        with pytest.warns(UserWarning):
            g["DATA"] = np.zeros([32, 131072])
        assert np.all(g["DATA"] == 0)
        # test some slicing
        assert np.shape(g["DATA"][:, 0:10]) == (32, 10)
        # assignment via slicing doesn't work.
        # I think because copies of the data are being made
        # g["DATA"][0][10:20] = np.random.rand(10)

        # File with multiple BinTableHDUs
        f = d / "TGBT17A_506_11/TGBT17A_506_11.raw.vegas.A_truncated_rows.fits"
        g = SDFITSLoad(f)
        # The binary tables have different shapes, so setting and getting is not allowed.
        with pytest.raises(Exception):
            g["DATA"]
        with pytest.warns(UserWarning):
            with pytest.raises(Exception):
                g["DATA"] = np.random.rand(1024)

    def test_udata(self):
        """Test that `udata` is working."""

        f = util.get_project_testdata() / "AGBT05B_047_01/AGBT05B_047_01.raw.acs/AGBT05B_047_01.raw.acs.fits"
        g = SDFITSLoad(f)
        assert g.udata("SCAN") == [51, 52, 53, 54, 55, 56, 57, 58]

    def test_ushow(self):
        """Test that `ushow` works."""

        f = util.get_project_testdata() / "AGBT05B_047_01/AGBT05B_047_01.raw.acs/AGBT05B_047_01.raw.acs.fits"
        g = SDFITSLoad(f)
        g.ushow("SCAN")

    def test_column_mask(self):
        """
        Check that _column_mask behaves.
        """

        f = util.get_project_testdata() / "TRCO_230413_Ka/TRCO_230413_Ka_scan43.fits"
        g = SDFITSLoad(f)

        expect = len(g.udata("IFNUM")) * len(g.udata("PLNUM")) * len(g.udata("FDNUM")) / len(g.udata("CAL"))
        mask = g._column_mask(
            {"OBJECT": "1256-0547", "DATE-OBS": "2023-04-13T04:10:30.25", "TUNIT7": "Counts", "CAL": "T"}
        )
        assert mask.sum() == expect
        assert list(mask[0:4]) == [False, True, False, True]
        _mask = (
            (g["OBJECT"] == "1256-0547")
            & (g["DATE-OBS"] == "2023-04-13T04:10:30.25")
            & (g["TUNIT7"] == "Counts")
            & (g["CAL"] == "T")
        )
        assert np.all(mask == _mask)
