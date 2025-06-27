import os

import numpy as np
import pytest
from astropy import units as u
from astropy.io import fits

from dysh import util
from dysh.spectra import core

LOCALDIR = os.path.dirname(os.path.realpath(__file__))


class TestMeanTsys:
    """
    Tests for `dysh.spectra.core.dcmeantsys` function.
    """

    def setup_method(self):
        self.data_dir = util.get_project_testdata()

    def test_tsys(self):
        expected = np.array([17.24000345, 17.17140405, 17.15663698])

        path_to_file = f"{self.data_dir}/TGBT21A_501_11"
        filename = "TGBT21A_501_11_ifnum_0_int_0-2.fits"
        sdf_file = f"{path_to_file}/{filename}"

        # Open and select data.
        hdu_sdf = fits.open(sdf_file)
        table = hdu_sdf[1].data
        table_pl0 = table[table["PLNUM"] == 0]
        table_pl0_off = table_pl0[table_pl0["SCAN"] == 153]
        tcal = table_pl0_off["TCAL"][0]
        tsys_dysh = np.empty(table_pl0_off["DATA"].shape[0] // 2, dtype=float)
        for i in range(len(tsys_dysh)):
            tsys_dysh[i] = core.mean_tsys(
                calon=table_pl0_off["DATA"][1::2][i],
                caloff=table_pl0_off["DATA"][0::2][i],
                tcal=tcal,
            )
        # Compare.
        assert tsys_dysh == pytest.approx(expected)

    def test_tsys2(self):
        path_to_file = f"{self.data_dir}/TGBT21A_501_11"
        filein = f"{path_to_file}/TGBT21A_501_11.raw.vegas.fits"
        gbtidl_file = f"{path_to_file}/TGBT21A_501_11_getps_scan_152_intnum_0_ifnum_0_plnum_0.fits"

        hdu = fits.open(filein)
        table = hdu[1].data
        mask = (table["SCAN"] == 153) & (table["IFNUM"] == 0) & (table["PLNUM"] == 0)
        mask_on = table[mask]["CAL"] == "T"
        mask_off = table[mask]["CAL"] == "F"
        table_on = table[mask][mask_on]
        table_off = table[mask][mask_off]
        nchan = table["DATA"].shape[1]  # noqa: F841
        tsys_dysh = core.mean_tsys(table_on["DATA"][0], table_off["DATA"][0], table_on["TCAL"][0])

        hdu = fits.open(gbtidl_file)
        gbtidl_table = hdu[1].data
        gbtidl_tsys = gbtidl_table["TSYS"]

        assert tsys_dysh == gbtidl_tsys

    def test_smooth(self):
        data0 = np.array([0, 0, 0, 1, 0, 0, 0])
        data0h = np.array([0, 0, 0.25, 0.5, 0.25, 0, 0])
        data0b = np.array([0, 0.2, 0.2, 0.2, 0.2, 0.2, 0])
        data0g = np.array([0, 1.3563e-5, 0.055554, 0.888865, 0.055554, 1.3563e-5, 0])
        data1h, meta = core.smooth(data0, "hanning")
        assert data1h == pytest.approx(data0h)
        data1b, meta = core.smooth(data0, "boxcar", 5)
        assert data1b == pytest.approx(data0b)
        data1g, meta = core.smooth(data0, "gaussian", 1 / 2.35482)

        assert data1g == pytest.approx(data0g)

    def test_decimate(self):
        a1 = np.arange(100)
        a2, _meta = core.decimate(a1, 2, None)
        assert len(a2) == 50
        a2, _meta = core.decimate(a1, 3, None)
        assert len(a2) == 34
        with pytest.raises(ValueError):
            # n must be an integer
            a2, _meta = core.decimate(a1, 3.1415, None)
        # now test with meta
        meta = {"NAXIS1": 100, "CDELT1": 1, "CRPIX1": 1, "CRVAL1": 0, "FREQRES": 0.5}
        for i in range(1, 6):
            a2, _meta = core.decimate(a1, i, meta)
            assert _meta["NAXIS1"] == len(a2)
            assert _meta["CDELT1"] == meta["CDELT1"] * float(i)
            assert _meta["FREQRES"] == meta["FREQRES"]  # should not change since no smoothing
            assert _meta["CRPIX1"] == 1.0 + (meta["CRPIX1"] - 1) / i + 0.5 * (i - 1) / i
            assert _meta["CRVAL1"] == meta["CRVAL1"] + 0.5 * (i - 1) * meta["CDELT1"]

    def test_tsys_weight(self):
        """Test that `dysh.spectra.core.tsys_weight` works."""

        # Keys are the expected values, the rest the inputs.
        pairs = {
            1: {"exposure": 1, "delta_freq": 1, "tsys": 1},
            1: {"exposure": 1 * u.s, "delta_freq": 1 * u.Hz, "tsys": 1 * u.K},  # noqa: F601
            1: {"exposure": 1, "delta_freq": -1, "tsys": 1},  # noqa: F601
            2: {"exposure": 1, "delta_freq": -1, "tsys": 0.5**0.5},
        }

        for k, v in pairs.items():
            assert k == pytest.approx(core.tsys_weight(v["exposure"], v["delta_freq"], v["tsys"]))
