import os

import numpy as np
import pytest
from astropy.io import fits

from dysh import util
from dysh.spectra import core

LOCALDIR = os.path.dirname(os.path.realpath(__file__))


class TestMeanTsys:
    """
    Tests for `dysh.spectra.core.dcmeantsys` function.
    """

    def setup_method(self):
        self.root_dir = util.get_project_root()
        self.data_dir = f"{self.root_dir}/testdata"

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
                calon=table_pl0_off["DATA"][1::2][i], caloff=table_pl0_off["DATA"][0::2][i], tcal=tcal
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
        nchan = table["DATA"].shape[1]
        tsys_dysh = core.mean_tsys(table_on["DATA"][0], table_off["DATA"][0], table_on["TCAL"][0])

        hdu = fits.open(gbtidl_file)
        gbtidl_table = hdu[1].data
        gbtidl_tsys = gbtidl_table["TSYS"]

        assert tsys_dysh == gbtidl_tsys
