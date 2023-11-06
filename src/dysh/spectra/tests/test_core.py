
import os
import pytest
import numpy as np

from astropy.io import fits

from dysh import util
from dysh.spectra import core


LOCALDIR = os.path.dirname(os.path.realpath(__file__))


class TestMeanTsys():
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
        table_pl0 = table[table["PLNUM"]==0]
        table_pl0_off = table_pl0[table_pl0["SCAN"]==153]
        tcal = table_pl0_off["TCAL"][0]
        tsys_dysh = np.empty(table_pl0_off["DATA"].shape[0]//2, dtype=float)
        for i in range(len(tsys_dysh)):
            tsys_dysh[i] = core.mean_tsys(calon=table_pl0_off["DATA"][1::2][i],  
                                         caloff=table_pl0_off["DATA"][0::2][i],
                                         tcal=tcal)
        # Compare.
        assert tsys_dysh == pytest.approx(expected)
