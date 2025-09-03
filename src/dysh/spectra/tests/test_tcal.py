"""
Test for TCal, including testing from GBTFITSLoad.gettcal.
"""

import numpy as np
import pytest

import dysh.util as util
from dysh.fits import gbtfitsload


class TestTCal:
    def setup_method(self):
        sdf_file = (
            util.get_project_testdata() / "AGBT04A_008_02/AGBT04A_008_02.raw.acs/AGBT04A_008_02.raw.acs.testrim.fits"
        )
        self.sdf = gbtfitsload.GBTFITSLoad(sdf_file)
        self.tcal = self.sdf.gettcal(scan=226, ifnum=0, plnum=0, fdnum=0, zenith_opacity=0.08)

    def test_name(self):
        assert self.tcal.name == "3C286"

    def test_snu(self):
        assert np.all(self.tcal.snu != 0)

    def test_plot(self):
        # Test that we can plot the TCal object.
        import matplotlib.pyplot as plt

        plt.ioff()
        self.tcal.plot()

    def test_get_tcal(self):
        assert self.tcal.get_tcal() == pytest.approx(20.48841200180447)

    def test_smooth(self):
        # By default, do not decimate.
        assert len(self.tcal.smooth("box", 16).data) == len(self.tcal.data)
        # Decimate if asked to.
        assert len(self.tcal.smooth("box", 16, decimate=0).data) == len(self.tcal.data) // 16

    def test_nchan(self):
        assert self.tcal.nchan == 2**13

    def test_invalid_method(self):
        with pytest.raises(TypeError):
            self.sdf.gettcal(scan=226, ifnum=0, plnum=0, fdnum=0, zenith_opacity=0.08, method="not real")
