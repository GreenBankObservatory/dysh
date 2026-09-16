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
        self.tcal = self.sdf.gettcal(scan=227, ref=226, ifnum=0, plnum=0, fdnum=0, zenith_opacity=0.08)
        self.tcal_l = self.sdf.gettcal(
            scan=227, ref=226, ifnum=0, plnum=0, fdnum=0, zenith_opacity=0.08, method="linear"
        )

    def test_name(self):
        assert self.tcal.name == "3C286"

    def test_snu(self):
        assert np.all(self.tcal.snu != 0)

    def test_plot(self):
        # Test that we can plot the TCal object.
        self.tcal.plot()

    def test_get_tcal(self):
        assert self.tcal.get_tcal() == pytest.approx(18.611146926879883)
        assert self.tcal_l.get_tcal() == pytest.approx(18.622)

    def test_smooth(self):
        # By default, do not decimate.
        assert len(self.tcal.smooth("box", 16).data) == len(self.tcal.data)
        # Decimate if asked to.
        assert len(self.tcal.smooth("box", 16, decimate=0).data) == len(self.tcal.data) // 16

    def test_nchan(self):
        assert self.tcal.nchan == 2**13

    def test_oshow(self):
        p = self.tcal.plot()
        p.oshow(self.tcal)

    def test_mathod_match(self):
        tcal = self.sdf.gettcal(scan=227, ref=226, ifnum=0, plnum=0, fdnum=0, zenith_opacity=0.08, method="Quad")
        assert tcal.get_tcal() == self.tcal.get_tcal()
        tcal = self.sdf.gettcal(scan=227, ref=226, ifnum=0, plnum=0, fdnum=0, zenith_opacity=0.08, method="Lin")
        assert tcal.get_tcal() == self.tcal_l.get_tcal()

    def test_invalid_method(self):
        with pytest.raises(ValueError):
            self.sdf.gettcal(scan=227, ref=226, ifnum=0, plnum=0, fdnum=0, zenith_opacity=0.08, method="Triple")
