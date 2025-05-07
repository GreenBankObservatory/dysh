""" """

import pytest

from dysh.fits import gb20mfitsload
from dysh.util import get_project_testdata


class TestGB20mFITSload:
    def setup_method(self):
        self.data_path = get_project_testdata() / "20m/Skynet_60476_DR21_118886_68343.cyb.fits"
        with pytest.warns(UserWarning):
            self.sdfits = gb20mfitsload.GB20MFITSLoad(self.data_path)

    def test_names(self):
        """
        Test basic filename
        """
        assert self.data_path == self.sdfits.filename

    def test_find_cal(self):
        """
        Test that the indices for finding the calibration data
        are correct.
        """
        self.sdfits.select(ifnum=0, plnum=0)
        self.sdfits._find_cal_idx()
        assert self.sdfits.cal_beg_idx == 10
        assert self.sdfits.cal_end_idx == 128

    def test_obsmode(self):
        """
        Test that the 'OBSMODE' column has been updated.
        """
        assert set(self.sdfits._index["OBSMODE"]) == set(("onoff:on", "onoff:off"))

    def test_getps(self):
        """
        Test that we can use `getps`.
        """
        spec = self.sdfits.getps(ifnum=0, plnum=0)  # noqa: F841
