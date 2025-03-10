"""Tests for specplot."""

from pathlib import Path

from dysh.fits import GBTFITSLoad
from dysh.util import get_project_testdata


class TestSpecplot:
    """ """

    def test_savefig(self, tmp_path):
        """
        Test that plots are saved.
        """
        p = get_project_testdata()
        sdf = GBTFITSLoad(p / "AGBT20B_014_03.raw.vegas/")
        tp = sdf.gettp(scan=6, plnum=0).timeaverage()
        tp.plot()
        of = tmp_path / "test_savefig.png"
        tp.savefig(of)
        assert of.is_file()
