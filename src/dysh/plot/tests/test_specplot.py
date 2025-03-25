"""Tests for specplot."""

from pathlib import Path
from unittest.mock import patch

from dysh.fits import GBTFITSLoad
from dysh.util import get_project_testdata


class TestSpecplot:
    """ """

    @patch("dysh.plot.specplot.plt.show")
    def test_savefig(self, mock_show, tmp_path):
        """
        Test that plots are saved.
        """
        p = get_project_testdata()
        sdf = GBTFITSLoad(p / "AGBT20B_014_03.raw.vegas/AGBT20B_014_03.raw.vegas.A.fits")
        tp = sdf.gettp(scan=6, plnum=0, ifnum=0, fdnum=0).timeaverage()
        tp.plot()
        of = tmp_path / "test_savefig.png"
        tp.savefig(of)
        assert of.is_file()

        tp.plot(show_header=False)
        of = tmp_path / "test_savefig_noheader.png"
        tp.savefig(of)
        assert of.is_file()
