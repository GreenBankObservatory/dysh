"""Tests for specplot."""

from dysh.fits import GBTFITSLoad
from dysh.util import get_project_testdata


class TestSpecplot:
    """ """

    def test_savefig(self, tmp_path):
        """
        Test that plots are saved.
        """

        # Disable interactive plotting.
        import matplotlib.pyplot as plt

        plt.ioff()

        p = get_project_testdata()
        sdf = GBTFITSLoad(p / "AGBT20B_014_03.raw.vegas/AGBT20B_014_03.raw.vegas.A6.fits")
        tp = sdf.gettp(scan=6, plnum=0, ifnum=0, fdnum=0).timeaverage()
        tp.plot(show=False)
        of = tmp_path / "test_savefig.png"
        tp.savefig(of)
        assert of.is_file()

        tp.plot(show_header=False, show=False)
        of = tmp_path / "test_savefig_noheader.png"
        tp.savefig(of)
        assert of.is_file()
