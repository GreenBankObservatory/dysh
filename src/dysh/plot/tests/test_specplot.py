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
        tpplot = tp.plot(interactive=False)
        of = tmp_path / "test_savefig.png"
        tpplot.savefig(of)
        assert of.is_file()

        tpplot = tp.plot(show_header=False, interactive=False)
        of = tmp_path / "test_savefig_noheader.png"
        tpplot.savefig(of)
        assert of.is_file()
