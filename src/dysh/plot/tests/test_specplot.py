"""Tests for specplot."""

import matplotlib.pyplot as plt
import pytest

from dysh.fits import GBTFITSLoad
from dysh.util import get_project_testdata

# Disable interactive plotting.
plt.ioff()


class TestSpecplot:
    """ """

    def setup_method(self):
        """
        Set up a plotter.
        """

        p = get_project_testdata()

        self.sdf = GBTFITSLoad(p / "AGBT20B_014_03.raw.vegas/AGBT20B_014_03.raw.vegas.A6.fits")
        self.tp = self.sdf.gettp(scan=6, plnum=0, ifnum=0, fdnum=0).timeaverage()
        self.tpplot = self.tp.plot()

    def test_savefig(self, tmp_path):
        """
        Test that plots are saved.
        """

        of = tmp_path / "test_savefig.png"
        self.tpplot.savefig(of)
        assert of.is_file()

        # Plot without header.
        tpplot = self.tp.plot(show_header=False)
        of = tmp_path / "test_savefig_noheader.png"
        tpplot.savefig(of)
        assert of.is_file()

    def test_bline(self):
        tp = self.sdf.gettp(scan=6, plnum=0, ifnum=1, fdnum=0).timeaverage()
        tp_plt = tp.plot()  # noqa: F841
        tp.baseline(2, remove=False)
        tp.baseline(2, remove=True)

    def test_oshows(self):
        tp1 = self.sdf.gettp(scan=6, plnum=0, ifnum=1, fdnum=0).timeaverage()
        tp2 = self.sdf.gettp(scan=6, plnum=0, ifnum=2, fdnum=0).timeaverage()
        self.tpplot = self.tp.plot(oshow=tp1)
        self.tpplot.oshow(tp2, color="r", linestyle="--")

    def test_kwargs(self):
        self.tp.plot(vel_frame="lsrk", doppler_convention="radio")
        self.tp.plot(title="testing")


class TestScanplot:
    def setup_method(self):
        """
        Set up a plotter.
        """

        p = get_project_testdata()

        self.sdf = GBTFITSLoad(p / "AGBT20B_014_03.raw.vegas/AGBT20B_014_03.raw.vegas.A6.fits")
        self.tp = self.sdf.gettp(scan=6, plnum=0, ifnum=0, fdnum=0)
        self.tpplot = self.tp.plot()

    def test_savefig(self, tmp_path):
        """
        Test that waterfall plots are saved.
        """

        of = tmp_path / "test_savefig.png"
        self.tpplot.savefig(of)
        assert of.is_file()

    def test_axis2_limits(self):
        """
        Test that the frequency axis matches that of the Spectrum object.
        """

        units = self.tpplot._axis2.get_ylabel()[-4:-1]
        assert self.tpplot._spectrum.spectral_axis.quantity.to(units).value.min() == pytest.approx(
            min(self.tpplot._axis2.get_ylim())
        )
