"""Tests for scanplot."""

import numpy as np

from dysh.fits import GBTFITSLoad
from dysh.util import get_project_testdata


class TestScanPlot:
    def setup_method(self):
        """
        Set up a plotter.
        """

        fnm = get_project_testdata() / "AGBT18B_354_03/AGBT18B_354_03.raw.vegas/AGBT18B_354_03.raw.vegas.A.fits"

        self.sdfits = GBTFITSLoad(fnm)
        self.scan_block_tp = self.sdfits.gettp(scan=[6, 7], ifnum=2, plnum=0, fdnum=0)

    def test_set_clim(self):
        """
        Test that set_clim sets the limits properly,
        and updates the label. See issue #1082.
        """
        lims = (1e7, 1e8)
        plot = self.scan_block_tp.plot()
        plot.set_clim(*lims)

        exponent = np.floor(np.log10(max(lims)))
        expected_label = f"Counts($\\times10^{{{exponent:.0f}}}$)"
        assert plot._colorbar.ax.get_ylabel() == expected_label
        assert plot._colorbar.ax.get_ylim() == lims
