"""Tests for scanplot."""

import astropy.units as u
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

    def test_spectral_unit_kHz(self):
        """
        Test that the y axis range of the image does not extend past
        the number of channels in the spectra. See issue #901.
        """
        plot = self.scan_block_tp.plot(spectral_unit="kHz")
        assert plot.axes.get_ylim() == (self.scan_block_tp.nchan, 0)

        plot = self.scan_block_tp.plot(spectral_unit="Hz")
        assert plot.axes.get_ylim() == (self.scan_block_tp.nchan, 0)

    def test_spectral_unit(self):
        """
        Test that valid spectral units do not result in exceptions.
        """

        _ = self.scan_block_tp.plot(spectral_unit="kHz")
        _ = self.scan_block_tp.plot(spectral_unit="m")
        _ = self.scan_block_tp.plot(spectral_unit=u.cm)
        _ = self.scan_block_tp.plot(spectral_unit=u.meV)
        _ = self.scan_block_tp.plot(spectral_unit=u.km / u.s)
