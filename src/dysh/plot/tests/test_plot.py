from pathlib import Path

from dysh import test_data_path
from dysh.fits.gbtfitsload import GBTFITSLoad
from dysh.plot import plot_spectrum
from dysh.plot.canvas import SpectrumPlot


class TestPlot:
    """Test various plotting things"""

    def setup_method(self):
        self.filename = test_data_path() / "TGBT21A_501_11/TGBT21A_501_11.raw.vegas.fits"
        self.sdfits = GBTFITSLoad(self.filename)
        self.psscan = self.sdfits.getps(scan=152, ifnum=0, plnum=0)
        self.spectrum_to_plot = self.psscan.timeaverage(weights="tsys")
        self.outfiles = []

    def teardown_method(self):
        del self.filename
        del self.sdfits
        del self.psscan
        del self.spectrum_to_plot
        for o in self.outfiles:
            if o.exists():
                o.unlink()
        del self.outfiles

    def test_plot_spectrum_static(self):
        test_plot = plot_spectrum(self.spectrum_to_plot, display=False, interactive=False)
        assert isinstance(test_plot, SpectrumPlot)

    def test_plot_spectrum_static_savefig(self):
        """Plot a spectrum and save it to a plot"""
        # Plot the spectrum
        test_plot = plot_spectrum(self.spectrum_to_plot, display=False, interactive=False)
        # Define where to save the plot
        savepath = Path("test_plot.png")
        # Delete the file if it already exists
        if savepath.exists():
            savepath.unlink()
        # Save the path of the image to delete after testing
        self.outfiles.append(savepath)
        # Try to save the plot as an image
        test_plot.savefig(savepath)
        # Assert that the image exists
        assert savepath.exists()
