import numpy as np

from dysh.fits.gbtfitsload import GBTFITSLoad
from dysh.util import get_project_testdata


class TestSpectrum:
    def setup_method(self):
        data_dir = get_project_testdata() / "AGBT05B_047_01"
        sdf_file = data_dir / "AGBT05B_047_01.raw.acs"
        sdf = GBTFITSLoad(sdf_file)
        getps0 = sdf.getps(51, plnum=0)
        self.ps0 = getps0.timeaverage()
        getps1 = sdf.getps(51, plnum=1)
        self.ps1 = getps1.timeaverage()

    def test_add(self):
        """Test that we can add two `Spectrum`."""
        addition = self.ps0 + self.ps1

        assert addition.meta["EXPOSURE"] == (self.ps0.meta["EXPOSURE"] + self.ps1.meta["EXPOSURE"])
        assert np.all(addition.flux.value == (self.ps0.flux.value + self.ps1.flux.value))

    def test_add_scalar(self):
        """Test that we can add a scalar to a `Spectrum`."""
        addition = self.ps0 + 10.0

        assert addition.meta["EXPOSURE"] == self.ps0.meta["EXPOSURE"]
        assert np.all(addition.flux.value == (self.ps0.flux.value + 10.0))

    def test_sub(self):
        """Test that we can subtract two `Spectrum`."""
        subtraction = self.ps0 - self.ps1

        assert subtraction.meta["EXPOSURE"] == (self.ps0.meta["EXPOSURE"] + self.ps1.meta["EXPOSURE"])
        assert np.all(subtraction.flux.value == (self.ps0.flux.value - self.ps1.flux.value))

    def test_sub_scalar(self):
        """Test that we can subtract a scalar from a `Spectrum`."""
        subtraction = self.ps0 - 10.0

        assert subtraction.meta["EXPOSURE"] == self.ps0.meta["EXPOSURE"]
        assert np.all(subtraction.flux.value == (self.ps0.flux.value - 10.0))

    def test_mul(self):
        """Test that we can multiply two `Spectrum`."""
        multiplication = self.ps0 * self.ps1

        assert np.all(multiplication.flux.value == (self.ps0.flux.value * self.ps1.flux.value))

    def test_div(self):
        """Test that we can divide two `Spectrum`."""
        division = self.ps0 / self.ps1

        assert np.all(division.flux.value == (self.ps0.flux.value / self.ps1.flux.value))
