from unittest.mock import patch

import astropy.units as u
import numpy as np

from dysh.fits.gbtfitsload import GBTFITSLoad
from dysh.spectra.spectrum import Spectrum
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
        assert addition.flux.unit == self.ps0.flux.unit
        assert addition.velocity_frame == self.ps0.velocity_frame

    def test_add_scalar(self):
        """Test that we can add a scalar to a `Spectrum`."""
        addition = self.ps0 + 10.0

        assert addition.meta["EXPOSURE"] == self.ps0.meta["EXPOSURE"]
        assert np.all(addition.flux.value == (self.ps0.flux.value + 10.0))
        assert addition.flux.unit == self.ps0.flux.unit
        assert addition.velocity_frame == self.ps0.velocity_frame

    def test_radd_scalar(self):
        """Test that we can add a scalar to a `Spectrum`."""
        addition = 10 + self.ps0

        assert addition.meta["EXPOSURE"] == self.ps0.meta["EXPOSURE"]
        assert np.all(addition.flux.value == (self.ps0.flux.value + 10.0))
        assert addition.flux.unit == self.ps0.flux.unit
        assert addition.velocity_frame == self.ps0.velocity_frame

    def test_sub(self):
        """Test that we can subtract two `Spectrum`."""
        subtraction = self.ps0 - self.ps1

        assert subtraction.meta["EXPOSURE"] == (self.ps0.meta["EXPOSURE"] + self.ps1.meta["EXPOSURE"])
        assert np.all(subtraction.flux.value == (self.ps0.flux.value - self.ps1.flux.value))
        assert subtraction.flux.unit == self.ps0.flux.unit
        assert subtraction.velocity_frame == self.ps0.velocity_frame

    def test_sub_scalar(self):
        """Test that we can subtract a scalar from a `Spectrum`."""
        subtraction = self.ps0 - 10.0

        assert subtraction.meta["EXPOSURE"] == self.ps0.meta["EXPOSURE"]
        assert np.all(subtraction.flux.value == (self.ps0.flux.value - 10.0))
        assert subtraction.flux.unit == self.ps0.flux.unit
        assert subtraction.velocity_frame == self.ps0.velocity_frame

    def test_rsub_scalar(self):
        """Test that we can subtract a scalar from a `Spectrum`."""
        subtraction = 10.0 - self.ps0

        assert subtraction.meta["EXPOSURE"] == self.ps0.meta["EXPOSURE"]
        assert np.all(subtraction.flux.value == (10.0 - self.ps0.flux.value))
        assert subtraction.flux.unit == self.ps0.flux.unit
        assert subtraction.velocity_frame == self.ps0.velocity_frame

    def test_mul(self):
        """Test that we can multiply two `Spectrum`."""
        multiplication = self.ps0 * self.ps1

        assert np.all(multiplication.flux.value == (self.ps0.flux.value * self.ps1.flux.value))
        assert multiplication.flux.unit == self.ps0.flux.unit * self.ps1.flux.unit
        assert multiplication.velocity_frame == self.ps0.velocity_frame

    def test_mul_scalar(self):
        """Test that we can multiply a `Spectrum` and a scalar."""
        multiplication = self.ps0 * 1.0

        assert np.all(multiplication.flux.value == (self.ps0.flux.value))
        assert multiplication.flux.unit == self.ps0.flux.unit
        assert multiplication.velocity_frame == self.ps0.velocity_frame

    def test_rmul_scalar(self):
        """Test that we can multiply a `Spectrum` and a scalar."""
        multiplication = 1.0 * self.ps0

        assert np.all(multiplication.flux.value == (self.ps0.flux.value))
        assert multiplication.flux.unit == self.ps0.flux.unit
        assert multiplication.velocity_frame == self.ps0.velocity_frame

    def test_div(self):
        """Test that we can divide two `Spectrum`."""
        division = self.ps0 / self.ps1

        assert np.all(division.flux.value == (self.ps0.flux.value / self.ps1.flux.value))
        assert division.flux.unit == self.ps0.flux.unit / self.ps1.flux.unit

    def test_div_scalar(self):
        """Test that we can divide a `Spectrum` by a scalar."""
        division = self.ps0 / 1.0

        assert np.all(division.flux.value == (self.ps0.flux.value))
        assert division.flux.unit == self.ps0.flux.unit
        assert division.velocity_frame == self.ps0.velocity_frame

    def test_write_read_fits(self):
        s = self.ps1
        s.write("test_spectrum_write.fits", format="fits", overwrite=True)
        s2 = Spectrum.read("test_spectrum_write.fits", format="fits")
        assert np.all(s.data == s2.data)
        assert s.target == s2.target
        assert np.all(s2.spectral_axis == s.spectral_axis)
        # This test will generally fail because SITELONG, SITELAT, SITEELEV
        # don't have enough precision to match exactly our known GBT coordinates.
        # @todo make a close_enough comparison by differencing the observer
        # attributes
        # if s2.observer is not None:
        #    assert s.observer == s2.observer

    def test_write_read_ascii(self):
        fmt = [
            "basic",
            "ascii.commented_header",
            "commented_header",
            "ascii.fixed_width",
            "ascii.ipac",
            "ipac",
            "votable",
            "ecsv",
            "mrt",
        ]
        s = self.ps1
        for f in fmt:
            file = f"testwrite.{f}"
            s.write(file, format=f, overwrite=True)
            # ECSV is the only ascii format that can
            # complete a roundtrip unscathed.
            # (See https://docs.astropy.org/en/latest/io/unified.html#table-io)
            if f == "ecsv":
                s2 = Spectrum.read(file, format=f)
                assert np.all(s.data == s2.data)
                assert np.all(s.spectral_axis == s2.spectral_axis)
                assert s.target == s2.target
        # Test reading in a GBTIDL ascii file
        gbtidl_file = get_project_testdata() / "gbtidl_spectra/onoff-L_gettp_156_intnum_0_LSR.ascii"
        s2 = Spectrum.read(gbtidl_file, format="gbtidl")
        assert s2.meta["SCAN"] == 156
        assert s2.meta["OBJECT"] == "NGC2782"
        # veldef can't be determined from header.
        assert s2.flux.unit == u.ct
        assert s2.flux[0].value == 3608710.0
        assert s2.spectral_axis.unit == u.GHz
        # Now try a gbtidl file that is gzipped and Ta units
        gbtidl_file = get_project_testdata() / "gbtidl_spectra/onoff-L_getps_152_OPTI-HEL.ascii.gz"
        s2 = Spectrum.read(gbtidl_file, format="gbtidl")
        assert s2.meta["SCAN"] == 152
        assert s2.meta["OBJECT"] == "NGC2415"
        assert s2.meta["VELDEF"] == "OPTI-HEL"
        assert s2.flux.unit == u.K
        assert s2.flux[0].value == -0.1042543
        assert s2.spectral_axis.unit == u.Unit("km/s")
        # @todo remove the temporary files.  This should be done in a teardown() method

    @patch("dysh.plot.specplot.plt.show")
    def test_slice(self, mock_show, tmp_path):
        """
        Test that we can slice a `Spectrum` using channels or units.
        For units we only consider frequencies for now.
        """
        meta_ignore = ["CRPIX1", "CRVAL1"]
        spec_pars = ["_target", "_velocity_frame", "_observer", "_obstime", "_observer_location"]
        s = slice(1000, 1100, 1)

        trimmed = self.ps0[s]
        assert trimmed.flux[0] == self.ps0.flux[s.start]
        assert trimmed.flux[-1] == self.ps0.flux[s.stop - 1]
        assert np.all(trimmed.flux == self.ps0.flux[s])
        # The slicing changes the values at the micro Hz level.
        assert np.all(trimmed.spectral_axis.value - self.ps0.spectral_axis[s].value < 1e-5)
        # Check meta values. The trimmed spectrum has an additional
        # key: 'original_wcs'.
        for k, v in self.ps0.meta.items():
            if k not in meta_ignore:
                assert trimmed.meta[k] == v
        # Check additional object properties.
        # Not all of them make sense, since their shapes will be different.
        for k in spec_pars:
            assert vars(trimmed)[k] == vars(self.ps0)[k]
        # Check that we can plot.
        trimmed.plot(xaxis_unit="km/s", yaxis_unit="mK")
        # Check that we can write.
        o = tmp_path / "sub"
        o.mkdir()
        out = o / "test_spec_slice_write.fits"
        trimmed.write(out, format="fits", overwrite=True)
        # Check that we can read it back.
        trimmed_read = Spectrum.read(out, format="fits")
        assert np.all(trimmed.flux == trimmed_read.flux)
        assert np.all(trimmed.spectral_axis == trimmed_read.spectral_axis)
        assert trimmed.target == trimmed_read.target

        # Now slice using units.
		# Hz.
        spec_ax = self.ps0.spectral_axis
        trimmed_nu = self.ps0[spec_ax[s.start].to("Hz"):spec_ax[s.stop].to("Hz")]
        assert np.all(trimmed_nu.flux == self.ps0.flux[s])
        assert np.all(trimmed_nu.spectral_axis.value - self.ps0.spectral_axis[s].value < 1e-5)
        for k, v in self.ps0.meta.items():
            if k not in meta_ignore:
                assert trimmed_nu.meta[k] == v
        for k in spec_pars:
            assert vars(trimmed_nu)[k] == vars(self.ps0)[k]
        trimmed_nu.plot(xaxis_unit="km/s", yaxis_unit="mK")

		# km/s.
        spec_ax = self.ps0.spectral_axis.to("km/s")
        trimmed_vel = self.ps0[spec_ax[s.start]:spec_ax[s.stop]]
        assert np.all(trimmed_vel.flux == self.ps0.flux[s])
        assert np.all(trimmed_vel.spectral_axis.value - self.ps0.spectral_axis[s].value < 1e-5)
        for k, v in self.ps0.meta.items():
            if k not in meta_ignore:
                assert trimmed_vel.meta[k] == v
        for k in spec_pars:
            assert vars(trimmed_vel)[k] == vars(self.ps0)[k]
        trimmed_vel.plot(xaxis_unit="MHz", yaxis_unit="mK")
