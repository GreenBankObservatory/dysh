import warnings
from unittest.mock import patch

import astropy.units as u
import numpy as np
import pytest
from astropy.io import fits

from dysh.fits.gbtfitsload import GBTFITSLoad
from dysh.spectra.spectrum import IGNORE_ON_COPY, Spectrum, average_spectra
from dysh.util import get_project_testdata


def fit_gauss(spectrum):
    """
    Fit a Gaussian.
    """

    from astropy.modeling import models
    from specutils.fitting import fit_lines

    g_init = models.Gaussian1D(
        amplitude=spectrum.flux.max(), mean=spectrum.spectral_axis.mean(), stddev=spectrum.meta["FREQRES"] * u.Hz
    )
    g_fit = fit_lines(spectrum, g_init)

    return g_fit


def compare_spectrum(one, other, ignore_history=False, ignore_comments=False):
    """ """

    for k, v in vars(one).items():
        if k in IGNORE_ON_COPY:
            continue
        if ignore_history and k == "_history":
            continue
        if ignore_history and k == "_comments":
            continue
        elif k in ["_wcs"]:
            v.to_header() == vars(other)[k].to_header()
        elif k in ["_spectral_axis"]:
            for k_, v_ in vars(v).items():
                assert v_ == vars(vars(other)[k])[k_]
        else:
            assert v == vars(other)[k]


def loadfits(fname):
    hdu = fits.open(fname)
    spec = hdu[1].data["DATA"][0]
    hdu.close()
    return spec


class TestSpectrum:
    def setup_method(self):
        data_dir = get_project_testdata() / "AGBT05B_047_01"
        sdf_file = data_dir / "AGBT05B_047_01.raw.acs"
        sdf = GBTFITSLoad(sdf_file)
        self.getps0 = sdf.getps(scan=51, plnum=0)
        self.ps0 = self.getps0.timeaverage()
        self.getps1 = sdf.getps(scan=51, plnum=1)
        self.ps1 = self.getps1.timeaverage()
        self.ss = self.ps0._copy()  # Synthetic one.
        x = np.arange(0, len(self.ss.data))
        fwhm = 5
        stdd = fwhm / 2.35482
        mean = int(x.mean())
        self.ss._data = 1 * np.exp(-0.5 * (x - mean) ** 2 / stdd**2)
        self.ss.meta["FREQRES"] = abs(self.ss.meta["CDELT1"])
        self.ss.meta["FWHM"] = fwhm
        self.ss.meta["CENTER"] = self.ss.spectral_axis[mean].value
        self.ss.meta["STDD"] = stdd

    def test_add(self):
        """Test that we can add two `Spectrum`."""
        addition = self.ps0 + self.ps1

        assert addition.meta["EXPOSURE"] == (self.ps0.meta["EXPOSURE"] + self.ps1.meta["EXPOSURE"])
        assert np.all(addition.flux.value == (self.ps0.flux.value + self.ps1.flux.value))
        assert addition.flux.unit == self.ps0.flux.unit
        assert addition.velocity_frame == self.ps0.velocity_frame
        compare_spectrum(self.ps0, addition)

    def test_add_scalar(self):
        """Test that we can add a scalar to a `Spectrum`."""
        addition = self.ps0 + 10.0

        assert addition.meta["EXPOSURE"] == self.ps0.meta["EXPOSURE"]
        assert np.all(addition.flux.value == (self.ps0.flux.value + 10.0))
        assert addition.flux.unit == self.ps0.flux.unit
        assert addition.velocity_frame == self.ps0.velocity_frame
        compare_spectrum(self.ps0, addition)

    def test_radd_scalar(self):
        """Test that we can add a scalar to a `Spectrum`."""
        addition = 10 + self.ps0

        assert addition.meta["EXPOSURE"] == self.ps0.meta["EXPOSURE"]
        assert np.all(addition.flux.value == (self.ps0.flux.value + 10.0))
        assert addition.flux.unit == self.ps0.flux.unit
        assert addition.velocity_frame == self.ps0.velocity_frame
        compare_spectrum(self.ps0, addition)

    def test_sub(self):
        """Test that we can subtract two `Spectrum`."""
        subtraction = self.ps0 - self.ps1

        assert subtraction.meta["EXPOSURE"] == (self.ps0.meta["EXPOSURE"] + self.ps1.meta["EXPOSURE"])
        assert np.all(subtraction.flux.value == (self.ps0.flux.value - self.ps1.flux.value))
        assert subtraction.flux.unit == self.ps0.flux.unit
        assert subtraction.velocity_frame == self.ps0.velocity_frame
        compare_spectrum(self.ps0, subtraction)

    def test_sub_scalar(self):
        """Test that we can subtract a scalar from a `Spectrum`."""
        subtraction = self.ps0 - 10.0

        assert subtraction.meta["EXPOSURE"] == self.ps0.meta["EXPOSURE"]
        assert np.all(subtraction.flux.value == (self.ps0.flux.value - 10.0))
        assert subtraction.flux.unit == self.ps0.flux.unit
        assert subtraction.velocity_frame == self.ps0.velocity_frame
        compare_spectrum(self.ps0, subtraction)

    def test_rsub_scalar(self):
        """Test that we can subtract a scalar from a `Spectrum`."""
        subtraction = 10.0 - self.ps0

        assert subtraction.meta["EXPOSURE"] == self.ps0.meta["EXPOSURE"]
        assert np.all(subtraction.flux.value == (10.0 - self.ps0.flux.value))
        assert subtraction.flux.unit == self.ps0.flux.unit
        assert subtraction.velocity_frame == self.ps0.velocity_frame
        compare_spectrum(self.ps0, subtraction)

    def test_mul(self):
        """Test that we can multiply two `Spectrum`."""
        multiplication = self.ps0 * self.ps1

        assert np.all(multiplication.flux.value == (self.ps0.flux.value * self.ps1.flux.value))
        assert multiplication.flux.unit == self.ps0.flux.unit * self.ps1.flux.unit
        assert multiplication.velocity_frame == self.ps0.velocity_frame
        compare_spectrum(self.ps0, multiplication)

    def test_mul_scalar(self):
        """Test that we can multiply a `Spectrum` and a scalar."""
        multiplication = self.ps0 * 1.0

        assert np.all(multiplication.flux.value == (self.ps0.flux.value))
        assert multiplication.flux.unit == self.ps0.flux.unit
        assert multiplication.velocity_frame == self.ps0.velocity_frame
        compare_spectrum(self.ps0, multiplication)

    def test_rmul_scalar(self):
        """Test that we can multiply a `Spectrum` and a scalar."""
        multiplication = 1.0 * self.ps0

        assert np.all(multiplication.flux.value == (self.ps0.flux.value))
        assert multiplication.flux.unit == self.ps0.flux.unit
        assert multiplication.velocity_frame == self.ps0.velocity_frame
        compare_spectrum(self.ps0, multiplication)

    def test_div(self):
        """Test that we can divide two `Spectrum`."""
        division = self.ps0 / self.ps1

        assert np.all(division.flux.value == (self.ps0.flux.value / self.ps1.flux.value))
        assert division.flux.unit == self.ps0.flux.unit / self.ps1.flux.unit
        compare_spectrum(self.ps0, division)

    def test_div_scalar(self):
        """Test that we can divide a `Spectrum` by a scalar."""
        division = self.ps0 / 1.0

        assert np.all(division.flux.value == (self.ps0.flux.value))
        assert division.flux.unit == self.ps0.flux.unit
        assert division.velocity_frame == self.ps0.velocity_frame
        compare_spectrum(self.ps0, division)

    def test_write_read_fits(self, tmp_path):
        """Test that we can read fits files written by dysh"""
        s = self.ps1
        o = tmp_path / "sub"
        o.mkdir()
        file = o / "test_spectrum_write.fits"
        s.write(file, format="fits", overwrite=True)
        s2 = Spectrum.read(file, format="fits")
        assert np.all(s.data == s2.data)
        assert s.target == s2.target
        assert np.all(s2.spectral_axis == s.spectral_axis)
        # This test will generally fail because SITELONG, SITELAT, SITEELEV
        # don't have enough precision to match exactly our known GBT coordinates.
        # @todo make a close_enough comparison by differencing the observer
        # attributes
        # if s2.observer is not None:
        #    assert s.observer == s2.observer

    def test_write_read_ascii(self, tmp_path):
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
        o = tmp_path / "sub"
        o.mkdir()
        for f in fmt:
            file = o / f"testwrite.{f}"
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

    def test_history_and_comments(self):
        s = self.ps1
        s.baseline(2, remove=True)
        print(s.comments)
        s.add_comment("I removed a baseline")
        # This tests that the baseline command self-logged to history
        # AND that order is preserved because
        # the baseline history should be the last entry
        # print(s.history)
        assert "baseline(2,remove=True,)" in s.history[-1]
        print(s.comments)
        assert "I removed a baseline" in s.comments

    @patch("dysh.plot.specplot.plt.show")
    def test_slice(self, mock_show, tmp_path):
        """
        Test that we can slice a `Spectrum` using channels or units.
        For units we only consider frequencies for now.
        """
        meta_ignore = ["CRPIX1", "CRVAL1"]
        spec_pars = ["_target", "_velocity_frame", "_observer", "_obstime"]
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
        trimmed.plot(xaxis_unit="km/s", yaxis_unit="mK", vel_frame="itrs")
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
        trimmed_nu = self.ps0[spec_ax[s.start].to("Hz") : spec_ax[s.stop].to("Hz")]
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
        trimmed_vel = self.ps0[spec_ax[s.start] : spec_ax[s.stop]]
        assert np.all(trimmed_vel.flux == self.ps0.flux[s])
        assert np.all(trimmed_vel.spectral_axis.value - self.ps0.spectral_axis[s].value < 1e-5)
        for k, v in self.ps0.meta.items():
            if k not in meta_ignore:
                assert trimmed_vel.meta[k] == v
        for k in spec_pars:
            assert vars(trimmed_vel)[k] == vars(self.ps0)[k]
        trimmed_vel.plot(xaxis_unit="MHz", yaxis_unit="mK")

        # m.
        spec_ax = self.ps0.spectral_axis.to("m")
        trimmed_wav = self.ps0[spec_ax[s.start] : spec_ax[s.stop]]
        assert np.all(trimmed_wav.flux == self.ps0.flux[s])
        assert np.all(trimmed_wav.spectral_axis.value - self.ps0.spectral_axis[s].value < 1e-5)

    def test_smooth(self):
        """Test for smooth with `decimate=0`"""
        width = 10
        ss = self.ps0.smooth("gauss", width)
        assert ss.meta["CDELT1"] == self.ps0.meta["CDELT1"] * width
        assert ss.meta["FREQRES"] == pytest.approx(abs(self.ps0.meta["CDELT1"]) * width)
        assert np.diff(ss.spectral_axis).mean().value == ss.meta["CDELT1"]
        assert ss._resolution == pytest.approx(1)
        assert ss.velocity_frame == self.ps0.velocity_frame
        assert ss.doppler_convention == self.ps0.doppler_convention
        assert ss.observer.frame_attributes == self.ps0.observer.frame_attributes

        # Now, change the reference frame and see if it still works.
        from dysh.coordinates import astropy_frame_dict

        s = Spectrum.fake_spectrum()
        for frame in astropy_frame_dict.keys():
            try:
                s.set_frame(frame)
            except Exception:
                print(f"set_frame fails for: {frame}")
                continue
            print(f"frame set to: {frame}")
            s.smooth("box", 3)

    def test_smooth_decimate(self):
        """Test for smooth with `decimate!=width`."""
        width = 10
        decimate = 8
        ss = self.ps0.smooth("gauss", width, decimate)
        assert ss.meta["CDELT1"] == self.ps0.meta["CDELT1"] * decimate
        assert ss.meta["FREQRES"] == pytest.approx(abs(self.ps0.meta["CDELT1"]) * width)
        assert np.diff(ss.spectral_axis).mean().value == ss.meta["CDELT1"]
        assert ss._resolution == pytest.approx(width / decimate)

        # Now with synthetic data.
        sss = self.ss.smooth("gauss", width, decimate)
        assert sss.meta["CDELT1"] == self.ss.meta["CDELT1"] * decimate
        assert sss.meta["FREQRES"] == pytest.approx(abs(self.ss.meta["CDELT1"]) * width, abs=100)
        assert np.diff(sss.spectral_axis).mean().value == sss.meta["CDELT1"]
        assert sss._resolution == pytest.approx(width / decimate, abs=1e-2)
        # Also check the line properties.
        g_fit = fit_gauss(sss)
        fwhm = g_fit.stddev.value * 2.35482
        assert g_fit.mean.value == pytest.approx(self.ss.meta["CENTER"])
        assert np.sqrt(fwhm**2 - sss.meta["FREQRES"] ** 2) == pytest.approx(
            abs(self.ss.meta["CDELT1"]) * self.ss.meta["FWHM"], abs=abs(self.ss.meta["CDELT1"]) / 9.0
        )
        assert ss.velocity_frame == self.ps0.velocity_frame
        assert ss.doppler_convention == self.ps0.doppler_convention
        assert ss.observer.frame_attributes == self.ps0.observer.frame_attributes

    def test_smooth_nodecimate(self):
        """Test for smooth without decimation."""
        width = 10
        decimate = -1
        ss = self.ps0.smooth("gauss", width, decimate)
        assert ss.meta["CDELT1"] == self.ps0.meta["CDELT1"]
        assert ss.meta["FREQRES"] == pytest.approx(abs(self.ps0.meta["CDELT1"]) * width)
        assert np.diff(ss.spectral_axis).mean().value == ss.meta["CDELT1"]
        assert ss._resolution == pytest.approx(width / abs(decimate))

    def test_smooth_multi(self):
        """Test for multiple passes of smooth."""
        widths = [10, 15, 15.1]
        decimate = -1

        # Check fitter first.
        g_fit = fit_gauss(self.ss)
        assert g_fit.stddev.value * 2.35482 == pytest.approx(
            abs(self.ss.meta["CDELT1"]) * self.ss.meta["FWHM"], g_fit.stddev.std * 2.35482
        )

        # Now smooth the same Spectrum multiple times.
        sss = self.ss._copy()
        for w in widths:
            sss = sss.smooth("gauss", w, decimate=decimate)
            g_fit = fit_gauss(sss)
            fwhm = g_fit.stddev.value * 2.35482
            assert sss.meta["FREQRES"] == pytest.approx(abs(self.ss.meta["CDELT1"]) * w)
            assert np.sqrt(fwhm**2 - sss.meta["FREQRES"] ** 2) == pytest.approx(
                abs(self.ss.meta["CDELT1"]) * self.ss.meta["FWHM"], abs=abs(self.ss.meta["CDELT1"]) / 9.0
            )
            assert g_fit.mean.value == pytest.approx(self.ss.meta["CENTER"])

    def test_smooth_and_slice(self):
        """Test for slicing after smoothing."""
        width = 10
        decimate = 8
        s = slice(-500 * u.km / u.s, 500 * u.km / u.s)
        sss = self.ss.smooth("gauss", width, decimate)
        ssss = sss[s]
        g_fit = fit_gauss(sss)
        fwhm = g_fit.stddev.value * 2.35482
        assert ssss.meta["CDELT1"] == self.ss.meta["CDELT1"] * decimate
        assert ssss.meta["FREQRES"] == pytest.approx(abs(self.ss.meta["CDELT1"]) * width, abs=100)
        assert np.diff(sss.spectral_axis).mean().value == sss.meta["CDELT1"]
        assert sss._resolution == pytest.approx(width / decimate, abs=1e-2)
        assert g_fit.mean.value == pytest.approx(self.ss.meta["CENTER"])
        assert np.sqrt(fwhm**2 - sss.meta["FREQRES"] ** 2) == pytest.approx(
            abs(self.ss.meta["CDELT1"]) * self.ss.meta["FWHM"], abs=abs(self.ss.meta["CDELT1"]) / 9.0
        )

    def test_shift(self):
        """Test the shift method against the results produced by GBTIDL"""

        # Prepare test data.
        filename = get_project_testdata() / "TGBT21A_501_11/TGBT21A_501_11_ifnum_0_int_0-2_getps_152_plnum_0.fits"
        sdf = GBTFITSLoad(filename)
        nchan = sdf["DATA"].shape[-1]
        spec = Spectrum.fake_spectrum(nchan=nchan, seed=1)
        spec.data[nchan // 2 - 5 : nchan // 2 + 6] = 10
        org_spec = spec._copy()
        # The next two lines were used to create the input for GBTIDL.
        # sdf["DATA"] = [spec.data]
        # sdf.write("shift_testdata.fits")

        # Apply method to be tested.
        shift = 5.5
        spec = spec.shift(shift)

        # Internal tests.
        assert np.all(np.isnan(spec[: int(np.round(shift))].data))

        # Load GBTIDL answer.
        with fits.open(get_project_testdata() / "gshift_box.fits") as hdu:
            table = hdu[1].data
        gbtidl = table["DATA"][0]

        # Compare.
        # Ignore the edge channels to avoid edge effects.
        diff = (spec.data - gbtidl)[10:-10]
        assert np.all(abs(diff) < 5e-4)

        assert spec.meta["CRPIX1"] == org_spec.meta["CRPIX1"] + shift
        assert spec.spectral_axis[0].to("Hz").value == (
            org_spec.spectral_axis[0].to("Hz").value - spec.meta["CDELT1"] * shift
        )

    def test_find_shift(self):
        """
        Test the find_shift method.
        * Test that the shift with respect to itself is zero.
        * Test that it can find an integer shift.
        * Test that it can find a fractional shift.
        * Test that it can find a shift in velocity units.
        * Test that it can find a shift in a different frame.
        """
        spec = Spectrum.fake_spectrum(seed=1)

        # Shift should be zero.
        assert spec.find_shift(spec) == pytest.approx(0)

        chan_wid = np.mean(np.diff(spec._spectral_axis))

        # Shift by one channel.
        spec2 = spec._copy()
        spec2._spectral_axis = spec2.spectral_axis.replicate(value=spec2.spectral_axis + chan_wid)
        assert spec.find_shift(spec2) == pytest.approx(-1)

        # Shift by one and a half channels.
        spec3 = spec._copy()
        spec3._spectral_axis = spec3.spectral_axis.replicate(value=spec3.spectral_axis + chan_wid * 1.5)
        assert spec.find_shift(spec3) == pytest.approx(-1.5)

        # Shift in velocity.
        spec4 = spec._copy()
        velo = spec4.spectral_axis.replicate(value=spec4.spectral_axis.to("km/s"))
        dvel = velo[1] - velo[0]
        spec4._spectral_axis = spec4.spectral_axis.replicate(value=velo + dvel * 0.5)
        assert spec.find_shift(spec4) == pytest.approx(-0.5)

        # Shift in a different frame.
        assert spec.find_shift(spec4, frame="lsrk") == pytest.approx(-0.5)

    def test_align_to(self):
        """
        Tests for align_to method.
        * Test that align_to itself does not change the spectrum.
        * Test that aligning to a spectrum with a integer shift preserves the data.
        * Test that aligning to a spectrum with a fractional shift preserves signal amplitude within errorbars.
        """

        spec = Spectrum.fake_spectrum(nchan=1024, seed=1)
        org_spec = spec._copy()

        # Align to itself.
        spec = spec.align_to(spec)
        compare_spectrum(spec, org_spec, ignore_history=True)
        assert np.all((spec - org_spec).data == 0)

        # Align to a shifted version.
        shift = 5
        spec = spec.shift(shift)
        assert np.all((spec.data[shift:] - org_spec.data[:-shift]) == 0.0)

        # Align to a shifted version with signal.
        fshift = 0.5
        spec = self.ss._copy()
        org_spec = spec._copy()
        spec = spec.shift(shift + fshift)
        # The amplitude of the signal will decrease because of the sampling.
        tol = np.sqrt(
            (1 - np.exp(-0.5 * (fshift) ** 2 / spec.meta["STDD"] ** 2)) ** 2.0
            + (np.nanstd(spec.data[: len(spec.data) - 50])) ** 2.0
        )
        assert spec.max().value == pytest.approx(org_spec.max().value, abs=3 * tol)

    def test_average_spectra(self):
        """
        Tests for average_spectra.
        Although not a class method of `Spectra` it is included here to reuse the setup method of the test.
        * Test that it does not crash.
        * Test that it does not crash whit alignment.
        """

        ps0_org = self.ps0._copy()
        ps1_org = self.ps1._copy()

        avg = average_spectra((self.ps0, self.ps1))
        avg2 = self.ps0.average((self.ps1))
        compare_spectrum(avg, avg2, ignore_history=True, ignore_comments=True)

        avg = average_spectra((self.ps0, self.ps1), align=True)
        compare_spectrum(ps0_org, self.ps0, ignore_history=True, ignore_comments=True)
        compare_spectrum(ps1_org, self.ps1, ignore_history=True, ignore_comments=True)

    def test_spectrum_with_frame(self):
        """Regression test for issue #401 to ensure Spectrum.with_frame functions as advertised.
        https://github.com/GreenBankObservatory/dysh/issues/401
        """
        spec = Spectrum.fake_spectrum()
        # Ensure that repeated changes of frame to the same frame do note
        # change after the first transform
        s1 = spec.with_frame("lsrk")
        s2 = s1.with_frame("lsrk")
        assert all(s1.spectral_axis == s2.spectral_axis)
        assert s2.velocity_frame == "lsrk"
        assert s2.meta["VELDEF"][4:] == "-LSR"

        # Test that topographic results in an Exception because
        # users must provide an ITRS coordinate instance in that case.
        # First change the spectrum frame to somthing else because
        # if it is already topo/itrs, then the error is circumvented
        spec2 = Spectrum.fake_spectrum()
        spec2.set_frame("gcrs")
        with pytest.raises(ValueError):
            spec2.set_frame("topo")

        # Setting a new frame to the old frame does NOT result in an
        # identical observer attribute on the resultant SpectralAxis.
        # See https://github.com/astropy/astropy/issues/17506
        # This test ensures that the difference remains small.

        location_diff = np.sqrt(
            (s1.observer.x - s2.observer.x) ** 2
            + (s1.observer.y - s2.observer.y) ** 2
            + (s1.observer.z - s2.observer.z) ** 2
        )
        velocity_diff = np.sqrt(
            (s1.observer.v_x - s2.observer.v_x) ** 2
            + (s1.observer.v_y - s2.observer.v_y) ** 2
            + (s1.observer.v_z - s2.observer.v_z) ** 2
        )
        assert location_diff < 1.0e-5 * u.m
        assert velocity_diff < 2e-8 * u.km / u.s

    def test_velocity_axis_to(self):
        """Regression test for issue 372 https://github.com/GreenBankObservatory/dysh/issues/372
        Calling velocity_axis_to should not change the object spectral axis.
        """
        spec = Spectrum.fake_spectrum()
        id1 = id(spec.spectral_axis)
        sa = spec.velocity_axis_to(toframe="lsrk")
        assert id(spec.spectral_axis) == id1
        assert id(sa) != id1
        # converting to identical frame should not change the values in the spectral axis
        sa = spec.velocity_axis_to(toframe="topo")
        assert all(sa == spec.spectral_axis)

    def test_baseline(self):
        """Test for comparing GBTIDL baseline to Dysh baselines"""

        def test_single_baseline(sdf, ex_reg, order, model, gbtidl_bmodel, exclude_region_upper_bounds=True):
            """For use with TestSpectrum.test_baseline()"""
            dysh_spec = sdf.getspec(0)
            temp_bmodel = np.copy(dysh_spec.data)
            dysh_spec.baseline(
                order, ex_reg, remove=True, model=model, exclude_region_upper_bounds=exclude_region_upper_bounds
            )
            dysh_bmodel = temp_bmodel - np.copy(dysh_spec.data)
            diff = np.sum(np.abs(dysh_bmodel - gbtidl_bmodel))
            assert diff < 1.5e-6

        data_dir = get_project_testdata() / "AGBT17A_404_01"
        sdf_file = data_dir / "AGBT17A_404_01_scan_19_prebaseline.fits"
        sdf = GBTFITSLoad(sdf_file)
        gbtidl_two_reg = loadfits(data_dir / "AGBT17A_404_01_scan_19_bmodel.fits")
        gbtidl_no_reg = loadfits(data_dir / "AGBT17A_404_01_scan_19_noregion_bmodel.fits")

        order = 3
        ex_reg = [(0, 99), (381, 449), (721, 819)]

        test_single_baseline(sdf, ex_reg, order, "chebyshev", gbtidl_two_reg)
        test_single_baseline(sdf, ex_reg, order, "legendre", gbtidl_two_reg)
        test_single_baseline(sdf, ex_reg, order, "hermite", gbtidl_two_reg)
        test_single_baseline(sdf, ex_reg, order, "polynomial", gbtidl_two_reg)

        ex_reg = None
        test_single_baseline(sdf, ex_reg, order, "chebyshev", gbtidl_no_reg)

        # Test that no warnings are issued.
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            self.ps0.baseline(model="poly", degree=2)

        # Test baseline removal with velocity units and known anwser.
        s = Spectrum.fake_spectrum()
        s.spectral_axis.to("km/s")
        s._spectral_axis = s.spectral_axis.to("km/s")
        # Add a polynomial.
        pcoef = np.array([1, 5, 10])
        s += np.poly1d(pcoef)(np.linspace(1, -1, len(s.spectral_axis)))
        s.baseline(degree=2, model="poly", remove=True)
        for pn in s._baseline_model.unitless_model.param_names:
            fit_val = getattr(s._baseline_model.unitless_model, pn).value
            in_val = pcoef[::-1][int(pn[-1])]
            assert (in_val - fit_val) / in_val < 0.1

        # Test with units. Nothing to compare to though.
        dysh_spec = sdf.getspec(0)
        kms = u.km / u.s
        ex_reg = [(0 * kms, 4200 * kms), (6000 * kms, 7000 * kms), (8800 * kms, 90000 * kms)]
        dysh_spec.baseline(order, ex_reg, model="chebyshev")
        ex_reg = [(1 * u.GHz, 1.38 * u.GHz), (1.388 * u.GHz, 1.392 * u.GHz), (1.4 * u.GHz, 2 * u.GHz)]
        dysh_spec.baseline(order, ex_reg, model="chebyshev")
        ex_reg = [21 * u.cm, 21.5 * u.cm]
        dysh_spec.baseline(order, ex_reg, model="chebyshev")

    def test_baseline_include(self):
        """
        Test for comparing GBTIDL baseline to dysh baselines
        using inclusion regions.
        """

        def test_single_baseline(sdf, in_reg, order, model, gbtidl_bmodel, exclude_region_upper_bounds=True):
            """For use with TestSpectrum.test_baseline()"""
            dysh_spec = sdf.getspec(0)
            temp_bmodel = np.copy(dysh_spec.data)
            dysh_spec.baseline(
                order, include=in_reg, remove=True, model=model, exclude_region_upper_bounds=exclude_region_upper_bounds
            )
            dysh_bmodel = temp_bmodel - np.copy(dysh_spec.data)
            diff = np.sum(np.abs(dysh_bmodel - gbtidl_bmodel))
            assert diff < 1.5e-6

        data_dir = get_project_testdata() / "AGBT17A_404_01"
        sdf_file = data_dir / "AGBT17A_404_01_scan_19_prebaseline.fits"
        sdf = GBTFITSLoad(sdf_file)
        gbtidl_two_reg = loadfits(data_dir / "AGBT17A_404_01_scan_19_bmodel.fits")
        gbtidl_no_reg = loadfits(data_dir / "AGBT17A_404_01_scan_19_noregion_bmodel.fits")

        order = 3
        in_reg = [(99, 381), (449, 721)]

        test_single_baseline(sdf, in_reg, order, "chebyshev", gbtidl_two_reg)
        test_single_baseline(sdf, in_reg, order, "legendre", gbtidl_two_reg)
        test_single_baseline(sdf, in_reg, order, "hermite", gbtidl_two_reg)

        # Test with units. Nothing to compare to though.
        dysh_spec = sdf.getspec(0)
        kms = u.km / u.s
        in_reg = [(4200 * kms, 6000 * kms), (7000 * kms, 8800 * kms)]
        dysh_spec.baseline(order, include=in_reg, model="chebyshev")
        in_reg = [(1.38 * u.GHz, 1.388 * u.GHz), (1.392 * u.GHz, 1.4 * u.GHz)]
        dysh_spec.baseline(order, include=in_reg, model="chebyshev")
        in_reg = [21 * u.cm, 21.5 * u.cm]
        dysh_spec.baseline(order, include=in_reg, model="chebyshev")
