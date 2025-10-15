#!/usr/bin/env python3
# A test of weather/gain classes that should only be run on the GBO network
# -*- coding: utf-8 -*-
import astropy.units as u
import numpy as np
import pytest

from dysh.spectra import Spectrum
from dysh.util.gaincorrection import GBTGainCorrection
from dysh.util.weatherforecast import GBTForecastScriptInterface, GBTWeatherForecast


@pytest.mark.gbo_only
class TestWeatherForecastGBO:
    """Test the GBT Gain Correction functions on GBO Network"""

    def setup_method(self):
        self.gbtwf = GBTWeatherForecast()
        self.mjd = np.array([60722.0, 60723.0])
        self.freq = [25] * u.GHz
        # t = Time(mjd, format="mjd", scale="utc")
        # Values returned by the script to be compared with what we get from interpolation
        # /users/rmaddale/bin/getForecastValues -type Opacity -timeList 60722 60723  -freq 25
        # Opacity(25,60722) = 0.2750
        # Opacity(25,60723) = 0.0694
        self.tau0_interp = np.array([[6.07220000e04, 2.50000000e01, 0.275], [6.07230000e04, 2.50000000e01, 6.94e-02]])
        # /users/rmaddale/bin/getForecastValues -type Tatm -timeList 60722 60723 -freq 25
        # Tatm(25,60722) = 274.4376
        # Tatm(25,60723) = 258.4968
        self.tatm_interp = np.array(
            [[6.07220000e04, 2.50000000e01, 274.4376], [6.07230000e04, 2.50000000e01, 258.4968]]
        )

    def test_script_interface(self):
        z = self.gbtwf.fetch(specval=self.freq, vartype="Opacity", mjd=self.mjd, coeffs=True)
        print(f"opacity true {z=}")
        # print(self.tau0_interp - z)
        assert self.tau0_interp == pytest.approx(z, abs=1e-4)
        z = self.gbtwf.fetch(specval=self.freq, vartype="Tatm", mjd=self.mjd, coeffs=True)
        print(f"tatm true {z=}")
        assert self.tatm_interp == pytest.approx(z, abs=1e-4)
        z = self.gbtwf.fetch(specval=self.freq, vartype="Tatm", mjd=self.mjd, coeffs=False)
        print(f"tatm false {z=}")
        # assert z.shape == (230, 3)
        # idx = np.where(z[:, 1] == self.freq.value)
        assert self.tatm_interp == pytest.approx(z, abs=1e-4)
        z = self.gbtwf.fetch(specval=self.freq, vartype="Opacity", mjd=self.mjd, coeffs=False)
        print(f"opacity false {z=}")
        # assert z.shape == (230, 3)
        # idx = np.where(z[:, 1] == self.freq.value)
        assert self.tau0_interp == pytest.approx(z, abs=1e-4)
        # Now test without only one MJD
        z = self.gbtwf.fetch(specval=self.freq, vartype="Opacity", mjd=self.mjd[0], coeffs=False)
        print(f"opacity false 1D {z=}")
        z = self.gbtwf.fetch(specval=self.freq, vartype="Opacity", mjd=self.mjd[0], coeffs=True)
        print(f"opacity true 1D {z=}")

        # get a variable that does not require frequency
        z = self.gbtwf.fetch(vartype="Winds", mjd=self.mjd, coeffs=False)
        ans = np.array([[6.0722e04, 0.0000e00, 2.6800e00], [6.0723e04, 0.0000e00, 1.5120e01]])

        assert ans == pytest.approx(z, abs=1e-3)

        # interpolation only allowed for Opacity and Tatm
        with pytest.raises(ValueError):
            self.gbtwf.fetch(vartype="Winds", mjd=self.mjd, coeffs=True)

    def test_pr_issues(self):
        # test the pull request 487 issues noted by pedro
        g = GBTForecastScriptInterface()
        with pytest.raises(TypeError):
            # freq and mjd must be lists
            g(vartype="Opacity", mjd=60733.2916667, freq=[50], coeffs=False)
        with pytest.raises(TypeError):
            # freq and mjd must be lists
            g(vartype="Opacity", mjd=[60733.2916667], freq=50, coeffs=False)

        g(vartype="Opacity", mjd=[60733.2916667], freq=[50], coeffs=False)

        self.gbtwf.fetch(vartype="Opacity", mjd=60733.2916667, specval=50 * u.GHz, coeffs=False)
        self.gbtwf.fetch(vartype="Opacity", mjd=60733.2916667, specval=[50 * u.GHz], coeffs=False)
        self.gbtwf.fetch(vartype="Opacity", mjd=None, specval=[50 * u.GHz], coeffs=False)

        f = Spectrum.fake_spectrum()
        with pytest.raises(ValueError):
            # channel widths less than 100 kHz will be rejected.  Too close together because we must round
            # before calling the GBO script. And honestly, there will be no difference in return values for
            # nu and nu+100kHz so user should think about what they are doing.
            self.gbtwf.fetch(vartype="Opacity", mjd=60733.2916667, specval=f.spectral_axis * 10, coeffs=False)

        # 1 MHz channels and freq>2GHz
        f = Spectrum.fake_spectrum(cdelt1=1e6, crval1=100e9)

        ggc = GBTGainCorrection()
        ggc.zenith_opacity(specval=f.spectral_axis)  # vartype=Opacity, mjd=now, coeff=True

        ggc.get_weather(vartype="Opacity", specval=f.spectral_axis)  # default mjd=now

        z = self.gbtwf.fetch(vartype="Opacity", mjd=f.obstime, specval=f.spectral_axis, coeffs=False)
        assert z.shape == (1024, 3)
        # this should raise an error because mjd=[Time] is wrong.
        # One should use a time object with multiple dates, i.e. Time([mjdlist],format="mjd",scale="utc")
        with pytest.raises(TypeError):
            self.gbtwf.fetch(vartype="Opacity", mjd=[f.obstime], specval=f.spectral_axis, coeffs=False)

        # check coeffs less than 2 GHz
        ans = np.array(
            [
                [5.92553185e04, 9.00000000e-01, 9.06443144e-03],
                [5.92553185e04, 1.40000000e00, 9.06408900e-03],
                [5.92553185e04, 2.10000000e00, 8.73946199e-03],
                [5.92553185e04, 3.00000000e00, 6.80079334e-03],
                [5.92553185e04, 4.00000000e00, 6.38185744e-03],
            ]
        )

        z = ggc.zenith_opacity(
            specval=[3.0, 4.0, 1.4, 2.1, 0.9] * u.GHz, mjd=f.obstime
        )  # default coeff=True, vartype=Opacity
        assert ans == pytest.approx(z, abs=1e-4)

        # check that bash can handle a ridiculous argument string
        nchan = 50000
        bigf = Spectrum.fake_spectrum(cdelt1=1e6, crval1=30e9, nchan=nchan)
        z = self.gbtwf.fetch(vartype="Opacity", mjd=bigf.obstime, specval=bigf.spectral_axis, coeffs=False)
        assert z.shape == (nchan, 3)


if __name__ == "__main__":
    tw = TestWeatherForecastGBO()
    tw.setup_method()
    tw.test_script_interface()
    tw.test_pr_issues()
