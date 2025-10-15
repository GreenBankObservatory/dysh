#!/usr/bin/env python3
"""
Created on Mon Feb 17 15:39:26 2025

@author: mpound
"""

import astropy.units as u
import numpy as np
import pytest
from astropy.time import Time

from dysh.util.weatherforecast import GBTWeatherForecast


class TestWeatherForecast:
    """Test the GBT Gain Correction functions"""

    def setup_method(self):
        self.gbtwf = GBTWeatherForecast(testmode=True)
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
        # print(f"{z=}")
        # print(self.tau0_interp - z)
        assert self.tau0_interp == pytest.approx(z, abs=1e-4)
        z = self.gbtwf.fetch(specval=self.freq, vartype="Tatm", mjd=self.mjd, coeffs=True)
        assert self.tatm_interp == pytest.approx(z, abs=1e-4)
        z = self.gbtwf.fetch(specval=self.freq, vartype="Tatm", mjd=self.mjd, coeffs=False)
        assert z.shape == (230, 3)
        idx = np.where(z[:, 1] == self.freq.value)
        assert self.tatm_interp == pytest.approx(z[idx], abs=1e-4)
        z = self.gbtwf.fetch(specval=self.freq, vartype="Opacity", coeffs=False)
        assert z.shape == (230, 3)
        idx = np.where(z[:, 1] == self.freq.value)
        assert self.tau0_interp == pytest.approx(z[idx], abs=1e-4)

        # test that astropy Time works and also that Opacity is the default
        t = Time(self.mjd, format="mjd", scale="utc")
        z = self.gbtwf.fetch(specval=self.freq, mjd=t, coeffs=True)
        assert self.tau0_interp == pytest.approx(z, abs=1e-4)

        # test that catching dates before 05-may-2004 works
        with pytest.raises(ValueError):
            self.gbtwf.fetch(specval=self.freq, vartype="Opacity", mjd=53000)

        # interpolation only allowed for Opacity and Tatm
        with pytest.raises(ValueError):
            self.gbtwf.fetch(vartype="Winds", mjd=self.mjd, coeffs=True)
