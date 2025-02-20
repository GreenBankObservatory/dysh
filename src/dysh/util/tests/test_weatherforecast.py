#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 15:39:26 2025

@author: mpound
"""
import astropy.units as u
import numpy as np
from astropy.time import Time

from dysh.util.weatherforecast import GBTWeatherForecast

print("TEST OPACITY COEFFS=TRUE")
freq = [25] * u.GHz
mjd = np.array([60722.0])
g = GBTWeatherForecast()
g.fetch(specval=freq, vartype="Opacity", mjd=mjd, coeffs=True)
print("TEST TATM COEFFS=TRUE")
g.fetch(specval=freq, vartype="Tatm", mjd=mjd, coeffs=True)

print("TEST TAMT COEFFS= FALSE")
g.fetch(specval=freq, vartype="Tatm", mjd=mjd, coeffs=True)
print("TEST OPACITY COEFFS= FALSE")
g.fetch(specval=freq, vartype="Opacity", coeffs=False)

times = times = ["2024-11-11T00:00:00.123456789", "2025-01-02T10:00:00"]
t = Time(times, format="isot", scale="utc")

print("test default is Opacity and astropy Time works")
g.fetch(specval=freq, mjd=t, coeffs=True)
