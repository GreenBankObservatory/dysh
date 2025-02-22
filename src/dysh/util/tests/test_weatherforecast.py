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
mjd = np.array([60722.0, 60723.0])
g = GBTWeatherForecast(testmode=True)
z = g.fetch(specval=freq, vartype="Opacity", mjd=mjd, coeffs=True)
print(z)
print("TEST TATM COEFFS=TRUE")
z = g.fetch(specval=freq, vartype="Tatm", mjd=mjd, coeffs=True)
print(z)

print("TEST TAMT COEFFS= FALSE")
g.fetch(specval=freq, vartype="Tatm", mjd=mjd, coeffs=False)
print("TEST OPACITY COEFFS= FALSE")
z = g.fetch(specval=freq, vartype="Opacity", coeffs=False)
print(z)
t = Time(mjd, format="mjd", scale="utc")

print("test default is Opacity and astropy Time works")
z = g.fetch(specval=freq, mjd=t, coeffs=True)
print(z)
