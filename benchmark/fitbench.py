#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fit benchmark data
Created on Mon Apr 21 09:47:46 2025

@author: mpound
"""
from astropy.table import Table

# import astropy.units as units

import numpy as np
from scipy.optimize import curve_fit


def obj_func(params, a, b, c, offset):
    nfile, size, nflag = params
    return a * nfile + b * size + c * nflag + offset


def obj_func2(nfile, a, offset):
    return a * nfile + offset


x = None
y = None
for i in np.arange(1, 5):
    file = f"benchtest{i}.tab"
    t = Table.read(file, format="ascii.ecsv")
    df = t.to_pandas()
    # df1 = df[(df["# files"] == i) & (df["time"]>0)]
    df1 = df[(df["name"] != "start") & (df["name"] != "end")]
    load = df1["time"].to_numpy()  # ms
    size = df1["file size"].to_numpy()  # MN
    nfile = df1["# files"].to_numpy()
    nflag = df1["flag (lines)"].to_numpy()
    sf = df1.skipflags.to_numpy()
    sw = np.where(sf == "True")
    nflag[sw] = 0
    if x is None:
        x = np.array([load, nfile, size, nflag])
    else:
        x = np.hstack([x, [load, nfile, size, nflag]])

popt, pcov = curve_fit(obj_func, x[1:], x[0])

# popt, pcov = curve_fit(obj_func2, x[1], x[0])
