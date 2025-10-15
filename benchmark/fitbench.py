#!/usr/bin/env python3
"""
Fit benchmark data
Created on Mon Apr 21 09:47:46 2025

@author: mpound
"""

# import astropy.units as units
import argparse
import sys

import numpy as np
from astropy.table import Table, vstack
from scipy.optimize import curve_fit


def _oldobj_func(params, b, c, d, offset):
    size, nflag, nchan = params
    # return a * nfile +
    return b * size + c * nflag + d * nchan + offset


def obj_func1(params, a, offset):
    return a * params[0] + offset


def obj_func2(params, a, b, offset):
    return a * params[0] + b * params[1] * offset


def obj_func3(params, a, b, c, offset):
    return a * params[0] + b * params[1] + c * params[2] + offset


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog=sys.argv[0])
    parser.add_argument(
        "--files",
        "-f",
        nargs="+",
        type=str,
        required=True,
        help="The input astropy table(s) with the benchmark data. Multiple tables will be concatenated before fit.",
    )
    parser.add_argument(
        "--params",
        "-p",
        nargs="+",
        type=str,
        required=True,
        help="The parameters to use in the fit.  These must be column names in the table.",
    )
    parser.add_argument(
        "--depvar",
        "-d",
        type=str,
        required=False,
        default="time",
        help="The dependent variable to fit. Default is 'time'",
    )
    args = parser.parse_args()
    obj_func = [None, obj_func1, obj_func2, obj_func3]
    tab = None
    pset = set(args.params)
    for file in args.files:
        # print(file)
        t = Table.read(file, format="ascii.ecsv")
        cset = set(t.columns)
        if not pset.issubset(cset):
            raise ValueError(f"Did not find all parameters in table columns: {cset}")

        if tab is None:
            tab = t
        else:
            tab = vstack([tab, t])
    npar = len(args.params)
    nfunc = len(obj_func) - 1
    if npar > len(obj_func) - 1:
        raise ValueError(f"Can only handle {nfunc} parameters.  Add a new obj_func{npar} to fitbench.py")
    y = tab[args.depvar].data  # default time in ms
    pars = {}
    x = None
    for p in args.params:
        pars[p] = tab[p].data
        if x is None:
            x = pars[p]
        else:
            x = np.vstack([x, pars[p]])
    if False:
        tab.write("allbench.tab", overwrite=True, format="ascii.ecsv")

    print(obj_func[len(pars)])
    # print(f"{pars=}")
    print(f"{np.shape(x)=}")
    popt, pcov = curve_fit(obj_func[len(pars)], x, y)
    print(f"{popt=}\n{pcov=}\n{np.diag(pcov)=}")

    # popt, pcov = curve_fit(obj_func2, x[1], x[0])
