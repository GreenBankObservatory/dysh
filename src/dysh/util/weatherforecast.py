#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 13:13:33 2025

@author: mpound
"""
import ast
import os
import re
import subprocess
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Union

import astropy.units as u
import numpy as np
from astropy.time import Time
from astropy.units.quantity import Quantity
from numpy.polynomial.polynomial import Polynomial
from pandas import DataFrame, concat

from ..log import logger

__all__ = ["BaseWeatherForecast", "GBTWeatherForecast"]


class BaseWeatherForecast(ABC):
    """A generic interface to get weather forecast values from an external source"""

    @abstractmethod
    def fetch(
        specval: Quantity, valueType: list = None, mjd: Union[Time, np.ndarray] = None, coeffs=None, **kwargs
    ) -> np.ndarray:
        pass


class GBTWeatherForecast(BaseWeatherForecast):
    def __init__(self, **kwargs):
        self._forecaster = GBTForecastScriptInterface(**kwargs)

    def fetch(
        self, specval: Quantity, vartype: str = "Opacity", mjd: Union[Time, float] = None, coeffs=None
    ) -> np.ndarray:
        frequency = specval.to(u.GHz, equivalencies=u.spectral()).value
        return self._forecaster(freq=frequency, vartype=vartype, mjd=mjd, coeffs=coeffs)


# ------------------------------------------------------------------------------
# Man page of CLEO script that fetches weather data at GBO
# Note: The format for any argument that is xxList is a space-separated
#       list of values
# ------------------------------------------------------------------------------
#  getForecastValues -- Returns either Opacities, Tsys_atm, Tatm, or Rest for a
#                   list of frequencies and MJD's or a list of polynomial
#                   coefficients that can then be used to generate values. Uses
#                   polynomial fits which are archived once per hour by the CLEO
#                   forecast utilities.  The results are generated much faster
#                   than when  using the full-blown cleo forecasts programs,
#                   though there is a small loss of accuracy over that provided
#                   by cleo forecasts. The archive goes back to 5 May 2004 and
#                   extends to up to 7 days into the future.
#
#  SYNOPSIS
#      getForecastValues [-h | -help]
#            [-type Opacity|Tatm|AtmTsys|Tsys|TotalTsys|Est|Rest|Trcvr|Tau0|Tatm0
#            |Winds|WindEffect|SurfaceEffect|TotalEffect] [-timeList listMJDs]
#            [-startMJD mjd] [-stopMJD mjd] [-incrMJD mjd] [-freqList listGHz]
#            [-startFreq GHz] [-stopFreq GHz] [-incrFreq GHz] [-elev elevation]
#            [-dir path] [-coeff] [-debug]
#
# DESCRIPTION
#
# If -coeff is NOT specified:
#       Returns a list of time and frequency-stamped values using the polynomial
#       coefficients found as subdirectories of ~rmaddale/Weather.  The values
#       in these files are calculated hourly by the CLEO forecasting package
#       from the average of the three local town forecasts.
#
#       The returned values will be interpolated between the nearest time stamps
#       in the archived coefficient data.  The output will be multiple lines,
#       each line will give the quantities 'type' (Opacity, Tsys, Tatm, or
#       Rest, ...), followed by the frequency and MJD of the calculated value
#       within parentheses, an equal sign, and then the calculated value.
#       For example:
#
#           Opacity(115.45,54553.8158912) = -9999
#           Opacity(22.2,54553.8158912)   = 0.185762649207
#           Opacity(30.0,54553.8158912)   = 0.113792336834
#           Opacity(45.0,54553.8158912)   = 0.270589812026
#
#       If a value cannot be calculated (either the users supplied a frequency or
#       time stamp that is not covered by the archive, or there is a hole > 4 hrs
#       in the archive for the desired time stamp), -9999 will be returned.
#
# If -coeff is specified:
#       Returns a list of coefficients so that the user can then calculate their
#       own values at any frequency.  The returned values will be  interpolated
#       between the nearest time stamps in the archived coefficient data.  The
#       output will be multiple lines, each line will give the quantities 'type'
#       (Opacity, Tsys, Tatm, or  Rest), followed by the MJD of the calculated
#       value within parentheses, an equal sign, and then a TCL list.  The list
#       will have three sublists, one for each of three frequency bands.  Each
#       sublist will itself contain a list of polynomial coefficients that are to
#       be used in that lists band, followed by thechi square of the fit.
#       For example:
#
#       # Frequency Bands: 6 - 22.2, 22.2 - 50, 67 - 116, in GHz
#       Opacity(55088.6778935) = {{0.1253202697 -0.04723900904 0.006761806394
#                 -0.0004079950568 9.012215005e-06} 6.053071176e-06
#                 {3.519785243 -0.3868765402 0.01625918622 -        debug = kwargs.get("debug", False)0.0003067177559
#                 2.209766347e-06} 0.000106007689 {2109.17448 -141.386304
#                 3.940260533 -0.05841611896 0.0004858081949 -2.14846586e-06
#                 3.947017435e-09} 7.022269956e-05}
#
#       If coefficients cannot be found (probably because time stamp is not
#       covered by the archive or there is a hole > 4 hrs in the archive),
#       coefficients and chi squares of -9999 will be returned.
#
# UNITS
#     Output units for the generated quantities are:
#         Kelvin              : System Temperature, Atmosphere Temperature, Rest,
#                               Receiver Temperature
#         Nepers/Air Mass     : Opacity
#
# OPTIONS
#       The user can use default options by not specifying any command-line
#       arguments or supply values to the following arguments:
#
#       -h or -help
#           Brings up this help page and all other options are ignored
#       -coeff
#           Whether coefficients, and not values, are to be         debug = kwargs.get("debug", False)returned.  The default
#           is to supply values.  Only available for -types of Opacity or Tatm
#       -typeList list
#           A list of the types of the values to be returned (Default: Opacity).
#               Opacity: the total zenith opacity
#               Tatm: the opacity-weighted (i.e., representative) temperature of
#                   the atmosphere and is the value that should be used when
#                   fitting traditional tipping curves.
#               AtmTsys:  the part of the system temperature due to just the
#                   atmosphere, doesn't include CMB, spillover, and electronics,
#                   and is calculated for the specified elevation.
#               Tsys: Deprecated -- the same as AtmTsys
#               TotalTsys: the above-described Tsys value, augmented by an
#                   estimate of the contributions from the receiver, spillover,
#                   and the CMB.
#               Est: the Effective System Temperature for the specified elevation.
#               Rest: the Relative Effective System Temperature, which includes
#                   the contributions from the CMB, spillover, and electronics,
#                   and is calculated for the specified elevation.  Essentially,
#                   the predicted loss in gain due to atmospheric opacity
#               Trcvr: An estimate of the receiver temperature for the given
#                   frequency and MJD
#               Tau0: The best possible opacity for the given frequency and MJD
#               Tau10, Tau25, Tau50, Tau75, Tau90: The opacity for various
#                   percentile weather conditions for the given frequency and MJD,
#                   based on multi-year statistical studies.
#               Tatm0: The Tatm for the given frequency and MJD at the time of the
#                   best possible opacity.
#               Tatm10, Tatm25, Tatm50, Tatm75, Tatm90: The Tatm for various
#                   percentile weather conditions for the given frequency and MJD,
#                   based on multi-year statistical studies.
#               Winds: The wind speed in MPH.  A specified freqList is ignored.
#               WindEffect: The predicted loss in point-surface efficiency due to
#                   winds
#               SurfaceEffect: The predicted loss in point-surface efficiency due
#                   to a deformed surface during PTCS daytime observing.
#              TotalEffect: The product of the Rest and the wind and surface
#                  effects.  That is, the predicted loss in gain due to all the
#                  various weather factors.
#              MinElev: Suggested minimum elevation for an object that rises to
#                  the elevation given by the value of -el        debug = kwargs.get("debug", False)ev.  Observing above
#                  the suggested elevation should keep the loss in gain due to
#                  atmospheric opacity to no more than 70% of the loss at transit
#                  (i.e., < factor of ~2 increase in observing time).
#      -elev deg
#          The elevation, in degrees, to use for Tsys, Est, Rest, TotalTsys,
#          TotalEffect, and MinElev calculations.  Ignored for all other
#          calculations.  Default is 30 deg.
#      -freqList list
#          List of frequencies in GHz that can range from 6 to 116 GHz.
#          Ignored if either a type of Winds or -coeff is specified
#      -startFreq, -stopFreq, incrFreq
#          Alternatively, one can enter a start, stop frequency and a frequency
#          increment, all in GHz.  These are ignored if freqList, Winds, or -coef
#          is specified.  Defaults are 6, 116, and 1 GHz
#      -timeList list
#          List of MJD's.  Default is the current time.  Must be between
#          5 May 2004 and within two hours of the last entry in the archive,
#          which is usually the current hour.
#      -startMJD, -stopMJD, incrMJD
#          Alternatively, one can enter a start, stop time and a time increment,
#          all specified as MJD's  These are ignored if timeList is specified.
#          Default is the current time and 0.04166666 (i.e., 1 hr).
#      -dir path
#          Optional path to the database files.  Normally, this is not needed
#          since the default path usually points to the best dataset.
#      -debug
#          Turns on debugging
#


class GBTForecastScriptInterface:
    def __init__(self, path: Union[Path, str] = "/users/rmaddale/bin/getForecastValues", **kwargs):
        debug = kwargs.get("debug", False)
        self._path = Path(path)
        self._fit = None
        # maximum number of fit polynomial coefficients as of 2/2025. Probably won't change
        self._MAX_COEFFICIENTS = 7
        # frequency ranges as of 2/2025.
        self.fr = [2.0, 22.0, 22.0 + 1e-9, 50.0, 67.0, 116.0]  # GHz
        ccols = [f"c{n}" for n in np.arange(self._MAX_COEFFICIENTS)]
        self._fitcols = ["MJD", "freqLoGHz", "freqHiGHz"] + ccols
        self._df = DataFrame(columns=self._fitcols)
        if not debug:
            if not self._path.exists() or not self._path.is_file():
                raise ValueError(f"{self._path} does not exist or is not a file")

    # For using the fit coefficients, we create a DataFrame that has columns
    # MJD freqLoGHz  freqHiGHz, coeff1, coeff2, ... coeffN
    # Where freqLo and freqHi are the frequency range over which the coefficients
    # are valid.
    # In order to have a regular grid, there will be as many coefficients as the ferquency
    # range that has the most coefficiets (currently 100GHz range with 7 coefficients),
    # and the high order coefficients for the other ranges will be set to zero.
    def __call__(
        self, coeffs: bool = True, vartype: str = "Opacity", freq: list = None, mjd: list = None
    ) -> np.ndarray:
        """set up and call the GBO weather script

        If polynomial coefficients `coeffs` are given then, the return values will be computed as a function of frequency:

            `:math:` \tau_0 = \sum_{i=0}^{n} C_i*\nu_i

        where `:math:`C_i are the coefficients and `:math:`\nu is the frequency **in GHz**.

        Parameters
        ----------
        coeffs : bool, optional
            Fetch the polynomial coefficients by passing '-coeffs' to the script.  The default is True.
        vartype : str optional
            Which weather variable to fetch. The default is "Opacity".
        freq : list, optional
           An input frequency list in GHz at which to evaluate the weather data. If `coeffs=True`, the polynomial
           is fetch first and then the `vartype` at each frequency is evaluated. The default is None.
        mjd : list, optional
           An input data list in MJD at which to evaluate the weather data. The default is None.

        Returns
        -------
        weather_data : np.ndarray
            The requested weather data evaluated at the input frequencies and MJDs.
        """
        #    Polynomial coefficients in order of increasing degree, including the constant term i.e.,
        #     ``(1, 2, 3)`` give ``1 + 2*x + 3*x**2``
        # if coeffs is not None:
        #    p = Polynomial(coeffs)
        print(f"{coeffs=}, {vartype=}, {freq=}, {mjd=}")
        _args = f"{self._path.as_posix()} "
        if coeffs:
            # call with -coeff
            _args += f"-coeff -type {vartype} "
            if mjd is not None:
                mjdformat = len(mjd) * "{:.2f} "
                timearg = f"-timeList {mjdformat.format(*mjd)}"
                _args += timearg

            script_output = self._call_script(_args)
            print(f"{script_output=}")
            self._df = concat(
                [self._df, DataFrame(data=self._parse_coefficients(script_output), columns=self._fitcols)],
                ignore_index=True,
            )
            if freq is None:
                raise ValueError(f"You must give a frequency list with {coeffs=}.")
            return self._eval_polynomial(freq, mjd)
        else:
            # call with other args and -type Opacity
            pass

    def _eval_polynomial(self, freq: list, mjd: list) -> np.ndarray:
        # freq is in GHz
        # returns array n_mjd x n_freq
        final = None
        for d in mjd:
            df = self._df[self._sdf.mjd == d]
            for f in freq:
                z = []
                df = df[(df.freqLoGHz <= f) & (df.freqHiGHz >= f)]
                coefficients = df.loc[:, df.columns.str.contains("^c")].to_numpy()
                p = Polynomial(coefficients)
                z.append(p(freq))
            ary = np.hstack([mjd, z])
            if final:
                final = np.vstack([final, ary])
            else:
                final = ary
        return final

    def _call_script(self, str_args: str) -> str:
        """call the script via python `subprocess` and return the output as a str. Lines will be separated by \n"""
        # thanks, Evan!
        print(f"Calling {str_args}")
        output = subprocess.run(str_args.split(), stdout=subprocess.PIPE).stdout
        return str(output.decode("utf-8"))

    def _parse_coefficients(self, script_output: str) -> np.ndarray:
        """Parse the coefficient list that comes out of `getForecastValues` script
        when `-coeff` is passed to it.
        """
        # a single line looks like this (including the # in front of Frequency)
        # # Frequency Bands: 2-22, 22-50, 67-116 GHz
        # Opacity(60719.88309) = {{a0 a1 a2 a3 a4} chisq_2-22} {{b0 b1 b2 b3 b4} 6.594244413740909e-6} {{c1 c2 c3 c4 c5 c6 c7} chisq_c}
        # where a_n are for 2-22GHz, b_n is for 22-50 GHz, and c_n are for 67-116 GHz
        # A bad fit will contain -9999 in the coefficients or values.
        #
        # If a date is given that is outside the range the script errors out with a long
        # error message that start's with 'can't read "Coeffs'
        # script_output.replace("{{", "\n{{")
        lines = script_output.split("\n")
        out = None
        for line in lines:
            if line.startswith("#") or line == "":
                continue
            row = self._parse_coeff_line(line)
            if out is None:
                out = row
            else:
                out = np.vstack([out, row])
        return out

    def _parse_coeff_line(self, line: str) -> tuple:
        """parse a single line of 'coefficient' script output
            Because the number of coefficients differs for the different
            frequency ranges (4 for low and mid frequency, 7 for high frequency range),
            this function will return a 3x7  np.ndarray  with zeros
            for the undefined low/mid frequency coefficients, which maybe I'll think about.

        Parameters
        ----------
        line : str
            a single line with `Opacity(MJD) = {{a0..an} chi_a {{b0..bn} chi_b {c0..cn} chi_c}`

        Returns
        -------
        coeffs - tuple
            The tuple consists of (mjd,coeff_list) where coeff_list is as described above.
        """
        # munge the string so it looks like a list of lists
        logger.debug(f"parsing ###{line=}###")
        line = re.sub(" +", " ", line)  # remove duplicate spaces
        x = line.split("=")
        mjd = float(x[0].replace("Opacity(", "").replace(")", "").strip())
        s = x[1].strip().replace("{", "[").replace("}", "]").replace(" ", ",")
        ary = ast.literal_eval(s)
        # now drop the chisq values which happen to not be contained in lists,
        # so easily found.
        # print(f"##{ary=}##")
        out = None
        p = []
        for a in ary:
            n = [q for q in a if isinstance(q, list)]
            p.append(n[0])  # get rid of double []
        # pad all out to maximum array length
        maxlen = np.max([len(z) for z in p])
        for z in p:
            z += [0] * (maxlen - len(z))
        # out = np.vstack([p])
        row = np.array(
            [
                np.hstack([mjd, self.fr[0], self.fr[1], p[0]]),
                np.hstack([mjd, self.fr[2], self.fr[3], p[1]]),
                np.hstack([mjd, self.fr[4], self.fr[5], p[2]]),
            ]
        )
        return row
        # return (mjd, out)
