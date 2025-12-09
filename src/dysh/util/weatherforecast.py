#!/usr/bin/env python3
"""
Created on Wed Feb 12 13:13:33 2025

@author: mpound
"""

import ast
import re
import subprocess
from abc import ABC, abstractmethod
from collections.abc import Sequence
from pathlib import Path

import astropy.units as u
import numpy as np
from astropy.time import Time
from astropy.units.quantity import Quantity
from numpy.polynomial.polynomial import Polynomial
from pandas import DataFrame

from ..log import logger
from .core import to_mjd_list, to_quantity_list

__all__ = ["BaseWeatherForecast", "GBTForecastScriptInterface", "GBTWeatherForecast"]


class BaseWeatherForecast(ABC):
    """A generic interface to get weather forecast values from an external source"""

    @abstractmethod
    def fetch(
        specval: Quantity,
        valueType: list = None,  # noqa: RUF013
        mjd: Time | np.ndarray = None,
        coeffs=None,
        **kwargs,  # noqa: RUF013, RUF100
    ) -> np.ndarray:
        pass


# @todo (maybe) Formally specval = None is allowed then the script returns data al all frequencies 2-116 GHz
class GBTWeatherForecast(BaseWeatherForecast):
    def __init__(self, **kwargs):
        self._testmode = kwargs.get("testmode", False)
        self._forecaster = GBTForecastScriptInterface(**kwargs)
        self.LOWER_MJD_LIMIT = 53130

    # this could just as easily be __call__
    def fetch(
        self,
        specval: Quantity = None,
        vartype: str = "Opacity",
        mjd: Time | float = None,
        coeffs=True,
    ) -> np.ndarray:
        r"""Call the GBO weather script and parse the results into numbers.
        For frequencies below 2 GHz, the value at 2 GHz will be returned since the `getForecastValues` does not
        cover < 2GHz.  Returned values will be sorted by frequency, low to high.

        Parameters
        ----------
        specval : `~astropy.units.quantity.Quantity`, optional
            The spectral value -- frequency or wavelength -- at which to compute `vartype`
            For data such as 'Winds' that don't depend on frequency, `specval` can be None.
        vartype : str, optional
            Which weather variable to fetch. See Notes for a description of valid values.
            **If the user is not on the GBO network , the only variable available is Opacity.**
        mjd : `~astropy.time.Time` or float
            The date at which to compute the opacity. If given as a float, it is interpreted as
            Modified Julian Day.  Default: None, meaning the data will be fetched at the most recent MJD available.
            If the user is not on the GBO network, this argument is ignored and the opacity will only be a function of frequency.
        coeffs : bool
            If True and at GBO, `getForecastValues` will be passed the `-coeffs` argument which returns
            polynomial coefficients to fit `vartype` as a function of frequency for each MJD.
            **This is only valid for `vartype` "Opacity" or "Tatm."**
            Because the polynomial is only defined from 2 GHz to 116 GHz, for values below 2 GHz the value for 2 GHz
            will be returned.

        Notes
        -----
        The `vartype` name for values that can be returned are described below. **These are case-sensitive.**
        This description comes from the help
        text for `getForecastValues` and may not fully describe the data that are returned, e.g., when other arguments
        are ignored.

        Opacity
            the total zenith opacity

        Tatm
            the opacity-weighted (i.e., representative) temperature of
            the atmosphere and is the value that should be used when
            fitting traditional tipping curves.

        AtmTsys
            the part of the system temperature due to just the
            atmosphere, doesn't include CMB, spillover, and electronics,
            and is calculated for the specified elevation.

        TotalTsys
            the above-described Tsys value, augmented by an
            estimate of the contributions from the receiver, spillover,
            and the CMB.

        Est
            the Effective System Temperature for the specified elevation.

        Rest
            the Relative Effective System Temperature, which includes
            the contributions from the CMB, spillover, and electronics,
            and is calculated for the specified elevation.  Essentially,
            the predicted loss in gain due to atmospheric opacity

        Trcvr
            An estimate of the receiver temperature for the given
            frequency and MJD

        Tau0
            The best possible opacity for the given frequency and MJD

        Tau10, Tau25, Tau50, Tau75, Tau90
            The opacity for various
            percentile weather conditions for the given frequency and MJD,
            based on multi-year statistical studies.

        Tatm0
            The atmospheric temperature `Tatm` for the given frequency and MJD at the time of the
            best possible opacity.

        Tatm10, Tatm25, Tatm50, Tatm75, Tatm90
            The `Tatm` for various
            percentile weather conditions for the given frequency and MJD,
            based on multi-year statistical studies.

        Winds
            The wind speed in MPH.  A specified `freqList` is ignored.

        WindEffect
            The predicted loss in point-surface efficiency due to
            winds

        SurfaceEffect
            The predicted loss in point-surface efficiency due
            to a deformed surface during PTCS daytime observing.

        TotalEffect
            The product of the Rest and the wind and surface
            effects.  That is, the predicted loss in gain due to all the
            various weather factors.

        MinElev
            Suggested minimum elevation for an object that rises to
            the elevation given by the value of `-elev`.  Observing above
            the suggested elevation should keep the loss in gain due to
            atmospheric opacity to no more than 70% of the loss at transit
            (i.e., < factor of ~2 increase in observing time).

        Examples
        --------

        Fetch the wind data for a range of dates.

        .. code:: python

            from dysh.util.weatherforecast import GBTWeatherForecast
            import numpy as np

            g = GBTWeatherForecast()
            wind = g.fetch(vartype="Winds", mjd=np.arange(60722,60732), coeffs=False)


        Returns
        -------
        weather_data : `~numpy.ndarray`
            The requested weather data evaluated at the input frequencies and MJDs.
            ** These will be sorted by frequency low to high **
        """
        frequency = None
        if specval is not None:
            specval = to_quantity_list(specval)
            frequency = specval.to(u.GHz, equivalencies=u.spectral()).value
        # If MJD was None, then it means the current MJD
        if mjd is None:
            mjd = Time.now()
        mjd_list = to_mjd_list(mjd)
        # catch dates before May 5, 2004, below which the weather archive does not exist.
        if mjd_list is not None:
            bad_date = [x < self.LOWER_MJD_LIMIT for x in mjd_list]
            if any(bad_date):
                raise ValueError("There is no weather date before MJD 53129/2004-May-05.")
        return self._forecaster(freq=frequency, vartype=vartype, mjd=mjd_list, coeffs=coeffs)


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
#                 {3.519785243 -0.3868765402 0.01625918622 -        debug = kwargs.get("testmode", False)0.0003067177559
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
#           Whether coefficients, and not values, are to be         debug = kwargs.get("testmode", False)returned.  The default
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
#                  the elevation given by the value of -el        debug = kwargs.get("testmode", False)ev.  Observing above
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
    """
    An interface to call the GBO weather forecast script.  Generally, users will not use this class directly, but
    rather use `~GBTWeatherForecast`.

    Parameters
    ----------
    path : str or `pathlib.Path`
        The script to run to get forecast values.
    debug : bool
        If True, don't check that `path` exists.  This is useful for testing when not on GBO network. Default: False
    """

    def __init__(self, path: Path | str = "/users/rmaddale/bin/getForecastValues", **kwargs):
        self._testmode = kwargs.get("testmode", False)
        self._path = Path(path)
        self._fit = None
        # maximum number of fit polynomial coefficients as of 2/2025. Probably won't change
        self._MAX_COEFFICIENTS = 7
        # frequency ranges as of 2/2025.
        self.fr = [2.0, 22.0, 22.0 + 1e-9, 50.0, 67.0, 116.0]  # GHz
        ccols = [f"c{n}" for n in np.arange(self._MAX_COEFFICIENTS)]
        self._fitcols = ["MJD", "freqLoGHz", "freqHiGHz"] + ccols  # noqa: RUF005
        self._valid_vartypes = [
            "Opacity",
            "Tatm",
            "AtmTsys",
            "TotalTsys",
            "Est",
            "Rest",
            "Trcvr",
            "Tau0",
            "Tau10",
            "Tau25",
            "Tau50",
            "Tau75",
            "Tau90",
            "Tatm0",
            "Tatm10",
            "Tatm25",
            "Tatm50",
            "Tatm75",
            "Tatm90",
            "Winds",
            "WindEffect",
            "SurfaceEffect",
            "TotalEffect",
            "MinElev",
        ]
        if not self._testmode:
            if not self._path.exists() or not self._path.is_file():
                raise ValueError(f"{self._path} does not exist or is not a file")

    # For using the fit coefficients, we create a DataFrame that has columns
    # MJD freqLoGHz  freqHiGHz, coeff1, coeff2, ... coeffN
    # Where freqLo and freqHi are the frequency range over which the coefficients
    # are valid.
    # In order to have a regular grid, there will be as many coefficients as the ferquency
    # range that has the most coefficiets (currently 100GHz range with 7 coefficients),
    # and the high order coefficients for the other ranges will be set to zero.
    # Note: input MJDs are rounded to the nearest ~5 minutes.  Data are only taken every hour
    def __call__(
        self,
        vartype: str = "Opacity",
        freq: list = None,  # noqa: RUF013
        mjd: list = None,  # noqa: RUF013
        coeffs: bool = True,  # noqa: RUF013, RUF100
    ) -> np.ndarray:
        r"""Call the GBO weather script and parse the results into numbers.

        Parameters
        ----------
        vartype : str, optional
            Which weather variable to fetch. See Notes for a description of valid values. The default is "Opacity".
        freq : list, optional
           An input frequency list in GHz at which to evaluate the weather data. If `coeffs=True`, the polynomial
           is fetched first and then the `vartype` at each frequency is evaluated. The default is None.
        mjd : list, optional
           An input data list in MJD at which to evaluate the weather data. The default is None.
        coeffs : bool, optional
            Fetch the polynomial coefficients by passing '-coeffs' to the script. **This is only valid for `vartype` "Opacity" or "Tatm."**
            The default is True.
            If polynomial coefficients are requested, the return values will be computed as a function of frequency:

                :math:`value = \sum_{i=0}^{n} C_i \nu^i`

            where :math:`C_i` are the coefficients and :math:`\nu` is the frequency **in GHz**.
            Because the polynomial is only defined from 2 GHz to 116 GHz, for values below 2 GHz the value for 2 GHz
            will be returned.

        Notes
        -----
        The `vartype` name for values that can be returned are described below. **These are case-sensitive.**
        This description comes from the help
        text for `getForecastValues` and may not fully describe the data that are returned, e.g., when other arguments
        are ignored.

        Opacity
            the total zenith opacity

        Tatm
            the opacity-weighted (i.e., representative) temperature of
            the atmosphere and is the value that should be used when
            fitting traditional tipping curves.

        AtmTsys
            the part of the system temperature due to just the
            atmosphere, doesn't include CMB, spillover, and electronics,
            and is calculated for the specified elevation.

        TotalTsys
            the above-described Tsys value, augmented by an
            estimate of the contributions from the receiver, spillover,
            and the CMB.

        Est
            the Effective System Temperature for the specified elevation.

        Rest
            the Relative Effective System Temperature, which includes
            the contributions from the CMB, spillover, and electronics,
            and is calculated for the specified elevation.  Essentially,
            the predicted loss in gain due to atmospheric opacity

        Trcvr
            An estimate of the receiver temperature for the given
            frequency and MJD

        Tau0
            The best possible opacity for the given frequency and MJD

        Tau10, Tau25, Tau50, Tau75, Tau90
            The opacity for various
            percentile weather conditions for the given frequency and MJD,
            based on multi-year statistical studies.

        Tatm0
            The atmospheric temperature `Tatm` for the given frequency and MJD at the time of the
            best possible opacity.

        Tatm10, Tatm25, Tatm50, Tatm75, Tatm90
            The `Tatm` for various
            percentile weather conditions for the given frequency and MJD,
            based on multi-year statistical studies.

        Winds
            The wind speed in MPH.  A specified `freqList` is ignored.

        WindEffect
            The predicted loss in point-surface efficiency due to
            winds

        SurfaceEffect
            The predicted loss in point-surface efficiency due
            to a deformed surface during PTCS daytime observing.

        TotalEffect
            The product of the Rest and the wind and surface
            effects.  That is, the predicted loss in gain due to all the
            various weather factors.

        MinElev
            Suggested minimum elevation for an object that rises to
            the elevation given by the value of `-elev`.  Observing above
            the suggested elevation should keep the loss in gain due to
            atmospheric opacity to no more than 70% of the loss at transit
            (i.e., < factor of ~2 increase in observing time).

        Examples
        --------

        Fetch the wind data for a range of dates.

        .. code:: python

            from dysh.util.weatherforecast import GBTForecastScriptInterface
            import numpy as np

            g = GBTForecastScriptInterface()
            winds = g(vartype="Winds", mjd=np.arange(60722,60732), coeffs=False)


        Returns
        -------
        weather_data : `~numpy.ndarray`
            The requested weather data evaluated at the input frequencies and MJDs.
            ** These will be sorted by frequency low to high **
        """
        if self._testmode:
            logger.warn("In debug mode, using test data.")
        self._check_vartype(vartype)
        logger.debug(f"{coeffs=}, {vartype=}, {freq=}, {mjd=}")
        _args = f"{self._path.as_posix()} -type {vartype} "
        if not isinstance(mjd, (Sequence, np.ndarray)):
            raise TypeError(f"mjd must be a list or numpy array, not {type(mjd)}")
        doctored_freq = None
        if freq is not None:
            if not isinstance(freq, (Sequence, np.ndarray)):
                raise TypeError(f"freq must be a list or numpy array, not {type(freq)}")
            # Ensure freq is a numpy array
            freq = np.array(freq)
            # sort the frequencies -- this is necessary because getForecastValues returns
            # the data sorted by frequency. So if we have to substitute low frequency values
            # the order must be guaranteed.
            freq.sort()
            self._check_deltafreq(freq)
            doctored_freq = np.copy(freq)
            # We have decided that if the frequecy is below 2 GHz, we will return the value at 2GHz.
            # Therefore we must send in a substitute list of frequencies, replacing anything below 2GHz
            # with 2GHz+epsilon, then replace that with the original list before returning the values.
            # The epsilon is necessary because the script removes duplicate frequencies.
            lo_freq_idx = np.where(doctored_freq < 2.0)
            lenlo = len(doctored_freq[lo_freq_idx])
            if lenlo != 0:
                # We have to add a tiny bit onto the 2.0 GHz because the script
                # will compress the return result, i.e. if freqList is 2 2 2 2 5 7,
                # the script will return only 3 values for 2,5,7.
                # So we have to trick it by making them within a few thousands of 2.
                a = np.round(np.random.rand(lenlo) * 0.001, 4)
                doctored_freq[lo_freq_idx] = 2.001 + a

        if mjd is not None:
            # round MJD to nearest 5 minutes. This helps to shorten the argument list so we don't run afoul of bash
            mjdformat = len(mjd) * "{:.4f} "
            timearg = f"-timeList {mjdformat.format(*mjd)} "
            _args += timearg
            n_mjd = len(mjd)
        else:
            n_mjd = 1

        if coeffs:
            if vartype != "Opacity" and vartype != "Tatm":
                raise ValueError("You can only use coeff=True for vartype Opacity or Tatm")  # limitation of the script
            # call with -coeff
            _args += "-coeff"

            if self._testmode:
                from .core import get_project_testdata

                script_file = get_project_testdata() / f"calibration/getForecastValues{vartype}_coeff.txt"
                logger.warn(f"Using testfile {script_file}")
                script_output = script_file.open("r").read()
            else:
                script_output = self._call_script(_args)
            # mjd needs float64
            self._df = DataFrame(
                data=self._parse_coefficients(vartype, script_output), columns=self._fitcols, dtype=float
            )
            if freq is None:
                raise ValueError(f"You must give a frequency list with {coeffs=}.")
            values = self._eval_polynomial(doctored_freq, mjd)
            # Now replace the original frequencies
            # We have the additional complication that if N MJDs were given,
            # the frequency array will be repeated N times. np.tile does this.
            if not self._testmode:
                # NB: The slice will fail with a ValueError array broadcast exception if the
                # any of the frequencies the user input happens to match any of
                # the doctored frequencies because the script removes duplicate frequencies.
                values[:, 1] = np.tile(freq, n_mjd)
        else:
            # call with other args and -type vartype  [-freqList -timeList]
            if doctored_freq is not None:
                # round freq to nearest 100 kHz (needed 4 digits for the low frequency fakeout)
                freqformat = len(doctored_freq) * "{:.4f} "
                freqarg = f"-freqList {freqformat.format(*doctored_freq)} "
                _args += freqarg

            if self._testmode:
                from .core import get_project_testdata

                script_file = get_project_testdata() / f"calibration/getForecastValues{vartype}.txt"
                logger.warn(f"Using testfile {script_file}")
                script_output = script_file.open("r").read()
            else:
                script_output = self._call_script(_args)
            logger.debug(f"{script_output=}")
            values = self._parse_list_values(vartype, script_output)
            # Now replace the original frequencies
            # We have the additional complication that if N MJDs were given,
            # the frequency array will be repeated N times. np.tile does this.
            if not self._testmode and freq is not None:
                # see caveat above
                values[:, 1] = np.tile(freq, n_mjd)

        # warn if any values returned are -9999 which
        # is what the script gives if it can't determine a value.
        if np.any(values == -9999.0):
            logger.warn(f"In fetching {vartype} a value of -9999 was detected. Be careful.")
        # remove any extra zero length dimensions added in parsing.
        return np.atleast_2d(np.squeeze(values))

    @property
    def valid_vartypes(self) -> list:
        "List of the valid weather variable type names that can be retrieved."
        return self._valid_vartypes

    def _check_vartype(self, vartype: str) -> None:
        if vartype is None or vartype not in self._valid_vartypes:
            raise ValueError(
                f"Unrecognized vartype {vartype}.  Valid choices are: {self._valid_vartypes} (case-sensitive)."
            )

    def _check_deltafreq(self, freq: np.ndarray) -> None:
        """When we pass a frequency list to script we round it to 100 kHz, so
        raise an error if CDELT is less than 100 kHz
        """
        if len(freq) < 2:
            return
        if np.abs(freq[0] - freq[1]) * u.GHz < 100 * u.kHz:
            raise ValueError("Frequencies in input list are too close, must be separated by at least 100 kHz")

    def _call_script(self, str_args: str) -> str:
        """call the script via python `subprocess` and return the output as a str. Lines will be separated by \n"""
        logger.debug(f"Calling {str_args}")
        # thanks, Evan!
        output = subprocess.run(str_args.split(), stdout=subprocess.PIPE).stdout
        return str(output.decode("utf-8"))

    def _parse_list_values(self, vartype: str, script_output: str) -> np.ndarray:
        """parse script output when values are returned instead of coefficients"""
        lines = script_output.split("\n")
        # tau = []
        out = None
        for line in lines:
            if line.startswith("#") or line == "":
                continue
            x = line.split("=")
            tau = float(x[1])
            vals = [float(q) for q in x[0].replace(f"{vartype}(", "").replace(")", "").strip().split(",")]
            row = np.array([vals[1], vals[0], tau])
            if out is None:
                out = row
            else:
                out = np.vstack([out, row])
        # because we will be slicing, ensure that out is 2D even if one dimension is zero.
        return np.atleast_2d(out)

    def _parse_coefficients(self, vartype: str, script_output: str) -> np.ndarray:
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
            row = self._parse_coeff_line(vartype, line)
            if out is None:
                out = row
            else:
                out = np.vstack([out, row])
        return out

    def _parse_coeff_line(self, vartype: str, line: str) -> tuple:
        """parse a single line of 'coefficient' script output
            Because the number of coefficients differs for the different
            frequency ranges (4 for low and mid frequency, 7 for high frequency range),
            this function will return a 3x7  np.ndarray  with zeros
            for the undefined low/mid frequency coefficients, which maybe I'll think about.

        Parameters
        ----------
        line : str
            a single line with e.g., `Opacity(MJD) = {{a0..an} chi_a {{b0..bn} chi_b {c0..cn} chi_c}`

        Returns
        -------
        coeffs : `~numpy.ndarray`
            A array of coefficient values for each frequency range.  MJD is the same in each row.
            [
                [mjd, lowfreq_0, hifreq_0, coefficient_0...coefficient_n]
                [mjd, lowfreq_1, hifreq_1, coefficient_0...coefficient_n]
                [mjd, lowfreq_2, hifreq_2, coefficient_0...coefficient_n]
             ]
        """
        # munge the string so it looks like a list of lists
        logger.debug(f"parsing ###{line=}###")
        line = re.sub(" +", " ", line)  # remove duplicate spaces
        x = line.split("=")
        mjd = float(x[0].replace(f"{vartype}(", "").replace(")", "").strip())
        s = x[1].strip().replace("{", "[").replace("}", "]").replace(" ", ",")
        ary = ast.literal_eval(s)
        # now drop the chisq values which happen to not be contained in lists,
        # so easily found.
        p = []
        for a in ary:
            n = [q for q in a if isinstance(q, list)]
            p.append(n[0])  # get rid of double []
        # pad all out to maximum array length
        maxlen = np.max([len(z) for z in p])
        for z in p:
            z += [0] * (maxlen - len(z))
        row = np.array(
            [
                np.hstack([mjd, self.fr[0], self.fr[1], p[0]]),
                np.hstack([mjd, self.fr[2], self.fr[3], p[1]]),
                np.hstack([mjd, self.fr[4], self.fr[5], p[2]]),
            ]
        )
        return row

    def _eval_polynomial(self, freq: list, mjd: list) -> np.ndarray:
        """Evaluate the polynomial at the given frequencies and MJDs

        Parameters
        ----------
        freq : list
            Frequencies in GHz.
        mjd : list
            Modified Julian Date.

        Returns
        -------
        polynomial : `~numpy.ndarray`
            Polynomial evaluated at `mjd` and `freq`.
        """
        # freq is in GHz
        # returns array n_mjd x n_freq
        # with values [mjd, freq, tau0]
        logger.debug(f"eval {freq=} {mjd=}")
        final = None
        if freq is None or mjd is None:
            raise ValueError("freq and mjd cannot be None")
        for d in mjd:
            df = self._df[np.round(self._df.MJD, 4) == np.round(d, 4)]
            for f in freq:
                z = []
                df = df[(df.freqLoGHz <= f) & (df.freqHiGHz >= f)]
                coefficients = df.loc[:, df.columns.str.contains("^c")].to_numpy()[0]
                p = Polynomial(coefficients)
                z.append(p(f))
                ary = np.hstack([d, f, z])
                if final is None:
                    final = ary
                else:
                    final = np.vstack([final, ary])
        # because we will be slicing, ensure that out is 2D even if one dimension is zero.
        return np.atleast_2d(final)
