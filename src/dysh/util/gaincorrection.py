from abc import ABC, abstractmethod
from pathlib import Path
from typing import Union

import astropy.units as u
import numpy as np
from astropy.coordinates import Angle
from astropy.table import QTable
from astropy.time import Time
from astropy.units.quantity import Quantity

from ..log import logger
from ..util import get_project_configuration
from .core import to_mjd_list
from .weatherforecast import GBTWeatherForecast

__all__ = ["BaseGainCorrection", "GBTGainCorrection"]


class BaseGainCorrection(ABC):
    r"""This class is the base class for gain corrections. It is intended to be subclassed for
    specific antennas. Subclasses will be used to calculate various gain corrections to go
    from antenna temperature to other scales like brightness temperature of flux density.
    Subclasses can implement the following attributes:

    * `ap_eff_0`
        long wavelength aperture efficiency (number between 0 and 1), i.e., in the absence of surface errors,
        $\lambda >> \epsilon_0$.

    * `epsilon_0`
        default rms surface accuracy with units of length (`~astropy.units.quanitity.Quantity`)

    * `physical_aperture`
        antenna physical aperture with units of length**2 (`~astropy.units.quanitity.Quantity`)

    """

    def __init__(self):
        self.ap_eff_0 = 1.0
        self.epsilon_0 = 100 * u.micron
        self.physical_aperture = 1.0 * u.m * u.m

    @abstractmethod
    def airmass(self, angle: Union[Angle, Quantity], zd: bool, **kwargs) -> Union[float, np.ndarray]:
        """Computes the airmass at given elevation(s) or zenith distance(s).
        Subclasses should implement an airmass function specific to their application.

        Parameters
        ----------
        angle : `~astropy.coordinates.Angle` or `~astro.units.quantity.Quantity`
            The elevation(s) or zenith distance(s) at which to compute the airmass

        zd: bool
            True if the input value is zenith distance, False if it is elevation.

        **kwargs : Any
            Other possible parameters that affect the airmass, e.g. weather data.

        Returns
        -------
            airmass - float or `~numpy.ndarray`
            The value(s) of the airmass at the given elevation(s)/zenith distance(s)

        """
        pass

    @abstractmethod
    def aperture_efficiency(self, specval: Quantity, **kwargs) -> Union[float, np.ndarray]:
        """
        Calculate the antenna aperture efficiency.

        Parameters
        ----------
        specval : `~astro.units.quantity.Quantity`
            The spectral value -- frequency or wavelength -- at which to compute the efficiency

        **kwargs : Any
            Other possible parameters that affect the aperture efficiency, e.g., elevation angle.

        Returns
        -------
            aperture_efficiency  - float or `~numpy.ndarray`
            The value(s) of the aperture efficiency at the given frequency/wavelength.
            The return value(s) are float(s) between zero and one.

        """
        pass

    def zenith_opacity(self, specval: Quantity, **kwargs) -> Union[float, np.ndarray]:
        """
        Compute the zenith opacity.

        Parameters
        ----------
        specval : `~astro.units.quantity.Quantity`
            The spectral value -- frequency or wavelength -- at which to compute the opacity
        **kwargs : Any
            Other possible parameters that affect the output opacity, e.g., MJD

        Returns
        -------
        float or `~numpy.ndarray`
            The values of the zenith opacity at the given frequency/wavelength.
            The return value(s) are non-negative float(s).
        """
        pass


class GBTGainCorrection(BaseGainCorrection):
    """Gain correction class and functions specific to the Green Bank Telescope.

    Parameters
    -----------
    gain_correction_table : str or `pathlib.Path`
         File to read that contains the parameterized gain correction as a function
         of zenith distance and time (see  `GBT Memo 301 <https://library.nrao.edu/public/memos/gbt/GBT_301.pdf>`_).
         Must be in an `~astropy.table.QTable` readable format.
         Default None will use dysh's internal GBT gain correction table.
    """

    def __init__(self, gain_correction_table: Path = None):

        if gain_correction_table is None:
            gain_correction_table = get_project_configuration() / "gaincorrection.tab"
        self._gct = QTable.read(gain_correction_table, format="ascii.ecsv")
        self._gct.sort("Date")  # just in case it ever is written unsorted.
        self.app_eff_0 = 0.71
        self.epsilon_0 = 230 * u.micron
        self.physical_aperture = 7853.9816 * u.m * u.m

    @property
    def gain_correction_table(self):
        """The table containing the parameterized gain correction as a fucntion of zenith distance and time"""
        return self._gct

    def airmass(self, angle: Union[Angle, Quantity], zd: bool = False, **kwargs) -> Union[float, np.ndarray]:
        """
        Computes the airmass at given elevation(s) or zenith distance(s).  The formula used is

        :math:`A = -0.0234 + 1.014/sin(El+5.18/(El+3.35))`

        for elevation in degrees. This function is specific for the GBT location derived
        from vertical weather data. Source: `(Maddalena 2007) <https://www.gb.nrao.edu/~rmaddale/GBT/Maddalena_HighPrecisionCalibration.pdf>`_

        Parameters
        ----------
        angle :  `~astropy.coordinates.Angle` or `~astro.units.quantity.Quantity`
            The elevation(s) or zenith distance(s) at which to compute the airmass

        zd: bool
            True if the input value is zenith distance, False if it is elevation. Default: False

        Returns
        -------
            airmass - float or `~numpy.ndarray`
            The value(s) of the airmass at the given elevation(s)/zenith distance(s)

        """
        ang_deg = angle.to(u.degree)
        if zd:
            ang_deg = 90.0 * u.degree - ang_deg

        c0 = -0.0234
        c1 = 5.18
        c2 = 3.35
        c3 = 1.014
        d = ang_deg.value + c1 / (ang_deg.value + c2)

        return c0 + c3 / np.sin(np.radians(d))

    def _get_gct_index(self, date: Time) -> int:
        """
        locate the row in GC table that is applicable to the input date.
        Assumes table is sorted (happens in constructor)!

        Parameters
        ----------
        date : `~astropy.time.Time`
            Date of observation

        Returns
        -------
        int
            Index to use from gain correction table

        """
        tablen = len(self._gct)
        index = tablen - 1
        for i in range(tablen):
            if date <= self._gct["Date"][i]:
                index = i
                break
        return index

    def surface_error(self, date: Time) -> Quantity:
        """
        Lookup the applicable surface error in the gain correction table
        for the observation date.

        Parameters
        ----------
        date : `~astropy.time.Time`
            Date of observation

        Returns
        -------
            `~astropy.units.quantity.Quantity`
                Surface error for the given date.

        """
        i = self._get_gct_index(date)
        return self._gct[i]["Surface Error"]

    def gain_correction(
        self,
        angle: Union[Angle, Quantity],
        date: Time,
        zd: bool = True,
    ) -> Union[float, np.ndarray]:
        r"""
        Compute the gain correction scale factor, to be used in the aperture efficiency
        calculation. The factor is a float between zero and 1.
        (See `GBT Memo 301 <https://library.nrao.edu/public/memos/gbt/GBT_301.pdf>`_).
        The factor is determined by:

        :math:`G = A0 + A1*ZD + A2*ZD^2`

        where An are the time-dependent coefficients and ZD is the zenith distance angle in degrees.

        Parameters
        ----------
        angle :  `~astropy.coordinates.Angle` or `~astro.units.quantity.Quantity`
            The elevation(s) or zenith distance(s) at which to compute the gain correction factor

        date  : `~astropy.time.Time`
            The date at which to cmopute the gain correction factor

        zd: bool
            True if the input value is zenith distance, False if it is elevation. Default: False

        Returns
        -------
            gain_correction - float or `~numpy.ndarray`
            The gain correction scale factor(s) at the given elevation(s)/zenith distance(s)

        """
        i = self._get_gct_index(date)
        a0 = self._gct[i]["A0"]
        a1 = self._gct[i]["A1"]
        a2 = self._gct[i]["A2"]
        logger.debug(f"Using {a0=} {a1=} {a2=} for {date}")
        ang_deg = angle.to(u.degree)
        if not zd:
            z = 90.0 - ang_deg.value
        else:
            z = ang_deg.value

        return a0 + a1 * z + a2 * z * z

    def aperture_efficiency(
        self, specval: Quantity, angle: Union[Angle, Quantity], date: Time, zd=False, eps0=None, **kwargs
    ) -> Union[float, np.ndarray]:
        r"""
        Compute the aperture efficiency, as a float between zero and 1. The aperture
        efficiency :math:`\eta_a`, is determined by:

                :math:`\eta_a = \eta_0 G(ZD) \exp(-(4\pi\epsilon_0/\lambda)^2)`

        where :math:`\eta_0` is the long wavelength aperture efficiency, :math:`G(ZD)` is the gain correction factor
        at a zenith distance :math:`ZD, \epsilon_0` is the surface error, and :math:`\lambda` is the wavelength.

        Parameters
        ----------
        specval : `~astro.units.quantity.Quantity`
            The spectral value -- frequency or wavelength -- at which to compute the efficiency

        angle :  `~astropy.coordinates.Angle` or `~astro.units.quantity.Quantity`
            The elevation(s) or zenith distance(s) at which to compute the efficiency

        date  : `~astropy.time.Time`
            The date at which to cmopute the efficieyncy

        zd: bool
            True if the input value is zenith distance, False if it is elevation. Default: False

        eps0 :  `~astro.units.quantity.Quantity` or None
            The value of :math:`\epsilon_0` to use. If given, must have units of length (typically microns).
            If None, the measured value from observatory testing will be used (See :meth:`surface_error`).

        Returns
        -------
            float or `~numpy.ndarray`
            The aperture efficiency at the given inputs

        """
        coeff = self.app_eff_0 * self.gain_correction(angle, date, zd)
        if eps0 is None:
            eps0 = self.surface_error(date)
        _lambda = specval.to(eps0.unit, equivalencies=u.spectral())
        a = (4.0 * np.pi * eps0 / _lambda) ** 2
        eta_a = coeff * np.exp(-a)  # this will be a Quantity with units u.dimensionless
        return eta_a.value

    # @todo (maybe) Formally specval = None is allowed then the script returns data al all frequencies 2-116 GHz
    def get_weather(
        self, specval: Quantity, vartype: str, mjd: Union[Time, float] = None, coeffs=True, **kwargs
    ) -> np.ndarray:
        r"""
        Call the GBO `getForecastValues` script with the given inputs.

        Parameters
        ----------
        specval : `~astro.units.quantity.Quantity`
            The spectral value -- frequency or wavelength -- at which to compute the opacity
        vartype : str, optional
            Which weather variable to fetch. See Notes for a description of valid values.
        mjd : `~astropy.time.Time` or float
            The date at which to compute the opacity. If given as a float, it is interpreted as
            Modified Julian Day.  Default: None, meaning the data will be fetched at the most recent MJD available.
            If the user is not on the GBO network, this argument is ignored and the opacity will only be a function of frequency.
        coeffs: bool
            If True and at GBO, `getForecastValues` will be passed the `-coeffs` argument which returns
            polynomial coefficients to fit opacity as a function of frequency for each MJD.

        Returns
        -------
            `~numpy.ndarray`
            The requested value(s) the given input(s) as a :math:`N_{mjd} \times N_{freq}` array
        """
        if self._forecast is None:
            try:
                self._forecast = GBTWeatherForecast()
            except Exception as e:
                raise Exception(
                    f"Could not create GBTWeatherForecast object because {str(e)} . Are you on the GBO network?"
                )
            return self._forecast.fetch(specval=specval, vartype=vartype, mjd=mjd, coeffs=coeffs)

    # Question: should use_script default to False?
    def zenith_opacity(
        self, specval: Quantity, mjd: Union[Time, float] = None, coeffs=True, use_script=True, **kwargs
    ) -> np.ndarray:
        r"""
        Compute the zenith opacity, optionally interfacing with the GBO `getForecastValues` script.  If multiple `specval` are given, an array is returned otherwise a float is returned.


        Parameters
        ----------
        specval : `~astro.units.quantity.Quantity`
            The spectral value -- frequency or wavelength -- at which to compute the opacity
        mjd : `~astropy.time.Time` or float or list of float
            The date(s) at which to compute the opacity. If given as a float, it is interpreted as
            Modified Julian Day.  Default: None, meaning the data will be fetched at the most recent MJD available.
            If the user is not on the GBO network, this argument is ignored and the opacity will only be a function of frequency.
        coeffs: bool
            If True and at GBO, `getForecastValues` will be passed the `-coeffs` argument which returns
            polynomial coefficients to fit opacity as a function of frequency for each MJD.
        use_script: If at GBO, use the `getForecastValues` script to determine the opacity. This argument is
                    ignore if the user is not on the GBO network.

        Returns
        -------
            `~numpy.ndarray`
            The zenith opacity at the given input(s) as a :math:`N_{mjd} \times N_{freq}` array

        """
        # specval can but value*unit or [value...]*unit.  We want [value...]*unit
        if len(specval.shape) == 0:
            specval = [specval.value] * specval.unit
        if mjd is None:
            mjd = np.zeros(specval.shape)
        mjd_list = to_mjd_list(mjd)
        if use_script:
            return self.get_weather(specval=specval, type="Opacity", mjd=mjd, coeffs=coeffs)
        else:
            frequency = specval.to(u.GHz, equivalencies=u.spectral())
            out = None
            for d in mjd_list:
                for f in frequency:
                    value = [d]  # MJD is irrelevant but need it so output signature is the same
                    value.append(self._default_gbtidl_opacity(f))
                    if out is None:
                        out = np.array(value)
                    else:
                        out = np.vstack([out, value])
        return out

    def atm_temperature(
        self, specval: Quantity, mjd: Union[Time, float] = None, coeffs=True, use_script=True, **kwargs
    ) -> np.ndarray:
        r"""
        Compute the atmospheric temperature `Tatm`, optionally interfacing with the GBO `getForecastValues` script.  If multiple `specval` are given, an array is returned otherwise a float is returned.

        Parameters
        ----------
        specval : `~astro.units.quantity.Quantity`
            The spectral value -- frequency or wavelength -- at which to compute the opacity
        mjd : `~astropy.time.Time` or float
            The date at which to compute the opacity. If given as a float, it is interpreted as
            Modified Julian Day.  Default: None, meaning ignore this parameter. If the user is not on the GBO network,
            this argument is ignored and the opacity will only be a function of frequency.
        coeffs: bool
            If True and at GBO, `getForecastValues` will be passed the `-coeffs` argument which returns
            polynomial coefficients to fit opacity as a function of frequency for each MJD.
        use_script: If at GBO, use the `getForecastValues` script to determine the opacity. This argument is
                    ignore if the user is not on the GBO network.

        Returns
        -------
            `~numpy.ndarray`
            The atmostpheric temperature at the given input(s) as a :math:`N_{mjd} \times N_{freq}` array

        """
        if use_script:
            return self.get_weather(specval=specval, type="Tatm", mjd=mjd, coeffs=coeffs)
        else:
            raise NotImplementedError("Don't yet know how to get Tatm without the getForecasValues script")

    def _default_gbtidl_opacity(self, frequency: Quantity) -> float:
        """Implementation of the GBTIDL method of computing zenith opacity.
        This method is not recommended (even by GBTIDL!). It is implemented here for compatibility only.

        Parameters
        ----------
        frequency : `~astro.units.quantity.Quantity`
            The frequency at which to compute the opacity

        Returns
        -------
            float
            The zenith opacity at the input frequency
        """
        freq = frequency.to(u.GHz).value
        if freq > 52.0:
            return 0.2
        if freq > 18.0 and freq < 26.0:
            tau = 0.008 + np.exp(np.sqrt(freq)) / 8000.0 + np.exp(-((freq - 22.2) ** 2) / 2.0) / 40.0
        else:
            tau = 0.008 + np.exp(np.sqrt(freq)) / 8000.0
        return tau
