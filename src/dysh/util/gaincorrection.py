from abc import ABC, abstractmethod
from pathlib import Path

import astropy.constants as ac
import astropy.units as u
import numpy as np
from astropy.coordinates import Angle
from astropy.table import QTable
from astropy.time import Time
from astropy.units.quantity import Quantity

from ..log import logger
from ..util import get_project_data, to_quantity_list
from .core import to_mjd_list
from .weatherforecast import GBTWeatherForecast

__all__ = ["BaseGainCorrection", "GBTGainCorrection"]


class BaseGainCorrection(ABC):
    r"""This class is the base class for gain corrections. It is intended to be subclassed for
    specific antennas. Subclasses will be used to calculate various gain corrections to go
    from antenna temperature to other scales like brightness temperature of flux density.
    Subclasses can implement the following attributes:

    * `ap_eff_0`
        Long wavelength aperture efficiency (number between 0 and 1), :math:`\eta_{a}`, i.e., in the absence of surface errors,
        :math:`\lambda >> \epsilon_0`.

    * `epsilon_0`
        Default rms surface accuracy with units of length (`~astropy.units.quanitity.Quantity`)

    * `physical_aperture`
        Antenna physical aperture with units of length**2 (`~astropy.units.quanitity.Quantity`)

    * `loss_eff_0`
        The telescope efficiency combining radiation efficiency :math:`\eta_r` and rearward scattering and spillover efficiency, :math:`\eta_{rss}`.
        :math:`\eta_{loss} = \eta_r\eta_{rss}`.  This is the term :math:`\eta_l` as defined by Kutner & Ulich (1981) equation 12.
        See https://articles.adsabs.harvard.edu/pdf/1981ApJ...250..341K

    """

    def __init__(self):
        self.ap_eff_0 = 1.0
        self.epsilon_0 = 100 * u.micron
        self.physical_aperture = 1.0 * u.m * u.m
        self.loss_eff_0 = 1.0

    @property
    def jyperk(self):
        r"""The default Gain off the telescope in Jy/K, :math:`G = 2 k_B/A_p`, where `k_B` is Boltzmann's constant
        and `A_p` is the area of the physical aperture of the telescope.
        """
        return (2.0 * ac.k_B / self.physical_aperture.to("m^2")).to("Jy/K")

    @abstractmethod
    def airmass(self, angle: Angle | Quantity, zd: bool, **kwargs) -> float | np.ndarray:
        """Computes the airmass at given elevation(s) or zenith distance(s).
        Subclasses should implement an airmass function specific to their application.

        Parameters
        ----------
        angle : `~astropy.coordinates.Angle` or `~astropy.units.quantity.Quantity`
            The elevation(s) or zenith distance(s) at which to compute the airmass

        zd: bool
            True if the input value is zenith distance, False if it is elevation.

        **kwargs : Any
            Other possible parameters that affect the airmass, e.g. weather data.

        Returns
        -------
        airmass : float or `~numpy.ndarray`
            The value(s) of the airmass at the given elevation(s)/zenith distance(s)

        """
        pass

    @abstractmethod
    def aperture_efficiency(self, specval: Quantity, **kwargs) -> float | np.ndarray:
        """
        Calculate the antenna aperture efficiency.

        Parameters
        ----------
        specval : `~astropy.units.quantity.Quantity`
            The spectral value -- frequency or wavelength -- at which to compute the efficiency

        **kwargs : Any
            Other possible parameters that affect the aperture efficiency, e.g., elevation angle.

        Returns
        -------
        aperture_efficiency : float or `~numpy.ndarray`
            The value(s) of the aperture efficiency at the given frequency/wavelength.
            The return value(s) are float(s) between zero and one.

        """
        pass

    def zenith_opacity(self, specval: Quantity, **kwargs) -> float | np.ndarray:  # noqa: B027
        """
        Compute the zenith opacity.

        Parameters
        ----------
        specval : `~astropy.units.quantity.Quantity`
            The spectral value -- frequency or wavelength -- at which to compute the opacity
        **kwargs : Any
            Other possible parameters that affect the output opacity, e.g., MJD

        Returns
        -------
        tau : float or `~numpy.ndarray`
            The values of the zenith opacity at the given frequency/wavelength.
            The return value(s) are non-negative float(s).
        """
        pass


class GBTGainCorrection(BaseGainCorrection):
    """Gain correction class and functions specific to the Green Bank Telescope.

    Parameters
    ----------
    gain_correction_table : str or `pathlib.Path`
         File to read that contains the parameterized gain correction as a function
         of zenith distance and time (see  `GBT Memo 301 <https://library.nrao.edu/public/memos/gbt/GBT_301.pdf>`_).
         Must be in an `~astropy.table.QTable` readable format.
         Default None will use dysh's internal GBT gain correction table.

    Attributes
    ----------
    valid_scales : tuple
        Strings representing valid options for scaling spectral data, specifically

        - 'ta'  : Antenna Temperature in K
        - 'ta*' : Antenna temperature corrected to above the atmosphere in K
        - 'flux'  : flux density in Jansky

    """

    valid_scales: tuple[str, str, str] = ("ta", "ta*", "flux")

    def __init__(self, gain_correction_table: Path = None):  # noqa: RUF013
        if gain_correction_table is None:
            gain_correction_table = get_project_data() / "gaincorrection.tab"
        self._gct = QTable.read(gain_correction_table, format="ascii.ecsv")
        self._gct.sort("Date")  # just in case it ever is written unsorted.
        self.ap_eff_0 = 0.71
        self.epsilon_0 = 230 * u.micron
        self.physical_aperture = 7853.9816 * u.m * u.m
        self.loss_eff_0 = 0.99
        self._forecast = None

    @classmethod
    def is_valid_scale(cls, scale):
        """
        Check that a string represents a valid option for scaling spectral data.
        See: `valid_scales`.

        Parameters
        ----------
        scale : str
            temperature scale descriptive string.

        Returns
        -------
        bool
            True if `scale` is a valid scaling option, False otherwise

        """
        return scale.lower() in cls.valid_scales

    @property
    def gain_correction_table(self):
        """The table containing the parameterized gain correction as a fucntion of zenith distance and time"""
        return self._gct

    def airmass(self, angle: Angle | Quantity, zd: bool = False, **kwargs) -> float | np.ndarray:
        """
        Computes the airmass at given elevation(s) or zenith distance(s).  The formula used is

        :math:`A = -0.0234 + 1.014/sin(El+5.18/(El+3.35))`

        for elevation in degrees. This function is specific for the GBT location derived
        from vertical weather data. Source: `(Maddalena 2007) <https://www.gb.nrao.edu/~rmaddale/GBT/Maddalena_HighPrecisionCalibration.pdf>`_

        Parameters
        ----------
        angle :  `~astropy.coordinates.Angle` or `~astropy.units.quantity.Quantity`
            The elevation(s) or zenith distance(s) at which to compute the airmass

        zd: bool
            True if the input value is zenith distance, False if it is elevation. Default: False

        Returns
        -------
        airmass : float or `~numpy.ndarray`
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
        Locate the row in GC table that is applicable to the input date.
        Assumes table is sorted (happens in constructor)!

        Parameters
        ----------
        date : `~astropy.time.Time`
            Date of observation

        Returns
        -------
        index : int
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
        angle: Angle | Quantity,
        date: Time,
        zd: bool = True,
    ) -> float | np.ndarray:
        r"""
        Compute the gain correction scale factor, to be used in the aperture efficiency
        calculation. The factor is a float between zero and 1.
        (See `GBT Memo 301 <https://library.nrao.edu/public/memos/gbt/GBT_301.pdf>`_).
        The factor is determined by:

        :math:`G = A0 + A1*ZD + A2*ZD^2`

        where An are the time-dependent coefficients and ZD is the zenith distance angle in degrees.

        Parameters
        ----------
        angle :  `~astropy.coordinates.Angle` or `~astropy.units.quantity.Quantity`
            The elevation(s) or zenith distance(s) at which to compute the gain correction factor

        date  : `~astropy.time.Time`
            The date at which to compute the gain correction factor

        zd: bool
            True if the input value is zenith distance, False if it is elevation. Default: False

        Returns
        -------
        gain_correction : float or `~numpy.ndarray`
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
        self,
        specval: Quantity,
        angle: Angle | Quantity,
        date: Time,
        zd: bool = False,
        surface_error: None | Quantity = None,
        **kwargs,
    ) -> float | np.ndarray:
        r"""
        Compute the aperture efficiency, as a float between zero and 1. The aperture
        efficiency :math:`\eta_a`, is determined by:

        .. math::

            \eta_a = \eta_0 G(ZD) \exp(-(4\pi\epsilon/\lambda)^2)

        where :math:`\eta_0` is the long wavelength aperture efficiency, :math:`G(ZD)` is the gain correction factor
        at a zenith distance :math:`ZD` (`zd`), :math:`\epsilon` (`surface_error`) is the surface error, and :math:`\lambda` is the wavelength.

        **Rules for input of multiple dates, spectral values, and angles**

        For a single date:

            - If one spectral value and multiple angles, then the aperture efficiency at each angle is returned.
            - If multiple spectral values and one angle, then the aperture efficiency at each spectral value is returned
            - If multiple spectral values and multiple angles, then it is assumed they are to be paired
              and the aperture efficiency at each pair will be returned.  Therefore the lengths must be equal.

        For mutiple dates:

            - For one spectral value and one angle, the aperture efficiency at each date is returned.
            - For multiple spectral values and multiple angles, it is assumed they go together, so the lengths must match
              the number of dates. The aperture efficiency for each (spectral value, angle, date) tuple will be returned.

        Parameters
        ----------
        specval : `~astropy.units.quantity.Quantity`
            The spectral value(s) -- frequency or wavelength -- at which to compute the efficiency

        angle :  `~astropy.coordinates.Angle` or `~astropy.units.quantity.Quantity`
            The elevation(s) or zenith distance(s) at which to compute the efficiency

        date  : `~astropy.time.Time`
            The date(s) at which to compute the efficiency.

        zd : bool
            True if the input value is zenith distance, False if it is elevation. Default: False

        surface_error : `~astropy.units.quantity.Quantity` or None
            The value of :math:`\epsilon` to use, the surface rms error. If given, must have units of length (typically microns).
            If None, the measured value from observatory testing will be used (See :meth:`surface_error`).

        Returns
        -------
        eta_a : float or `~numpy.ndarray`
            The aperture efficiency at the given inputs

        """

        # The code easier to read if we make them all Quantities with dimensioned values.
        sp = to_quantity_list(specval)
        ang = to_quantity_list(angle)
        if date.isscalar:
            if len(sp) == 1 or len(ang) == 1 or (len(sp) == len(ang)):
                coeff = self.ap_eff_0 * self.gain_correction(ang, date, zd)
                if surface_error is None:
                    surface_error = self.surface_error(date)
            else:
                raise ValueError(f"Number of specvals {len(sp)} and angles {len(ang)} must be equal")
        else:
            if all([len(sp) == 1, len(ang) == 1]):
                coeff = []
                for i in range(len(date)):
                    coeff.append(self.ap_eff_0 * self.gain_correction(angle, date[i], zd))
                coeff = np.squeeze(np.array(coeff))
            elif len(sp) == len(date) and len(ang) == len(date):
                coeff = []
                for i in range(len(date)):
                    coeff.append(self.ap_eff_0 * self.gain_correction(ang[i], date[i], zd))
                coeff = np.squeeze(np.array(coeff))
            else:
                raise ValueError(
                    f"Number of specvals {len(sp)} and angles {len(ang)} must be equal to number of dates {len(date)}"
                )

            if surface_error is None:
                surface_error = self._surface_error_array(date)
        _lambda = sp.to(surface_error.unit, equivalencies=u.spectral())

        a = (4.0 * np.pi * surface_error / _lambda) ** 2
        eta_a = coeff * np.exp(-a)  # this will be a Quantity with units u.dimensionless
        return eta_a.value

    def _surface_error_array(self, date):
        return to_quantity_list([self.surface_error(x) for x in date])

    def scale_ta_to(
        self,
        tscale: str,
        specval: Quantity,
        angle: Angle | Quantity,
        date: Time,
        zenith_opacity,
        zd=False,
        ap_eff=None,
        surface_error=None,
    ) -> float | np.ndarray:
        r"""
        Scale the antenna temperature to a different brightness temperature unit.

        tscale : str
            The brightness scale unit for the output scan, must be one of (case-insensitive)
                - 'Ta'  : Antenna Temperature in K
                - 'Ta*' : Antenna temperature corrected to above the atmosphere in K
                - 'Flux'  : flux density in Jansky

            If 'Ta*' or 'Flux' the zenith opacity must also be given.

        specval : `~astropy.units.quantity.Quantity`
            The spectral value(s) -- frequency or wavelength -- at which to compute the efficiency

        angle :  `~astropy.coordinates.Angle` or `~astropy.units.quantity.Quantity`
            The elevation(s) or zenith distance(s) at which to compute the efficiency

        date  : `~astropy.time.Time`
            The date(s) at which to compute the efficiency.

        zenith_opacity: float
            The zenith opacity to use in calculating the scale factors for the integrations.

        zd : bool
            True if the input value is zenith distance, False if it is elevation. Default: False

        ap_eff : float or None
            Aperture efficiency to be used when scaling data to brightness temperature of flux. The provided aperture
            efficiency must be a number between 0 and 1.  If None, `dysh` will calculate it :meth:`aperture_efficiency`. Only one of `ap_eff` or `surface_errors`
            can be provided.

        surface_error : `~astropy.units.quantity.Quantity` or None
            The value of :math:`\epsilon_0` to use, the surface rms error. If given, must have units of length (typically microns).
            If None, the measured value from observatory testing will be used (See :meth:`surface_error`).
        """
        if not GBTGainCorrection.is_valid_scale(tscale):
            raise ValueError(
                f"Unrecognized temperature scale {tscale}. Valid options are {GBTGainCorrection.valid_scales} (case-insensitive)."
            )
        if ap_eff is not None and surface_error is not None:
            raise ValueError("Only one of ap_eff or surface_error should be specified")
        s = tscale.lower()
        if s == "ta":
            return 1.0
        am = self.airmass(angle, zd)
        if ap_eff is None:
            eta_a = self.aperture_efficiency(specval, angle, date, zd, surface_error)
        else:
            eta_a = ap_eff
        # Calculate Ta* because in both cases we need it
        # Ta* = T_a exp(tau*A)/(eta_a * eta_loss )
        # - the airmass as a function of elevation, A
        # - the aperture efficiency as a function of frequency and date, eta_a
        # - the loss efficiency, eta_loss
        factor = np.exp(zenith_opacity * am) / (eta_a * self.loss_eff_0)
        if s == "ta*":
            return factor
        # Snu = 2kT_a exp(tau*A)/(eta_a * eta_loss *  A_p )
        #     = 2kTa*/A_p
        # where
        # - k is Boltzmann's constant
        # - the physical aperture, A_p
        jyperk = factor * self.jyperk
        return jyperk.value

    def get_weather(
        self, specval: Quantity, vartype: str, mjd: Time | float = None, coeffs: bool = True, **kwargs
    ) -> np.ndarray:
        r"""
        Call the GBO `getForecastValues` script with the given inputs.
        For frequencies below 2 GHz, the value at 2 GHz will be returned since the `getForecastValues` does not
        cover < 2GHz.  Returned values will be sorted by frequency, low to high.

        See `GBTdocs <https://gbtdocs.readthedocs.io/en/latest/how-tos/data_reduction/calculate_opacity.html#calculate-sky-opacity-retrieving-more-granular-values-of-zenith-opacity>`_
        for more details.


        Parameters
        ----------
        specval : `~astropy.units.quantity.Quantity`, optional
            The spectral value -- frequency or wavelength -- at which to compute `vartype`.
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

        Returns
        -------
        weather_data : `~numpy.ndarray`
            The requested value(s) at the given input(s) as a :math:`N_{mjd} \times N_{freq}` array
        """
        if self._forecast is None:
            try:
                self._forecast = GBTWeatherForecast()
            except Exception as e:
                self._forecast = None
                raise Exception(  # noqa: B904
                    f"Could not create GBTWeatherForecast object because {e!s} . Are you on the GBO network?"
                )
        if mjd is None:
            mjd = Time.now()
        return self._forecast.fetch(specval=specval, vartype=vartype, mjd=mjd, coeffs=coeffs)

    # Question: should use_script default to False?
    def zenith_opacity(
        self, specval: Quantity, mjd: Time | float = None, coeffs: bool = True, use_script: bool = True, **kwargs
    ) -> np.ndarray:
        r"""
        Compute the zenith opacity, optionally interfacing with the GBO `getForecastValues` script.  If multiple `specval` are given, an array is returned otherwise a float is returned.

        For frequencies below 2 GHz, the value at 2 GHz will be returned since the `getForecastValues` does not
        cover < 2GHz.  Returned values will be sorted by frequency, low to high.

        Parameters
        ----------
        specval : `~astropy.units.quantity.Quantity`
            The spectral value -- frequency or wavelength -- at which to compute the zenith opacity.
        mjd : `~astropy.time.Time` or float or list of float
            The date(s) at which to compute the opacity. If given as a float, it is interpreted as
            Modified Julian Day.  Default: None, meaning the data will be fetched at the most recent MJD available.
            If the user is not on the GBO network, this argument is ignored and the opacity will only be a function of frequency.
        coeffs : bool
            If True and at GBO, `getForecastValues` will be passed the `-coeffs` argument which returns
            polynomial coefficients to fit opacity as a function of frequency for each MJD.
        use_script : bool
            If at GBO, use the `getForecastValues` script to determine the opacity. This argument is
            ignored if the user is not on the GBO network.

        Returns
        -------
        tau : `~numpy.ndarray`
            The zenith opacity at the given input(s) as a :math:`\ N_{mjd} \times N_{freq}\ ` array.

        """
        # specval can but value*unit or [value...]*unit.  We want [value...]*unit
        if len(specval.shape) == 0:
            specval = [specval.value] * specval.unit
        if use_script:
            return self.get_weather(specval=specval, vartype="Opacity", mjd=mjd, coeffs=coeffs)
        else:
            frequency = specval.to(u.GHz, equivalencies=u.spectral())
            out = None
            if mjd is None:
                mjd = Time.now()
            mjd_list = to_mjd_list(mjd)
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
        self, specval: Quantity, mjd: Time | float = None, coeffs: bool = True, use_script: bool = True, **kwargs
    ) -> np.ndarray:
        r"""
        Compute the atmospheric temperature `Tatm`, optionally interfacing with the GBO `getForecastValues` script.  If multiple `specval` are given, an array is returned otherwise a float is returned.

        For frequencies below 2 GHz, the value at 2 GHz will be returned since the `getForecastValues` does not
        cover < 2GHz.  Returned values will be sorted by frequency, low to high.

        Parameters
        ----------
        specval : `~astropy.units.quantity.Quantity`
            The spectral value -- frequency or wavelength -- at which to compute `Tatm`.
        mjd : `~astropy.time.Time` or float
            The date at which to compute `Tatm`. If given as a float, it is interpreted as
            Modified Julian Day.  Default: None, meaning ignore this parameter. If the user is not on the GBO network,
            this argument is ignored and the opacity will only be a function of frequency.
        coeffs : bool
            If True and at GBO, `getForecastValues` will be passed the `-coeffs` argument which returns
            polynomial coefficients to fit `Tatm` as a function of frequency for each MJD.
        use_script : bool
            If at GBO, use the `getForecastValues` script to determine `Tatm`. This argument is
            ignored if the user is not on the GBO network.

        Returns
        -------
        tatm : `~numpy.ndarray`
            The atmospheric temperature at the given input(s) as a :math:`\ N_{mjd} \times N_{freq}\ ` array

        """
        if use_script:
            return self.get_weather(specval=specval, type="Tatm", mjd=mjd, coeffs=coeffs)
        else:
            raise NotImplementedError("Don't yet know how to get Tatm without the getForecasValues script")

    def _default_gbtidl_opacity(self, frequency: Quantity) -> float:
        r"""Implementation of the GBTIDL method of computing zenith opacity.
        This method is not recommended (even by GBTIDL!). It is implemented here for compatibility only.

        Parameters
        ----------
        frequency : `~astropy.units.quantity.Quantity`
            The frequency at which to compute the zenith opacity.

        Returns
        -------
        tau : float
            The zenith opacity at the input frequency.
        """
        logger.warning("Using default time-independent zenith opacity. This is not recommended.")
        freq = frequency.to(u.GHz).value
        if freq > 52.0:
            return 0.2
        if freq > 18.0 and freq < 26.0:
            tau = 0.008 + np.exp(np.sqrt(freq)) / 8000.0 + np.exp(-((freq - 22.2) ** 2) / 2.0) / 40.0
        else:
            tau = 0.008 + np.exp(np.sqrt(freq)) / 8000.0
        return tau
