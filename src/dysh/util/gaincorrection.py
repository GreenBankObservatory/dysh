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
        angle : ~astropy.coordinates.Angle or ~astro.units.quantity.Quantity
            The elevation(s) or zenith distance(s) at which to compute the airmass

        zd: bool
            True if the input value is zenith distance, False if it is elevation.

        **kwargs : Any
            Other possible parameters that affect the airmass, e.g. weather data.

        Returns
        -------
            airmass - float or ~np.ndarray
            The value(s) of the airmass at the given elevation(s)/zenith distance(s)
        """
        pass

    @abstractmethod
    def aperture_efficiency(self, specval: Quantity, **kwargs) -> Union[float, np.ndarray]:
        """
        Calculate the antenna aperture efficiency.

        Parameters
        ----------
        frequency : Quantity
            The frequency at which to calculate the efficiency.

        **kwargs : Any
            Other possible parameters that affect the aperture efficiency, e.g., elevation angle.

        Returns
        -------
            aperture_efficiency  - float or ~np.ndarray
            The value(s) of the aperture efficiency at the given frequency.
            The return value(s) are float(s) between zero and one.

        """
        pass


class GBTGainCorrection(BaseGainCorrection):
    """Gain correction class and functions specific to the Green Bank Telescope"""

    def __init__(self, gain_correction_table: Path = None):
        """
        Parameters
        -----------
        gain_correction_table : str or `pathlib.Path`
             File to read that contains the parameterized gain correction as a function
             of zenith distance and time (see GBT Memo 301: https://library.nrao.edu/public/memos/gbt/GBT_301.pdf).
             Must be in an `~astropy.table.QTable` readable format.
             Default None will usedysh's internal GBT gain correction table.
        """
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

        A = -0.0234 + 1.014/sin(El+5.18/(El+3.35))

        for elevation in degrees. This function is specific for the GBT location derived
        from vertical weather data. Source: (Maddalena 2007)
        https://www.gb.nrao.edu/~rmaddale/GBT/Maddalena_HighPrecisionCalibration.pdf

        Parameters
        ----------
        angle :  ~astropy.coordinates.Angle or ~astro.units.quantity.Quantity
            The elevation(s) or zenith distance(s) at which to compute the airmass

        zd: bool
            True if the input value is zenith distance, False if it is elevation. Default: False

        Returns
        -------
            airmass - float or ~np.ndarray
            The value(s) of the airmass at the given elevation(s)/zenith distance(s)

        """
        ang_deg = angle.to(u.degree)
        if zd:
            ang_deg.value = 90.0 - ang_deg.value

        c1 = 5.18 * u.degree
        c2 = 3.35 * u.degree
        return -0.0234 + 1.014 / np.sin(ang_deg + c1 / (ang_deg + c2))

    def _get_gct_index(self, date: Time) -> int:
        """
        locate the row in GC table that is applicable to the input date.
        Assumes table is sorted (happens in constructor)!

        Parameters
        ----------
        date : ~astropy.time.Time
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
        date : ~astropy.time.Time
            Date of observation

        Returns
        -------
        ~astropy.units.quantity.Quantity
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
        """
        Compute the gain correction scale factor, to be used in the aperture efficiency
        calculation. The factor is a float between zero and 1.  See GBT Memo 301. The factor is
        determined by:

        G = A0 + A1*ZD + A2*ZD^2

        where An are the time-dependent coefficients and ZD is the zenith distance angle in degrees.

        Parameters
        ----------
        angle :  ~astropy.coordinates.Angle or ~astro.units.quantity.Quantity
            The elevation(s) or zenith distance(s) at which to compute the gain correction factor

        date  : ~astropy.time.Time
            The date at which to cmopute the gain correction factor

        zd: bool
            True if the input value is zenith distance, False if it is elevation. Default: False

        Returns
        -------
            gain_correction - float or np.ndarray
            The gain correction scale factor(s) at the given elevation(s)/zenith distance(s)
        """
        i = self._get_gct_index(date)
        a0 = self._gct[i]["A0"]
        a1 = self._gct[i]["A1"]
        a2 = self._gct[i]["A2"]
        print(f"Using {a0=} {a1=} {a2=} for {date}")
        ang_deg = angle.to(u.degree)
        if not zd:
            z = 90.0 - ang_deg.value
        else:
            z = ang_deg.value

        return a0 + a1 * z + a2 * z * z

    def aperture_efficiency(
        self, specval: Quantity, angle: Union[Angle, Quantity], date: Time, zd=False, **kwargs
    ) -> Union[float, np.ndarray]:
        r"""
        Compute the aperture efficiency, as a float between zero and 1. The aperture
        efficiency $\eta_a$, is determined by:

                $$\eta_a = \eta_0 G(ZD) \exp(-(4\pi\epsilon_0/\lambda)^2)$$

        where $\eta_0$ is the long wavelength aperture efficiency, $G(ZD)$ is the gain correction factor
        at a zenith distance $ZD$, $\epsilon_0$ is the surface error, and $\lambda$ is the wavelength.

        Parameters
        ----------
        specval : ~astro.units.quantity.Quantity
            The spectral value -- frequency or wavelength -- at which to compute the efficiency

        angle :  ~astropy.coordinates.Angle or ~astro.units.quantity.Quantity
            The elevation(s) or zenith distance(s) at which to compute the efficiency

        date  : ~astropy.time.Time
            The date at which to cmopute the efficieyncy

        zd: bool
            True if the input value is zenith distance, False if it is elevation. Default: False

        Returns
        -------
            float or np.ndarray
            The aperture efficiency at the given inputs

        """
        coeff = self.app_eff_0 * self.gain_correction(angle, date, zd)
        se = self.surface_error(date)
        _lambda = specval.to(se.unit, equivalencies=u.spectral())
        a = (4.0 * np.pi * se / _lambda) ** 2
        eta_a = coeff * np.exp(-a)  # this will be a Quantity with units u.dimensionless
        return eta_a.value
