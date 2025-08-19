"""
Calibrator class.
"""

import json
import warnings
from collections.abc import Callable
from pathlib import Path
from typing import Literal

import astropy.units as u
import numpy as np
import numpy.typing as npt
from astropy.units import Quantity

from dysh.log import logger
from dysh.util import get_project_data, minimum_string_match


class CalibratorTable:
    r"""
    This class is used to hold calibrator information.
    It uses as input a json file. By default it will look for
    dysh/data/calibrators.json, which is distributed with `dysh`.

    The structure of the input table should be:

    .. code-block::

        {
            "Objects" : {
                "name" : {
                    "LAS"   : number,
                    "cal"   : boolean,
                    "alias" : array
                    },
            "Scale name" : {
                "fluxscale" : string,
                "objects" : {
                    "name" : {
                        "coefs"  : array,
                        "nu_min" : number,
                        "nu_max" : number
                        }
                    },
                "method" : string
                }
            }
        }

    Where "name" is the object name. It should match under the "Objects" and "Scale name" "objects".

    "cal" is a boolean indicating if the object is suitable as a flux calibrator. This is
    mainly used to issue warnings.

    "alias" is a list of alises for the object name.

    "LAS" is the largest angular extent of the object in degrees.

    "Scale name" is the name of the scale, for example "Ott 1994".

    "fluxscale" is a repeat of "Scale name".

    "coefs" are the cofficients that define the spectral energy distribution of the object.
    The order and the meaning of the coefficients depends on the function defined by "method".

    "nu_min" is the minimum frequency, in GHz, where the coeffcients are valid.

    "nu_max" is the maximum frequency, in GHz, where the coeffcients are valid.

    "method" is the method used to compute the spectral energy distribution given the list of
    coefficients provided for each object. The method should be defined in this file, calibrators.py.

    Parameters
    ----------
    calibrator_table_file : `~pathlib.Path` or `str` or `None`
        Path to the json table with calibrator information.
        The contents of the file are defined above.
        If `None`, it will use dysh/data/calibrators.json

    """

    def __init__(self, calibrator_table_file: Path | str | None = None):
        if calibrator_table_file is None:
            self.calibrator_table_file = get_project_data() / "calibrators.json"
        else:
            self.calibrator_table_file = Path(calibrator_table_file)
        self.load()
        self._make_aliases()

    def load(self):
        """
        Load the json table with calibrator information.
        """
        with open(self.calibrator_table_file) as f:
            self.data = json.load(f)

    def valid_names(self, scale: None | str = None):
        """
        Recognized calibrator names and their alises
        defined in the json table.

        Parameters
        ----------
        scale : str or None
            If set to None will return all the objects names found in the table.
            If set to a known flux scale it will return the names of the objects
            defined for that flux scale.
        """
        names = []
        if scale is None:
            for k, v in self.data["Objects"].items():
                names.append(k)
                names += v["alias"]
        else:
            for k in self.data[scale]["objects"].keys():
                names.append(k)
                names += self.data["Objects"][k]["alias"]
        # This sorts alphanumerically, which is easier for humans
        # to read.
        return sorted(set(names))

    def valid_scales(self):
        """
        Valid flux scales defined in the json table.
        """
        keys = list(self.data.keys())
        keys.remove("Objects")
        return keys

    def _make_aliases(self):
        self.alias = {}
        for k, v in self.data["Objects"].items():
            self.alias[k] = k
            for vv in v["alias"]:
                self.alias[vv] = k


class Calibrator:
    """
    Holds calibrator parameters for a specific calibrator.
    Using these parameters it is possible to get the flux density
    for the calibrator using the `Calibrator.compute_sed` method.

    Parameters
    ----------
    name : str
        Calibrator name.
    scale : str
        Scale in which the calibrator flux density is defined.
    coefs : list
        List of coefficients that define the flux density as a function of frequency.
        The exact definition of the cofficients is determined by the `method` used.
    method : Callable
        Function used to compute the flux density of the calibrator as a function
        of frequency given the input `coefs`.
    calibrator_table : `~dysh.util.calibrator.CalibratorTable`
        Table with calibrator details. This is only used when creating a Calibrator using the
        `from_name` class method. Otherside it is kept for reference.
    las : float or None
        Largest angular size for the calibrator in degrees.
    cal : bool
        Is the calibrator suitable for flux density calibration?
        This is used to issue warnings.
    nu_min : `~astropy.units.Quantity` or None
        Minimum frequency over which the coefficients are valid.
    nu_max : `~astropy.units.Quantity` or None
        Maximum frequency over which the coefficients are valid.

    Examples
    --------

    Create a `~dysh.util.calibrator.Calibrator` for 3C123 using the Perley & Butler 2017 flux scale.

    >>> from dysh.util import calibrator
    >>> c = calibrator.Calibrator.from_name("3C123")

    Compute the flux density at 1 GHz

    >>> import astropy.units as u
    >>> c.compute_sed(1*u.GHz).value
    63.34320007697665

    """

    @u.quantity_input(nu_min=u.GHz, equivalencies=u.spectral())
    @u.quantity_input(nu_max=u.GHz, equivalencies=u.spectral())
    def __init__(
        self,
        name: str,
        scale: str,
        coefs: list,
        method: Callable,
        calibrator_table: CalibratorTable | None = None,
        las: float | None = None,
        cal: bool = True,
        nu_min: u.Quantity | None = None,
        nu_max: u.Quantity | None = None,
    ):
        if calibrator_table is None:
            self.calibrator_table = CalibratorTable()
        else:
            self.calibrator_table = calibrator_table

        self.name = name
        self.scale = scale
        self.coefs = coefs
        self.method = method
        self.las = las
        self.cal = cal
        self.nu_min = nu_min
        self.nu_max = nu_max

    @classmethod
    def from_name(
        cls,
        name: str,
        scale: Literal["Perley-Butler 2017", "Ott 1994"] | None = None,
        calibrator_table: CalibratorTable | None = None,
    ):
        """
        Create a `~dysh.util.calibrator.Calibrator` given a calibrator `name`.

        Parameters
        ----------
        name : str
            Calibrator name, case insensitive.
        scale : None or str or "Perley-Butler 2017" or "Ott 1994"
            The name of the flux scale to use, case insensitive. Minimum string match allowed.
            If set to None, the default, it will use the Perley & Butler 2017 scale
            if available, otherwise Ott et al. 1994.
            Only the Perley & Butler 2017 and Ott et al. 1994 scales
            are shipped with `dysh`.
            If a custom `calibrator_table` is provided this can be the name of a scale
            defined in that table.

        Raises
        ------
        ValueError
            If `scale` is not one of the keys in the `calibrator_table`,
            or if `name` is not one of the keys in `calibrator_table.data[scale]`.
        """

        if calibrator_table is None:
            calibrator_table = CalibratorTable()

        # Make them lower case for user friendliness.
        _name = name.lower()
        _valid_names = list(map(str.lower, calibrator_table.valid_names()))
        _valid_scales = list(map(str.lower, calibrator_table.valid_scales()))

        if _name not in _valid_names:
            raise ValueError(
                f"Unrecognized calibrator name {name}. Valid names: {', '.join(calibrator_table.valid_names())}"
            )

        if scale is not None and minimum_string_match(scale.lower(), _valid_scales) is None:
            raise ValueError(f"Unrecognized scale {scale}. Valid scales: {', '.join(calibrator_table.valid_scales())}")
        elif scale is not None:
            scale = minimum_string_match(scale.lower(), _valid_scales).title()
            if _name not in list(map(str.lower, calibrator_table.valid_names(scale.title()))):
                raise ValueError(f"Calibrator {name} is not defined in scale {scale}. Try using 'scale=None'")
        elif scale is None:
            for vs in calibrator_table.valid_scales():
                if _name in list(map(str.lower, calibrator_table.valid_names(vs))):
                    scale = vs
                    logger.info(f"Using {vs} for {name}.")
                    break

        alias = calibrator_table.alias[name.title()]
        coefs = calibrator_table.data[scale]["objects"][alias]["coefs"]
        method = globals()[calibrator_table.data[scale]["method"]]
        las = calibrator_table.data["Objects"][alias]["LAS"]
        cal = calibrator_table.data["Objects"][alias]["cal"]
        nu_min = calibrator_table.data[scale]["objects"][alias]["nu_min"] * u.GHz
        nu_max = calibrator_table.data[scale]["objects"][alias]["nu_max"] * u.GHz

        return cls(name, scale, coefs, method, calibrator_table, las=las, cal=cal, nu_min=nu_min, nu_max=nu_max)

    def compute_sed(self, nu: Quantity, hpbw: float | None = None):
        """
        Evaluate the calibrator flux density at `nu`.

        Parameters
        ----------
        nu : `~astropy.units.Quantity`
            Frequency values.
        hpbw : float or None
            Half Power Beam Width in degrees.

        Warns
        -----
        UserWarning
            If the `Calibrator` is not suitable as a flux density calibrator, if `nu` is above or below the range where `coefs` are defined,
            of if `hpbw` is smaller than the largest angular size of the calibrator.
        """

        if not self.cal:
            warnings.warn(
                "Source {self.name} deemed unsuitable as a flux calibrator. Proceed with caution.", stacklevel=2
            )

        with u.set_enabled_equivalencies(u.spectral()):
            if self.nu_min is not None and np.min(nu) < self.nu_min:
                warnings.warn(
                    f"Frequency ({np.min(nu)}) is below the minimum frequency where the flux scale is defined ({self.nu_min}). Proceed with caution.",
                    stacklevel=2,
                )

        with u.set_enabled_equivalencies(u.spectral()):
            if self.nu_max is not None and np.max(nu) > self.nu_max:
                warnings.warn(
                    f"Frequency ({np.max(nu)}) is above the maximum frequency where the flux scale is defined ({self.nu_max}). Proceed with caution.",
                    stacklevel=2,
                )

        if hpbw is not None and self.las is not None and hpbw < self.las:
            warnings.warn(
                f"Largest angular size ({self.las}) is larger than the HPBW ({hpbw}). Proceed with caution.",
                stacklevel=2,
            )

        # Calibrator flux density in Jy.
        snu = self.method(nu, self.coefs)

        return snu


def poly_pb(nu: Quantity, coefs: npt.ArrayLike):
    r"""
    Equation (1) in Perley & Butler 2017 [1]_.

    Parameters
    ----------
    nu : `~astropy.units.Quantity`
        Frequency values.
    coefs : array
        Coefficients that define the radio spectral energy distribution.

    Returns
    -------
    snu : `~astropy.units.Quantity`
        Radio flux density evaluated at `nu` in Jy.

    Notes
    -----
    The flux density, :math:`S`, is computed from

    .. math::

        \log(S)=a_{0}+a_{1}\log(\nu)+a_{2}\log^{2}(\nu)+\cdots

    with :math:`a_{0}` the last element of `coefs`, :math:`a_{1}` the
    second to last, and so on. :math:`\nu` is the frequency in GHz.

    References
    ----------

    .. [1] `Perley & Butler (2017) <https://ui.adsabs.harvard.edu/abs/2017ApJS..230....7P/abstract>`_
    """

    # Make sure the frequency is in GHz.
    nu = nu.to(u.GHz).value

    # Calibrator flux density in Jy.
    snu = np.power(10.0, np.polyval(coefs, np.log10(nu))) * u.Jy

    return snu


def poly_ott(nu: Quantity, coefs: list):
    r"""
    From Table 5 in Ott et al. 1994 [1]_.

    Parameters
    ----------
    nu : `~astropy.units.Quantity`
        Frequency values.
    coefs : array
        Coefficients that define the radio spectral energy distribution.
        For the Ott et al. 1994 model this must have three values.

    Returns
    -------
    snu : `~astropy.units.Quantity`
        Radio flux density evaluated at `nu` in Jy.

    Notes
    -----
    The flux density, :math:`S`, is computed from

    .. math::

        \log(S)=a+b\log(\nu)+c\log^{2}(\nu)

    with :math:`a` the first element of `coefs` and :math:`\nu` the frequency in MHz.

    References
    ----------

    .. [1] `Ott et al. (1994) <https://ui.adsabs.harvard.edu/abs/1994A%26A...284..331O/abstract>`_
    """

    nu = nu.to(u.MHz).value
    lognu = np.log10(nu)

    snu = np.power(10.0, (coefs[0] + coefs[1] * lognu + coefs[2] * np.power(lognu, 2.0))) * u.Jy

    return snu
