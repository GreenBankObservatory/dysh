"""
Calibrator class.
"""

import json
from collections.abc import Callable
from pathlib import Path
from typing import Literal

import astropy.units as u
import numpy as np
import numpy.typing as npt
from astropy.units import Quantity

from dysh.util import get_project_data


class CalibratorTable:
    r"""
    This class is used to hold calibrator information.
    It uses as input a json file. By default it will look for
    dysh/data/calibrators.json, which is distributed with `dysh`.

    The structure of the input table should be:
    {
        "Objects" : {
            "name" : {
                "LAS"   : float,
                "cal"   : bool,
                "alias" : list
            },
        "Scale name" : {
                "fluxscale" : str,
                "objects" : {
                    "name" : list of coefficients
                },
                "method" : method used to compute the SED
            }
        }
    }

    Where "name" is the object name. It should match under the "Objects" and "Scale name" "objects".
    "LAS" is the largest angular extent of the object in degrees.
    "alias" is a list of alises for the object name.
    "Scale name" is the name of the scale, for example "Ott 1994".
    "method" is the method used to compute the spectral energy distribution given the list of
    coefficients provided for each object. The method should be defined in this file.

    Parameters
    ----------
    calibrator_table_file : Path
        Path to the json table with calibrator information.
        The contents of the file are defined above.
    """

    def __init__(self, calibrator_table_file: Path | None = None):
        if calibrator_table_file is None:
            self.calibrator_table_file = get_project_data() / "calibrators.json"
        else:
            self.calibrator_table_file = calibrator_table_file
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
            for k, v in self.data["Objects"]:
                names.append(k)
                names += v["alias"]
        else:
            for k in self.data[scale]["objects"].keys():
                names.append(k)
                names += self.data["Objects"][k]["alias"]
        return sorted(set(names), key=names.index)

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
    Holds calibrator details for a specific calibrator.

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
    las : float
        Largest angular size for the calibrator in degrees.
    cal : bool
        Is the calibrator suitable for flux density calibration?
        This is used to issue warnings.

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

    def __init__(
        self,
        name: str,
        scale: str,
        coefs: list,
        method: Callable,
        calibrator_table: CalibratorTable | None = None,
        las: float | None = None,
        cal: bool = True,
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

    @classmethod
    def from_name(
        cls,
        name: str,
        scale: Literal["Perley-Butler 2017", "Ott 1994"] = "Perley-Butler 2017",
        calibrator_table: CalibratorTable | None = None,
    ):
        """
        Create a `~dysh.util.calibrator.Calibrator` given a calibrator `name`.

        Parameters
        ----------
        name : str
            Calibrator name.
        scale : str or "Perley-Butler 2017" or "Ott 1994"
            The name of the flux scale to use. By defauilt `dysh` only knows
            of the Perley & Butler 2017 and Ott et al. 1994 scales.
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

        if scale not in calibrator_table.valid_scales():
            raise ValueError(f"Unrecognized scale. Valid scales: {', '.join(calibrator_table.valid_scales())}")
        if name not in calibrator_table.valid_names(scale):
            raise ValueError(
                f"Unrecognized calibrator name. Valid names: {', '.join(calibrator_table.valid_names(scale))}"
            )

        alias = calibrator_table.alias[name]
        coefs = calibrator_table.data[scale]["objects"][alias]
        method = globals()[calibrator_table.data[scale]["method"]]
        las = calibrator_table.data["Objects"][alias]["LAS"]
        cal = calibrator_table.data["Objects"][alias]["cal"]

        return cls(name, scale, coefs, method, calibrator_table, las=las, cal=cal)

    def compute_sed(self, nu: Quantity):
        """
        Evaluate the calibrator flux density at `nu`.

        Parameters
        ----------
        nu : `~astropy.units.Quantity`
            Frequency values.
        """

        # Calibrator flux density in Jy.
        snu = self.method(nu, self.coefs)

        return snu


def poly_pb(nu: Quantity, coefs: npt.ArrayLike):
    """
    Equation (1) in Perley & Butler 2017.

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
    """

    # Make sure the frequency is in GHz.
    nu = nu.to(u.GHz).value

    # Calibrator flux density in Jy.
    snu = np.power(10.0, np.polyval(coefs, np.log10(nu))) * u.Jy

    return snu


def poly_ott(nu: Quantity, coefs: list):
    """
    From Table 5 in Ott et al. 1994.

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
    """

    nu = nu.to(u.MHz).value
    lognu = np.log10(nu)

    snu = np.power(10.0, (coefs[0] + coefs[1] * lognu + coefs[2] * np.power(lognu, 2.0))) * u.Jy

    return snu
