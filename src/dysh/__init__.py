"""Top-level package for dysh."""

import warnings

# Astropy Angle.to_string() triggers numpy vectorize floating-point warning (astropy#18989)
warnings.filterwarnings("ignore", message="invalid value encountered in do_format")
# ipympl multiple-inheritance MRO issue with traitlets (ipympl#488)
warnings.filterwarnings("ignore", message="Passing unrecognized arguments to super")

# Pre-load the IERS-B table into astropy's internal cache once at import time.
# This prevents astropy from re-reading the table from disk on every
# EarthLocation.get_itrs() call during coordinate transforms.
from astropy.utils.iers import IERS_B
from astropy.utils.iers import conf as iers_conf

iers_conf.auto_download = False
IERS_B.open()

__version__ = "0.14.0"

all = ["version"]


def version():
    """Version of the dysh code

    Returns
    -------
    version : str
        dysh version.
    """
    return __version__
