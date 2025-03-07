"""Classes and functions for importing SDFITS files"""

# Configuration options follow.
# These need to come after the global config variables, as some of the
# submodules use them.
from dysh import config as _config


class Conf(_config.ConfigNamespace):
    """
    Configuration parameters for `dysh.fits`.
    """

    summary_max_rows = _config.ConfigItem(
        1000,
        "Maximum number of rows to be displayed by summary",
    )


conf = Conf()

from dysh.fits.core import *
from dysh.fits.gb20mfitsload import GB20MFITSLoad  # noqa:F401
from dysh.fits.gbtfitsload import GBTFITSLoad  # noqa:F401
from dysh.fits.gbtfitsload import GBTOffline  # noqa:F401
from dysh.fits.gbtfitsload import GBTOnline  # noqa:F401
from dysh.fits.sdfitsload import SDFITSLoad  # noqa:F401

from . import core

__all__ = (
    ["Conf", "conf"]
    + core.__all__
    + [
        "GB20MFITSLoad",
        "GBTFITSLoad",
        "GBTOffline",
        "GBTOffline",
        "SDFITSLoad",
    ]
)
