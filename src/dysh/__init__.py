"""Top-level package for dysh."""

import warnings

# Astropy Angle.to_string() triggers numpy vectorize floating-point warning (astropy#18989)
warnings.filterwarnings("ignore", message="invalid value encountered in do_format")
# ipympl multiple-inheritance MRO issue with traitlets (ipympl#488)
warnings.filterwarnings("ignore", message="Passing unrecognized arguments to super")

__version__ = "0.12.0"


all = ["version"]


def version():
    """Version of the dysh code

    Returns
    -------
    version : str
        dysh version.
    """
    return __version__
