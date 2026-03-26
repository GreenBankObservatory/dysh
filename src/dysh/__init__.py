"""Top-level package for dysh."""

import warnings

# Astropy Angle.to_string() triggers numpy vectorize floating-point warning (astropy#18989)
warnings.filterwarnings("ignore", message="invalid value encountered in do_format")
# ipympl multiple-inheritance MRO issue with traitlets (ipympl#488)
warnings.filterwarnings("ignore", message="Passing unrecognized arguments to super")

__version__ = "1.1.0"

__all__ = ["version", "system_info"]


def version():
    """Version of the dysh code

    Returns
    -------
    version : str
        dysh version.
    """
    return __version__


from dysh.util.system_info import system_info
