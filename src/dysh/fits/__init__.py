"""Classes and functions for importing SDFITS files."""

from importlib import import_module

try:
    import fitsio  # noqa: F401

    HAS_FITSIO = True
except ImportError:
    HAS_FITSIO = False

from dysh import config as _config

from . import core


class Conf(_config.ConfigNamespace):
    """
    Configuration parameters for `dysh.fits`.
    """

    summary_max_rows = _config.ConfigItem(
        1000,
        "Maximum number of rows to be displayed by summary",
    )


conf = Conf()

_EXPORTS = {
    "GB20MFITSLoad": (".gb20mfitsload", "GB20MFITSLoad"),
    "GBTBackend": (".gbtfitsload", "GBTBackend"),
    "GBTFITSLoad": (".gbtfitsload", "GBTFITSLoad"),
    "GBTOffline": (".gbtfitsload", "GBTOffline"),
    "GBTOnline": (".gbtfitsload", "GBTOnline"),
    "LazyFlagArray": (".lazyflag", "LazyFlagArray"),
    "LazyFlagContainer": (".lazyflag", "LazyFlagContainer"),
    "SDFITSLoad": (".sdfitsload", "SDFITSLoad"),
}

_MODULE_EXPORTS = {
    "core": ".core",
    "gb20mfitsload": ".gb20mfitsload",
    "gbtfitsload": ".gbtfitsload",
    "lazyflag": ".lazyflag",
    "sdfitsload": ".sdfitsload",
}

__all__ = ["Conf", "conf", "HAS_FITSIO"] + core.__all__ + list(_EXPORTS) + list(_MODULE_EXPORTS)


def __getattr__(name):
    if name in core.__all__:
        return getattr(core, name)
    if name in _EXPORTS:
        module_name, attr_name = _EXPORTS[name]
        module = import_module(module_name, __name__)
        value = getattr(module, attr_name)
        globals()[name] = value
        return value
    if name in _MODULE_EXPORTS:
        module = import_module(_MODULE_EXPORTS[name], __name__)
        globals()[name] = module
        return module
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
