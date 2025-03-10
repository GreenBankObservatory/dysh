"""Configuration"""

# https://docs.astropy.org/en/stable/config/index.html#customizing-config-location-in-affiliated-packages

import astropy.config as astropyconfig

from dysh.config import core
from dysh.config.core import *


class ConfigNamespace(astropyconfig.ConfigNamespace):
    rootname = "dysh"


class ConfigItem(astropyconfig.ConfigItem):
    rootname = "dysh"


__all__ = [
    "core",
    "ConfigNamespace",
    "ConfigItem",
] + core.__all__
