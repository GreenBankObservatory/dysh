"""Configuration"""

# https://docs.astropy.org/en/stable/config/index.html#customizing-config-location-in-affiliated-packages

import astropy.config as astropyconfig


class ConfigNamespace(astropyconfig.ConfigNamespace):
    rootname = "dysh"


class ConfigItem(astropyconfig.ConfigItem):
    rootname = "dysh"
