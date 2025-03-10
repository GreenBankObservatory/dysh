"""This subpackage contains general utility classes and functions"""

__all__ = ["core", "selection", "download", "gaincorrection", "weatherforecast"]
from dysh.util.selection import Flag  # noqa:F401
from dysh.util.selection import Selection  # noqa:F401

from .core import *
