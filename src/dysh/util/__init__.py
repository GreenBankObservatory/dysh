"""This subpackage contains general utility classes and functions

isort:skip_file
"""

__all__ = ["core", "selection", "download", "gaincorrection", "weatherforecast"]
from .core import *  # this needs to be first to avoid circular imports
from dysh.util.selection import Flag  # noqa:F401
from dysh.util.selection import Selection  # noqa:F401
