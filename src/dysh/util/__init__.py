"""This subpackage contains general utility classes and functions"""

__all__ = ["core", "selection", "download", "gaincorrection", "timers", "weatherforecast"]
from .core import *  # this needs to be first to avoid circular imports
from .docstring_manip import (
    docstring_parameter,
    append_docstr_nosections,
    copy_docstring,
    insert_docstr_section,
    append_docstr_sections,
)
from dysh.util.selection import Flag  # noqa:F401
from dysh.util.selection import Selection  # noqa:F401
