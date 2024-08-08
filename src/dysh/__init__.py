"""Top-level package for dysh."""

from .log import init_logging

__version__ = "0.3.0b"

all = ["version"]

init_logging()


def version():
    """Version of the dysh code

    :rtype: str
    """
    return __version__
