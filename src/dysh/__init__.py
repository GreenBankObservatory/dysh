"""Top-level package for dysh."""

from pathlib import Path

__version__ = "0.3.1"

all = ["version"]


def version():
    """Version of the dysh code

    :rtype: str
    """
    return __version__


def test_data_path():
    """
    Gives the path to the testdata directory

    Returns
    -------
    data_path : pathlib.Path
        The path to the testdata directory
    """
    data_path = Path(__file__).parent.parent.parent / "testdata"
    return data_path
