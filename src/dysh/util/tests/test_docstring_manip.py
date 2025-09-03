"""
Tests for docstring_manip submodule.
"""

from dysh.util import docstring_manip as dm


def test_insert_docstr_section():
    # fmt: off
    @dm.insert_docstr_section("Parameters \n--------- \narg : str \nThis is a string.", section="Parameters")
    def fun():
        """
    This is a function.

    Parameters
    ----------
    {0}

    Returns
    -------
    Nothing
        """
        pass
    # fmt: on

    expected = "\n    This is a function.\n\n    Parameters\n    ----------\n    arg : str\n\tThis is a string.\n\n    Returns\n    -------\n    Nothing\n        "

    assert fun.__doc__ == expected
