import textwrap
from collections.abc import Callable
from typing import Any, TypeAlias, TypeVar

from astroquery.utils.docstr_chompers import remove_sections
from typing_extensions import ParamSpec

T = TypeVar("T")
P = ParamSpec("P")
WrappedFuncDeco: TypeAlias = Callable[[Callable[P, T]], Callable[P, T]]


def copy_docstring(copy_func: Callable[..., Any]) -> WrappedFuncDeco[P, T]:
    """Copies the doc string of the given function to another.
    This function is intended to be used as a decorator.

    .. code-block:: python3

        def foo():
            '''This is a foo doc string'''
            ...

        @copy_docstring(foo)
        def bar():
            ...
    """

    def dec(func: Callable[P, T]) -> Callable[P, T]:
        func.__doc__ = copy_func.__doc__
        return func

    return dec


def docstring_parameter(*args):
    """Decorator to pass variable value(s) into a docstring.

    Examples
    --------
    @docstring_parameter('Ocean', 'Sea')
    def foo():
        \"\"\"My Docstring Lies Over The {0}.
        My Docstring Lies Over The {1}.\"\"\"
        pass

    print(foo.__doc__)  # Output: My Docstring Lies Over The Ocean. My Docstring Lies Over The Sea.
    """

    def dec(obj):
        obj.__doc__ = textwrap.dedent(obj.__doc__).format(*args)
        return obj

    return dec


def keep_sections(doc, sections):
    """
    Given a numpy-formatted docstring, remove the section blocks provided in
    ``sections`` and dedent the whole thing.

    Returns
    -------
    List of lines
    """

    lines = iter(textwrap.dedent(doc).split("\n"))
    outlines = []
    keep_block = False
    for line in lines:
        lstrip = line.rstrip()
        if lstrip in sections:
            keep_block = True
            # Remove the section heading underlines.
            next(lines)
            continue
        elif keep_block:
            outlines.append(lstrip)
            if lstrip == "":
                keep_block = False
                continue
            else:
                continue

    return outlines


def append_docstr_nosections(doc, *, sections=None):
    """
    Decorator to append to the function's docstr after stripping out the
    list of sections provided (by default "Returns" only).
    """

    if sections is None:
        sections = ["Returns"]

    def dec(fn):
        fn.__doc__ = textwrap.dedent(fn.__doc__) + textwrap.indent(
            "\n".join(remove_sections(doc, sections)), prefix="\t"
        )
        return fn

    return dec


def append_docstr_sections(doc, *, sections=None, prefix="\t"):
    """
    Decorator to append to the function's docstr the list of
    `sections` provided from `doc` (by default "Parameters" only).

    Parameters
    ----------
    doc : str
        The docstr from which to extract the sections.
    sections : list
        The sections to keep from `doc`.
    """

    if sections is None:
        sections = ["Parameters"]

    def dec(fn):
        fn.__doc__ = textwrap.dedent(fn.__doc__) + textwrap.indent(
            "\n".join(keep_sections(doc, sections)), prefix=prefix
        )
        return fn

    return dec


def insert_docstr_section(doc, *, section="Parameters"):
    r"""
    Decorator to insert `section` from `doc`  into the function docstr.

    Examples
    --------
    >>> @insert_docstr_section("Parameters \n--------- \narg : str \nThis is a string.", section="Parameters")
    >>> def fun():
    >>>     \"\"\"This is a function.

    >>>     Parameters
    >>>     ----------
    >>>     {0}

    >>>     Returns
    >>>     -------
    >>>     Nothing
    >>>     \"\"\"
    >>>     pass

    >>> print(fun.__doc__)
    This is a function.

    Parameters
    ----------
    arg : str
        This is a string.

    Returns
    -------
    Nothing

    """

    def dec(obj):
        _doc = keep_sections(doc, sections=[section])
        _doc = (
            _doc[0]
            + "\n"
            + textwrap.indent("\n".join(_doc[1:-1]), prefix="\t")
            + textwrap.indent(_doc[-1], prefix="\t")
        )
        obj.__doc__ = obj.__doc__.format(_doc)
        return obj

    return dec
