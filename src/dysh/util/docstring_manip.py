import textwrap

from astroquery.utils.docstr_chompers import remove_sections


def docstring_parameter(*args):
    """Decorator to pass variable value(s) into a docstring.

    Example
    -------
    @docstring_parameter('Ocean', 'Sea')
    def foo():
        \"\"\"My Docstring Lies Over The {0}.
        My Docstring Lies Over The {1}.\"\"\"
        pass

    print(foo.__doc__)  # Output: My Docstring Lies Over The Ocean. My Docstring Lies Over The Sea.
    """

    def dec(obj):
        obj.__doc__ = obj.__doc__.format(*args)
        return obj

    return dec


def append_docstr_nosections(doc, *, sections=["Returns"]):  # noqa B006
    """
    Decorator to append to the function's docstr after stripping out the
    list of sections provided (by default "Returns" only).
    """

    def dec(fn):
        fn.__doc__ = textwrap.dedent(fn.__doc__) + textwrap.indent(
            "\n".join(remove_sections(doc, sections)), prefix="\t"
        )
        # fn.__doc__ = fn.__doc__ + "\n".join(remove_sections(doc, sections))
        return fn

    return dec
