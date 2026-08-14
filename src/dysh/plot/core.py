"""
Core plotting functions.
"""

from dysh.log import logger


def check_kwargs(known_kwargs, kwargs):
    """Check if `kwargs` are in `known_kwargs`"""
    diff = set(kwargs) - set(known_kwargs)
    if len(diff) > 0:
        logger.warning(f"Unknown kwargs: {', '.join(diff)}")


def catalog_html_to_latex(s):
    """Turn html-styled text from spectral line search to matplotlib mathtext"""

    # handle subscripts and superscripts
    s = s.replace("<sub>", "$_{")
    s = s.replace("</sub>", "}$")

    s = s.replace("<sup>", "$^{")
    s = s.replace("</sup>", "}$")

    # ge, le, etc.
    s = s.replace("le; ", "leq$")
    s = s.replace("ge; ", "geq$")

    # other formatting
    s = s.replace("<i>", "").replace("</i>", "")  # remove italics
    s = s.replace("&", "$\\").replace(";", "$")  # greek letters

    # strip everything else, maybe
    s = s.replace('<font color="red">', "").replace("</font>", "")
    # spacing not consistent in font labels, yeesh
    s = s.replace('<font color ="red">', "").replace("</font>", "")
    s = s.replace("<font face=monospace>", "")

    s = s.replace("<b>", "").replace("</b>", "")
    s = s.replace(" (TopModel)", "")
    # if there are multiple $ symbols, replace them with a single pair at
    # the beginning and end of the string
    # If there is space, which may be important, replace it with stretchable glue
    # (for instance " <sup>2</sup>" should be "~^2" because " ^2" fails latex)
    #
    count = sum(c == "$" for c in s)
    if count > 2:
        s = s.replace(" ", "~")
        s = s.replace("$", "")
        # ensure text is Roman font
        s = "$ {\\rm " + s + " }$"
    return s
