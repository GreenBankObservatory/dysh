"""
Core plotting functions.
"""

from dysh.log import logger


def check_kwargs(known_kwargs, kwargs):
    """Check if `kwargs` are in `known_kwargs`"""
    diff = set(kwargs) - set(known_kwargs)
    if len(diff) > 0:
        logger.warning(f"Unknown kwargs: {', '.join(diff)}")


def parse_html(s):
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
    s = s.replace("<i>", "$").replace("</i>", "$")  # italics
    s = s.replace("&", "$\\").replace(";", "$")  # greek letters

    # strip everything else, maybe
    s = s.replace('<font color="red">', "").replace("</font>", "")
    s = s.replace("<font face=monospace>", "")

    s = s.replace("<b>", "").replace("</b>", "")
    s = s.replace(" (TopModel)", "")

    return s
