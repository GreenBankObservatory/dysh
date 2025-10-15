"""
Methods for deling with non-ASCII strings
"""

import re

# From astropy.io.fits.Card:
# FSC commentary card string which must contain printable ASCII characters.
# Note: \Z matches the end of the string without allowing newlines
_ascii_text_re = re.compile(r"[ -~]*\Z")


def _ensure_ascii_str(text: str, check: bool = False) -> str:
    """does the actual cleaning of a text string"""
    clean_text = text.encode("ascii", "ignore").decode("ascii")
    clean_text = clean_text.replace("\n", " ")
    if check and _ascii_text_re.match(clean_text) is None:
        raise ValueError(f"Unable to fully clean string:{clean_text!r} of non-ASCII or non-printable characters.")

    return clean_text


def ensure_ascii(text: str | list[str], check: bool = False) -> str | list[str]:
    """
    Remove non-printable ASCII characters from a string or list of strings. This is to ensure that
    FITS cards conform to the standard

    Parameters
    ----------
    text : str
        The text to clean

    check: bool
        Check if the clean value is truly clean according to astropy FITS, raise ValueError if not
    Returns
    -------
    str or list[str]
        The cleaned text

    """
    if isinstance(text, str):
        return _ensure_ascii_str(text)
    else:
        clean_text = []
        for c in text:
            clean_text.append(_ensure_ascii_str(c))
        return clean_text
