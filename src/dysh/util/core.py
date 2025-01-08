"""
Core utility definitions, classes, and functions
"""

import hashlib
import re
import sys
from pathlib import Path
from typing import Union

# import astropy.units as u
import numpy as np

# import pandas as pd
from astropy.time import Time

ALL_CHANNELS = "all channels"


def select_from(key, value, df):
    """
    Select data where key=value.

    Parameters
    ----------
    key : str
        The key value (SDFITS column name)
    value : any
        The value to match
    df : `~pandas.DataFrame`
        The DataFrame to search

    Returns
    -------
    df : `~pandas.DataFrame`
        The subselected DataFrame

    """
    # nb this fails if value is None
    return df[(df[key] == value)]


def eliminate_flagged_rows(df, flag):
    """
    Remove rows from an index (selection) where all channels have been flagged.

    Parameters
    ----------
    df : `~pandas.DataFrame`
        The input dataframe from which flagged rows will be removed.
    flag : `~pandas.DataFrame`
        The flag dataframe.  Should be the result of e.g. `~util.Flag.final`

    Returns
    -------
        A data frame which is the input data frame with flagged rows removed.
    """
    if len(flag) > 0:
        # in the final flagging selection any rows that have CHAN=ALL_CHANNELS
        # indicate that the entire row is flagged
        ff = flag[flag["CHAN"].isin([ALL_CHANNELS])]
        flagged_rows = set(ff["ROW"])
        if len(flagged_rows) > 0:
            userows = list(set(df["ROW"]) - flagged_rows)
            if len(userows) > 0:
                return df[df["ROW"].isin(userows)]
            else:
                return df.iloc[0:0]  # all rows removed
    return df


def indices_where_value_changes(colname, df):
    """
    Find the `~pandas.DataFrame` indices where the value of the input column name changes.

    Parameters
    ----------
    colname : str
        The column name to query.
    df : `~pandas.DataFrame`
            The DataFrame to search

    Returns
    -------
    indices : ~numpy.ndarray
        The indices of the Dataframe where `colname` changes value.

    """
    # @todo add option to return changing values along with index
    # e.g., [["A",0],["B",125],["C",246]]
    # This is some super panda kung-fu.
    # See https://stackoverflow.com/questions/48673046/get-index-where-value-changes-in-pandas-dataframe-column
    if colname not in df:
        raise KeyError(f"Column {colname} not in input DataFrame")
    # df.shift() shifts the index by one, so we are then comparing df[N] to df[N-1]. This gets us
    # a truth table of where values change.  We filter on colname, then return a list of indices
    # where the value is true. Finally, we squeeze out the empty dimensions of the np array.
    ary = df.ne(df.shift()).filter(items=[colname]).apply(lambda x: x.index[x].tolist()).values
    return np.squeeze(ary, axis=1)


def gbt_timestamp_to_time(timestamp):
    """Convert the GBT sdfits timestamp string format to
    an :class:`~astropy.time.Time` object.  GBT SDFITS timestamps have the form
    YYYY_MM_DD_HH:MM:SS in UTC.

    Parameters
    ----------
    timestamp : str or list-like
        The GBT format timestamp as described above. If str, a Time object containing a single time is returned.
        If list-like, a Time object containing  multiple UTC times is returned.

    Returns
    -------
    time : `~astropy.time.Time`
        The time object
    """
    # convert to ISO FITS format  YYYY-MM-DDTHH:MM:SS(.SSS)
    if isinstance(timestamp, str):
        t = timestamp.replace("_", "-", 2).replace("_", "T")
    else:
        t = [ts.replace("_", "-", 2).replace("_", "T") for ts in timestamp]
    return Time(t, scale="utc")


def generate_tag(values, hashlen, add_time=True):
    """
    Generate a unique tag based on input values.  A hash object is
    created from the input values using SHA256, and a hex representation is created.
    The first `hashlen` characters of the hex string are returned.

    Parameters
    ----------
    values : array-like
        The values to use in creating the hash object
    hashlen : int, optional
        The length of the returned hash string.
    add_time: bool
        Add the time of the call to the values for hash generation.

    Returns
    -------
    tag : str
        The hash string

    """
    if add_time:
        values.append(Time.now().value)
    data = "".join(map(str, values))
    hash_object = hashlib.sha256(data.encode())
    unique_id = hash_object.hexdigest()
    return unique_id[0:hashlen]


def consecutive(data, stepsize=1):
    """Returns the indices of elements in `data`
    separated by less than stepsize separated into
    groups.

    Parameters
    ----------
    data : array
        Array with values to split.
    stepsize : int
        Maximum separation between elements of `data`
        to be considered a single group.

    Returns
    -------
    groups : `~numpy.ndarray`
        Array with values of `data` separated into groups.
    """
    return np.split(data, np.where(np.diff(data) >= stepsize)[0] + 1)


def sq_weighted_avg(a, axis=0, weights=None):
    # @todo make a generic moment or use scipy.stats.moment
    r"""Compute the mean square weighted average of an array (2nd moment).

    :math:`v = \sqrt{\frac{\sum_i{w_i~a_i^{2}}}{\sum_i{w_i}}}`

    Parameters
    ----------
    a : `~numpy.ndarray`
        The data to average
    axis : int
        The axis over which to average the data.  Default: 0
    weights : `~numpy.ndarray` or None
        The weights to use in averaging.  The weights array must be the
        length of the axis over which the average is taken.  Default:
        `None` will use equal weights.

    Returns
    -------
    average : `~numpy.ndarray`
        The average along the input axis
    """
    if weights is None:
        w = np.ones_like(a)
    else:
        w = weights
    v = np.sqrt(np.average(a * a, axis=axis, weights=w))
    return v


def get_project_root() -> Path:
    """
    Returns the project root directory.
    """
    return Path(__file__).parent.parent.parent.parent


def get_project_testdata() -> Path:
    """
    Returns the project testdata directory
    """
    return get_project_root() / "testdata"


def get_size(obj, seen=None):
    """Recursively finds size of objects.
    See https://goshippo.com/blog/measure-real-size-any-python-object/
    """
    size = sys.getsizeof(obj)
    if seen is None:
        seen = set()
    obj_id = id(obj)
    if obj_id in seen:
        return 0
    # Important mark as seen *before* entering recursion to gracefully handle
    # self-referential objects
    seen.add(obj_id)
    if isinstance(obj, dict):
        size += sum([get_size(v, seen) for v in obj.values()])
        size += sum([get_size(k, seen) for k in obj.keys()])
    elif hasattr(obj, "__dict__"):
        size += get_size(obj.__dict__, seen)
    elif hasattr(obj, "__iter__") and not isinstance(obj, (str, bytes, bytearray)):
        size += sum([get_size(i, seen) for i in obj])
    return size


def minimum_string_match(s, valid_strings):
    """
    return the valid string from a list, given a minimum string input

    Example:  minimum_string_match('a',['alpha','beta','gamma'])
    returns:  'alpha'

    Parameters
    ----------
    s : string
        string to use for minimum match
    valid_strings : list of strings
        list of full strings to min match on

    Returns
    -------
    string
        matched string, if one is found.
        Otherwise "None" is returned.

    """
    n = len(valid_strings)
    m = []
    for i in range(n):
        if valid_strings[i].find(s) == 0:
            m.append(i)
    if len(m) == 1:
        return valid_strings[m[0]]
    return None


def uniq(seq):
    """Remove duplicates from a list while preserving order.
    from http://stackoverflow.com/questions/480214/how-do-you-remove-duplicates-from-a-list-in-python-whilst-preserving-order
    """
    seen = set()
    seen_add = seen.add
    return [x for x in seq if x not in seen and not seen_add(x)]


def keycase(d, case="upper"):
    """
    Change the case of dictionary keys

    Parameters
    ----------
    d : dict
        The input dictionary
    case : str, one of 'upper', 'lower'
        Case to change keys to The default is "upper".

    Returns
    -------
    newDict : dict
        A copy of the dictionary with keys changed according to `case`

    """
    if case == "upper":
        newDict = {k.upper(): v for k, v in d.items()}
    elif case == "lower":
        newDict = {k.lower(): v for k, v in d.items()}
    return newDict


# Example of logging a function call
# @log_function_call(log_level="debug")
def powerof2(number):
    """
    Computes the closest power of 2 for a given `number`.

    Parameters
    ----------
    number : float
        number to determine the closest power of 2.

    Returns
    -------
    pow2 : int
        the closest power of 2.
    """

    return round(np.log10(number) / np.log10(2.0))


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


def ensure_ascii(text: Union[str, list[str]], check: bool = False) -> Union[str, list[str]]:
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


def convert_array_to_mask(a, length, value=True):
    """
    This method interprets a simple or compound array and returns a numpy mask
    of length `length`. Single arrays/tuples will be treated as element index lists;
    nested arrays will be treated as *inclusive* ranges, for instance:

    ``
    # mask elements 1 and 10
    convert_array_to_mask([1,10])
    # mask elements 1 thru 10 inclusive
    convert_array_to_mask([[1,10]])
    # mask ranges 1 thru 10 and 47 thru 56 inclusive, and element 75
    convert_array_to_mask([[1,10], [47,56], 75)])
    # tuples also work, though can be harder for a human to read
    convert_array_to_mask(((1,10), [47,56], 75))
    ``

    Parameters
    ----------
    a : number or array-like
        The
    length : int
        The length of the mask to return, e.g. the number of channels in a spectrum.

    value : bool
        The value to fill the mask with.  True to mask data, False to unmask.

    Returns
    -------
    mask : ~np.ndarray
        A numpy array where the mask is True according to the rules above.

    """

    if a == ALL_CHANNELS:
        return np.full(length, value)

    mask = np.full(length, False)

    for v in a:
        if isinstance(v, (tuple, list, np.ndarray)) and len(v) == 2:
            # If there are just two numbers, interpret is as an inclusive range
            mask[v[0] : v[1] + 1] = value
        else:
            mask[v] = value
    return mask


def abbreviate_to(length, value, squeeze=True):
    """
    Abbreviate a value for display in limited space. The abbreviated
    value will have initial characters, ellipsis, and final characters, e.g.
    '[(a,b),(c,d)...(w,x),(y,z)]'.

    Parameters
    ----------
    length : int
        Maximum string length.
    value : any
        The value to be abbreviated.
    squeeze : bool, optional
        Squeeze blanks. If True, replace ", " (comma space) with "," (comma). The default is True.

    Returns
    -------
    strv : str
        Abbreviated string representation of the input value

    """
    strv = str(value)
    if squeeze:
        strv = strv.replace(", ", ",")
    if len(strv) > length:
        bc = int(length / 2) - 1
        ec = bc - 1
        strv = strv[0:bc] + "..." + strv[-ec:]
    return strv
