"""
Core utility definitions, classes, and functions
"""

import hashlib
import importlib
import numbers
import sys
from collections.abc import Sequence
from itertools import zip_longest
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.table import Table
from astropy.time import Time
from astropy.units.quantity import Quantity
from IPython.display import HTML, display

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
    ary = df.ne(df.shift()).filter(items=[colname]).apply(lambda x: x.index[x].tolist()).to_numpy()
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


def to_mjd_list(time_val: Time | float) -> np.ndarray:
    """Convert an astropy Time, list of MJD, or single MJD to a list of MJD

    Parameters
    ----------
    time_val : `~astropy.time.Time` or float or list of float
        The time value to convert.

    Returns
    -------
    mjd : ~np.ndarray
        The Modified Julian Day values in an array. (or None if `time_val` was None)

    """
    if time_val is None:
        return None
    # check for Time first since it is also a Sequence
    if isinstance(time_val, Time):
        if time_val.isscalar:
            return np.array([time_val.mjd])
        else:
            return time_val.mjd
    if isinstance(time_val, (Sequence, np.ndarray)) and not isinstance(time_val, str):  # str is also a Sequence
        return time_val
    if isinstance(time_val, numbers.Number):
        return np.array([time_val])

    else:
        raise ValueError(f"Unrecognized type for time value: {type(time_val)}")


def to_quantity_list(q: Quantity | Sequence) -> Quantity:
    # if given quanity or [quanity], return [quanity.value]*quantity.units
    # handle quantities first
    if isinstance(q, Quantity):
        if q.isscalar:
            return [q.value] * q.unit
        else:
            return q
    # now handle lists of quantities
    if isinstance(q, Sequence):
        if len(set([x.unit for x in q])) != 1:
            raise ValueError("Units must all be the same in input list")
        return [x.value for x in q] * q[0].unit


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
    return importlib.resources.files("dysh")


def get_project_testdata() -> Path:
    """
    Returns the project testdata directory
    """
    return get_project_root().parent.parent / "testdata"


def get_project_data() -> Path:
    """
    Returns the directory where dysh configuration files are kept.

    Returns
    -------
    Path
        The project configuration directory.

    """
    return get_project_root() / "data"


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


def minimum_string_match(s, valid_strings, casefold=False):
    """
    return the valid string from a list, given a minimum string input

    Example:  minimum_string_match('a',['alpha','beta','gamma'])
    returns:  'alpha'

    Parameters
    ----------
    s : string
        string to use for minimum match
    valid_strings : list of strings
        list of full strings to minimum match on.
    casefold: bool
        If True, do a case insensitive match

    Returns
    -------
    string
        matched string, if one is found.  An exact match will
        also count as a match, even if others are present with
        longer match.
        Otherwise "None" is returned.

    """
    n = len(valid_strings)
    if casefold:
        vsfold = [a.casefold() for a in valid_strings]
        s = s.casefold()
    else:
        vsfold = valid_strings
    m = []
    for i in range(n):
        if vsfold[i].find(s) == 0:
            m.append(i)
    if len(m) >= 1:
        return valid_strings[m[0]]
    return None


def minimum_list_match(strings, valid_strings, casefold=False):
    """
    Return the list of valid strings given a list of minimum string inputs.

    Parameters
    ----------
    strings : str or list of str
        The strings to compare for minimum match
    valid_strings : list of str
        list of full strings to min match on.
    casefold: bool
        If True, do a case insensitive match

    Returns
    -------
    list
        List of all minimum matches or None if no matches found

    """
    valid = []
    # if user passes in a string instead of a list, it should act like minimum_string_match
    # Note: strings=list(strings) is not the same as [strings]!
    if isinstance(strings, str):
        strings = [strings]
    for s in strings:
        p = minimum_string_match(s, valid_strings, casefold)
        if p is not None:
            valid.append(p)
    if len(valid) == 0:
        return None
    else:
        return valid


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

    if str(a) == ALL_CHANNELS:
        return np.full(length, value)

    mask = np.full(length, False)
    for v in a:
        if isinstance(v, (tuple, list, np.ndarray)) and len(v) == 2:
            # If there are just two numbers, interpret is as an inclusive range
            mask[v[0] : v[1] + 1] = value
        else:
            mask[v] = value
    return mask


def abbreviate_to(length, value, squeeze=True) -> str:
    """
    Abbreviate a value for display in limited space. The abbreviated
    value will have initial characters, ellipsis, and final characters, e.g.
    '[(a,b),(c,d),...,(w,x),(y,z)]'.

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
        Abbreviated string representation of the input value.
    """
    strv = str(value)
    sep = ", "
    if squeeze:
        strv = strv.replace(", ", ",")
        sep = ","
    if len(strv) > length:
        try:
            bc = strv.rindex(sep, 0, length // 2 - 1)
        except ValueError:
            bc = strv.index(sep)
        try:
            ec = strv[-length // 2 + 1 :].index(sep)
            eci = -length // 2 + 1
        except ValueError:
            ec = strv.rindex(sep)
            eci = None
        strv = strv[:bc] + sep + "..." + strv[eci:][ec:]
    return strv


def merge_ranges(ranges):
    """
    Merge overlapping and adjacent ranges and yield the merged ranges
    in order. The argument must be an iterable of pairs (start, stop).

    Taken from: https://codereview.stackexchange.com/a/21333

    Parameters
    ----------
    ranges : iterable
        Pairs of (start, stop) ranges.

    Yields
    ------
    iterable
        Merged ranges.

    Examples
    --------
    >>> list(merge_ranges([(5,7), (3,5), (-1,3)]))
    [(-1, 7)]
    >>> list(merge_ranges([(5,6), (3,4), (1,2)]))
    [(1, 2), (3, 4), (5, 6)]
    >>> list(merge_ranges([]))
    []
    """
    ranges = iter(sorted(ranges))
    try:
        current_start, current_stop = next(ranges)
    except StopIteration:
        return
    for start, stop in ranges:
        if start > current_stop:
            # Gap between segments: output current segment and start a new one.
            yield current_start, current_stop
            current_start, current_stop = start, stop
        else:
            # Segments adjacent or overlapping: merge.
            current_stop = max(current_stop, stop)
    yield current_start, current_stop


def grouper(iterable, n, *, incomplete="fill", fillvalue=None):
    "Collect data into non-overlapping fixed-length chunks or blocks."
    # grouper('ABCDEFG', 3, fillvalue='x') → ABC DEF Gxx
    # grouper('ABCDEFG', 3, incomplete='strict') → ABC DEF ValueError
    # grouper('ABCDEFG', 3, incomplete='ignore') → ABC DEF
    iterators = [iter(iterable)] * n
    match incomplete:
        case "fill":
            return zip_longest(*iterators, fillvalue=fillvalue)
        case "strict":
            return zip(*iterators, strict=True)
        case "ignore":
            return zip(*iterators, strict=False)
        case _:
            raise ValueError("Expected fill, strict, or ignore")


def in_notebook() -> bool:
    """
    Check if the code is being run inside a notebook.
    """
    try:
        from IPython import get_ipython

        if "IPKernelApp" not in get_ipython().config:  # pragma: no cover
            return False
    except ImportError:
        return False
    except AttributeError:
        return False
    return True


def show_dataframe(df, show_index=False, max_rows=None, max_cols=None):
    """
    Function to show a `~pandas.DataFrame` in IPython or Jupyter.

    Parameters
    ----------
    df : `~pandas.DataFrame`
        The `~pandas.DataFrame` to be shown.
    show_index : bool
        Show the index of the `~pandas.DataFrame`.
    max_rows : int or None
        Maximum number of rows to display.
    max_cols : int or None
        Maximum number of columns to display.
    """

    kwargs = {"max_rows": max_rows, "max_cols": max_cols, "index": show_index}
    if in_notebook():
        display(HTML(df.to_html(**kwargs)))
    else:
        print(df.to_string(**kwargs))


def calc_vegas_spurs(
    vsprval: float | np.ndarray,
    vspdelt: float | np.ndarray,
    vsprpix: float | np.ndarray,
    maxchan: float,
    keep_central=False,
) -> np.ma.masked_array:
    """
    Calculate VEGAS spur channel locations.

    SPUR_CHANNEL = (J-VSPRVAL)*VSPDELT+VSPRPIX - 1

    where 0 <= J < 32.

    Spur channels are counted from zero.

    Parameters
    ----------
    vsprval : float or ~numpy.ndarray
        VEGAS spur channel offset
    vspdelt : float or ~numpy.ndarray
        VEGAS spur separation width in channels.
    vsprpix : float or ~numpy.ndarray
        VEGAS spur reference pixel.
    maxchan : float
        Maximum channel number (counting from zero), above which calculated spurs are masked.
    keep_central: bool
        Whether to keep the central VEGAS spur location in the returned array or not.
        The GBO SDFITS writer by default replaces the value at the central SPUR with the average of the
        two adjacent channels, and hence the central channel is not typically flagged.

    Note
    ----
    All input arrays must have the same shape

    Returns
    -------
    `~numpy.ma.masked_array`
        The array of channel numbers where spurs occur, with shape (`N_vsp`,31) where
        `N_vsp` is the length of the VSP arrays.  Invalid spur locations will be masked.

    """
    NSPURS = 32
    spurs = np.outer(np.arange(NSPURS), vspdelt) - vsprval * vspdelt + vsprpix - 1

    if not keep_central:
        # GBTIDL says: The central channel is defined by NCHAN/2 when counting from zero.
        # Which actually is ambiguous but we will take to mean e.g. 8192 when NCHAN=16384
        # This produces the same return array as GBTIDL dcspurschan.pro
        central = NSPURS // 2
        a = np.sort(spurs[np.arange(len(spurs)) != central].astype(int))
    else:
        a = np.sort(spurs.astype(int))
    # Mask any spur locations <0 or > maxchan.
    # Calling functions rely on the number of rows len(a[1]) to match the length of
    # the input VSP values, so mask instead of dropping rows.
    a = np.ma.masked_array(a, mask=np.logical_or(a < 0, a > maxchan))
    return a.T


def get_valid_channel_range(channel: list | np.ndarray) -> list:
    """
    Check that a channel range (e.g., that was given to Selection) defines a contiguous range of channels and
    return a list of length 2 if valid.  The returned list can be used in data calibration.
    Unlike :meth:`~dysh.util.selection.Selection.select_channel`,
    if the input `channel` list has
    two elements, it will be assumed that channel[0] is the first chan and channel[1] is the last channel.
    However, channel ranges specified identically to :meth:`~dysh.util.selection.Selection.select_channel` can also
    be used as input.

    Parameters
    ----------
    channel : list|np.ndarray
        List of beginning and end channel

    Raises
    ------
    ValueError
        If the input channel list cannot be converted to [first,last]

    Returns
    -------
    list
        An inclusive list of [first_channel, last_channel]
    """
    try:
        a = np.array(channel)
    except ValueError as ve:
        # it was some sort of compound list like [1,4,(30,40)]
        raise ValueError(f"Could not parse {channel} into a valid,contiguous channel range.") from ve
    ash = a.shape
    if ash == (1, 2):
        first = a[0][0]
        last = a[0][1]
    elif ash == (2,):
        first = a[0]
        last = a[1]
    if first > last:
        raise ValueError(f"In channel range {channel}, first channel is greater than last channel.")
    return [int(first), int(last)]


def isot_to_mjd(isot):
    """
    Convert an ISOT string to MJD.
    """
    EPOCH_MJD = 40587
    MSEC_IN_A_DAY = 86400e3
    return np.array(isot, dtype="datetime64[ms]").astype("int64") / MSEC_IN_A_DAY + EPOCH_MJD


def replace_col_astype(t: Table, colname: str, astype, fill_value):
    if hasattr(t[colname], "mask"):
        savemask = t[colname].mask.copy()
    else:
        savemask = False
    q = np.ma.masked_array(pd.to_numeric(t[colname]), savemask, fill_value=fill_value, dtype=astype)
    t[colname].fill_value = fill_value
    t.replace_column(colname, q)
    if hasattr(t, "mask"):
        t[colname].mask = savemask
