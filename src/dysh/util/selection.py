import datetime
import numbers
import warnings
from collections.abc import Sequence
from copy import deepcopy

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import Angle
from astropy.table import Table
from astropy.time import Time
from astropy.units.quantity import Quantity
from pandas import DataFrame

from ..log import logger

# from ..fits import default_sdfits_columns
from . import ALL_CHANNELS, abbreviate_to, generate_tag, keycase

default_aliases = {
    "freq": "crval1",
    "ra": "crval2",
    "dec": "crval3",
    "glon": "crval2",
    "glat": "crval3",
    "gallon": "crval2",
    "gallat": "crval3",
    "elevation": "elevatio",
    "source": "object",
    "pol": "plnum",
    "subref": "subref_state",  # subreflector state
}


# workaround to avoid circular import error in sphinx (and only sphinx)
def _default_sdfits_columns():
    from ..fits import default_sdfits_columns

    return default_sdfits_columns()


DEFAULT_COLUMN_WIDTH = 32  # char
DEFAULT_COLUMN_TYPE = f"<U{DEFAULT_COLUMN_WIDTH}"


class SelectionBase(DataFrame):
    """This class is the base class for selection and flagging. Selection and flagging are both kinds of
    data selection, so `SelectionBase` can encapsulate most necessary functionality.  Derived classes
    implement specific named methods e.g. select_channel, flag_channel that will simply call
    the base class methods.
    """

    def __init__(self, initobj, aliases=default_aliases, **kwargs):
        if hasattr(initobj, "_index"):  # it's an SDFITSLoad object
            super().__init__(initobj._index)
            DEFKEYS = list(initobj._index.keys())
        else:
            super().__init__(initobj)  # it's a Selection or DataFrame
            DEFKEYS = _default_sdfits_columns()
        # adding attributes that are not columns will result
        # in a UserWarning, which we can safely ignore.
        warnings.simplefilter("ignore", category=UserWarning)
        self._add_datetime_column()
        self["CHAN"] = None
        # if we want Selection to replace _index in sdfits
        # construction this will have to change. if hasattr("_index") etc
        self._idtag = ["ID", "TAG"]
        # Add channel, timestamp, and number rows selected.
        DEFKEYS.extend(["FITSINDEX", "CHAN", "UTC", "# SELECTED"])
        # Add ID and TAG as the first columns.
        DEFKEYS = self._idtag + DEFKEYS
        # Remove duplicates.
        DEFKEYS = sorted(set(DEFKEYS), key=DEFKEYS.index)
        self._defkeys = DEFKEYS
        # set up object types for the np.array
        dt = np.full(len(DEFKEYS) - 1, np.dtype(DEFAULT_COLUMN_TYPE))
        dt[0] = np.int32
        # add number selected column which is an int
        dt = np.insert(dt, len(dt), np.int32)
        # ID is also an int
        dt[0] = np.int32
        self._deftypes = dt
        self._make_table()
        self._valid_coordinates = [
            "RA",
            "DEC",
            "GALLON",
            "GALLAT",
            "GLON",
            "GLAT",
            "CRVAL2",
            "CRVAL3",
        ]
        self._selection_rules = {}
        self._aliases = {}
        self.alias(aliases)
        self._channel_selection = None
        self._flag_channel_selection = {}  # used in Flag only
        warnings.resetwarnings()

    def _add_datetime_column(self):
        """
        Add column to the selection/flag dataframe with a
        representation of the SDFITS DATE-OBS which is a string,
        as an `~np.datetime64`.

        Returns
        -------
        None.

        """
        # Do not add utc=True to this call, as later comparisons will not work.

        self["UTC"] = pd.to_datetime(self["DATE-OBS"])

    def _make_table(self):
        """Create the table for displaying the selection rules"""
        self._table = Table(data=None, names=self._defkeys, dtype=self._deftypes)
        for t in self._idtag:
            self._table.add_index(t)

    @property
    def aliases(self):
        """
        The aliases that may be used to refer to SDFITS columns.

        Returns
        -------
        dict
            The dictionary of aliases and SDFITS column names
        """
        return self._aliases

    def alias(self, aliases):
        """
        Alias a set of keywords to existing columns. Multiple aliases for
        a single column are allowed, e.g.,
        { 'glon':'crval2', 'lon':'crval2'}

        Parameters
        ----------
        aliases : {}
            The dictionary of keywords and column names
            where the new alias is the key and
            the column name is the value and , i.e., {alias:column}

        Returns
        -------
        None.

        Raises
        ------
            ValueError if the column name is not recognized.
        """
        self._check_keys(aliases.values())
        for k, v in aliases.items():
            self._alias(k, v)

    def _alias(self, key, column):
        """
        Alias a new keyword to an existing column, e.g..
        to alias the SDFITS column 'CRVAL2' as 'RA':

            `alias('RA','CRVAL2')`

        The map is case insensitive, so `alias('ra', 'crval2')` also works.

        Parameters
        ----------
        key : str
            The new keyword to use as an alias.
        column : str
            The existing SDFITS column name to alias

        Returns
        -------
        None.

        """
        self._aliases[key.upper()] = column.upper()

    def _set_pprint_exclude_names(self):
        """Use `~astropy.Table.pprint_exclude_names` to set the list
        columns that have no entries.
        """
        if len(self._table) > 0:
            emptycols = np.array(self._table.colnames)[
                [np.all([self._table[k].data[i] == "" for i in range(len(self._table))]) for k in self._table.colnames]
            ]
            self._table.pprint_exclude_names.set(emptycols)

    def columns_selected(self):
        """The names of any columns which were used in a selection rule

        Returns
        -------
        colnames - set
            A set of str column names. An empty set is returned if no selection rule has yet been made.
        """
        if len(self._table) == 0:
            return set()

        self._set_pprint_exclude_names()  # ensure __attributes__ gets set.
        return (
            set(self._table.colnames)
            - set(self._table.meta["__attributes__"]["pprint_exclude_names"])
            - set(["# SELECTED", "ID", "TAG"])
        )

    def _sanitize_input(self, key, value):
        """
        Sanitize a key-value pair for.
        Coordinate and boolean types are checked for.

        Parameters
        ----------
        key : str
            Upper case key value.
        value : any
            The value for the key.

        Returns
        -------
        sanitized_value : str
            The sanitized value
        """
        # @todo   Allow minimum match str for key?
        if key in self._aliases.keys():
            key = self._aliases[key]
        if key not in self:
            raise KeyError(f"{key} is not a recognized column name.")
        v = self._sanitize_boolean(key, value)
        v = self._sanitize_coordinates(key, v)
        return v

    def _sanitize_boolean(self, key, value):
        """
        Sanitize a boolean selection key-value pair. Boolean values
        will be converted to "T" or "F" characters if the key is
        "SIG" or "CAL".

        Parameters
        ----------
        key : str
            Upper case key value.
        value : bool or any
            The value for the key.

        Returns
        -------
        sanitized_value : str
            The sanitized value. Either "T" or "F" if `value` is True or False, respectively.
            Otherwise, return the input `value`.
        """
        TF = {True: "T", False: "F"}
        bool_char_cols = ["SIG", "CAL"]
        if key in bool_char_cols and isinstance(value, bool):
            value = TF[value]
        return value

    def _sanitize_coordinates(self, key, value):
        """
        Sanitize a coordinate selection key-value pair. Coordinates will be
        converted to floats before the final value is created.

        Parameters
        ----------
        key : str
            Upper case key value.
        value : float or `~astropy.coordinates.Angle` or str or any
            The value for the key. It can be a single float,
            a single Angle (Quantity), a tuple of Angles
            (a1,a2,a3) or an Angle tuple, e.g., (n1,n2)*u.degree

        Returns
        -------
        sanitized_value : float or `~astropy.coordinates.Angle` or `value` type
            The sanitized value if it is a number, `~astropy.coordinates.Angle` or str.
            If it is not any of those, then return the input `value`.
        """
        if key not in self._valid_coordinates and key not in self.aliases:
            return value
        # note Quantity is derivative of np.ndarray, so
        # need to filter that out in the recursive call.
        # This is to handle (q1,q2) as a range.
        # (n1,n2)*u.degree is handled below
        if isinstance(value, (tuple, np.ndarray, list)) and not isinstance(value, Quantity):
            return [self._sanitize_coordinates(key, v) for v in value]
        if isinstance(value, numbers.Number):
            a = Angle(value * u.degree)
        else:  # it should be a str or Quantity
            a = Angle(value)
        return a.degree

    def _generate_tag(self, values, hashlen=9):
        """
        Generate a unique tag based on row values.  A hash object is
        created from the input values using SHA256, and a hex representation is created.
        The first `hashlen` characters of the hex string are returned.

        Parameters
        ----------
        values : array-like
            The values to use in creating the hash object
        hashlen : int, optional
            The length of the returned hash string. The default is 9.

        Returns
        -------
        tag : str
            The hash string

        """
        return generate_tag(values, hashlen)

    @property
    def _next_id(self) -> int:
        """
        Get the next ID number in the table.

        Returns
        -------
        id : int
            The highest existing ID number plus one
        """
        ls = len(self._table)
        if ls == 0:
            return 0
        return max(self._table["ID"]) + 1

    def _check_keys(self, keys):
        """
        Check a dictionary for unrecognized keywords.  This method is called in any select method to check inputs.

        Parameters
        ----------
        keys : list or array-like
           Keyword arguments

        Returns
        -------
        None.

        Raises
        ------
        KeyError
            If one or more keywords are unrecognized

        """
        # ignorekeys = ["PROPOSED_CHANNEL_RULE"]
        unrecognized = []
        ku = [k.upper() for k in keys]
        # for k in ignorekeys:
        #    if k in ku:
        #        ku.remove(k)
        for k in ku:
            if k not in self and k not in self._aliases:
                unrecognized.append(k)
        if len(unrecognized) > 0:
            raise KeyError(f"The following keywords were not recognized: {unrecognized}")

    def _check_numbers(self, **kwargs):
        self._check_type(
            numbers.Number,
            "Expected numeric value for these keywords but did not get a number",
            **kwargs,
        )

    def _check_range(self, **kwargs):
        bad = []
        badtime = []
        for k, v in kwargs.items():
            ku = k.upper()
            if not isinstance(v, (tuple, list, np.ndarray)):
                raise ValueError(f"Invalid input for key {ku}={v}. Range inputs must be tuple or list.")
            for a in v:
                if a is not None:
                    if isinstance(a, Quantity):
                        a = self._sanitize_coordinates(ku, a)
                    try:
                        if ku == "UTC":
                            badtime = self._check_type(np.datetime64, "Expected np.datetime64", silent=True, **{ku: a})
                        else:
                            self._check_numbers(**{ku: a})
                    except ValueError:
                        bad.append(ku)
        if len(bad) > 0 or len(badtime) > 0:
            msg = "Expected"
            a = " "
            if len(badtime) > 0:
                msg += f" np.datetime64 object for {badtime}"
                a = " and "
            if len(bad) > 0:
                msg += a
                msg += f"numeric value(s) for {bad} "
            msg += " but did not get that."
            raise ValueError(msg)

    def _check_type(self, reqtype, msg, silent=False, **kwargs):
        # @todo allow Quantities
        """
        Check that a list of keyword arguments is all a specified type.

        Parameters
        ----------

        reqtype : type
            The object type to check against, e.g. numbers.Number, str, etc

        msg : str
            The exception message to show if the inputs are not the specific reqtype

        **kwargs : dict or key=value
           Keyword arguments

        Raises
        ------
        ValueError
            If one or more of the values is not numeric.

        Returns
        -------
        None or, if silent is True,s a list of keywords that raised errors.

        """
        # deal with potential arrays first by calling
        # this method recursively on each array member.
        kw = deepcopy(kwargs)  # prevent concurrent modification
        recursive_bad = []
        for k, v in kwargs.items():
            # if the input value is an array, then we want to
            # check the type for each member of the array, but
            # not raise an exception (silent=True) but collect
            # the bad ones and pass them on.
            if isinstance(v, (Sequence, np.ndarray)) and not isinstance(v, str):
                bad = [self._check_type(reqtype, msg, **{k: x}, silent=True) for x in v]
                # there's probably a smarter way to do this
                for i in range(len(bad)):
                    if len(bad[i]) != 0:
                        recursive_bad.extend(bad[i])
                kw.pop(k)
        # now the main check
        ku = np.ma.masked_array([k.upper() for k in kw.keys()])
        ku.mask = np.array([isinstance(x, reqtype) for x in kw.values()])
        if len(recursive_bad) != 0:
            ku = np.ma.append(ku, recursive_bad)
        if silent:
            return list(ku[~ku.mask])
        if not np.all(ku.mask):
            raise ValueError(f"{msg}: {np.squeeze(ku[~ku.mask])}")

    def _check_for_duplicates(self, df):
        """
        Check that the user hasn't already added a rule matching this one

        Parameters
        ----------
        df : `~pandas.DataFrame`
            The selection to check


        Returns
        -------
        bool
           True if a duplicate was found, False if not.

        """
        for _id, s in self._selection_rules.items():
            tag = self._table.loc[_id]["TAG"]
            if s.equals(df):
                tag = self._table.loc[_id]["TAG"]
                warnings.warn(  # noqa: B028
                    f"A rule that results in an identical selection has already been added: ID: {_id}, TAG:{tag}."
                    " Ignoring."
                )
                return True
        return False

    def _addrow(self, row, dataframe, tag=None, check=False):
        """
        Common code to add a tagged row to the internal table after the selection has been created.
        Should be called in select* methods.

        Parameters
        ----------
        row : dict
            key, value pairs of the selection
        dataframe : `~pandas.DataFrame`
            The dataframe created by the selection.
        tag : str, optional
            An identifying tag by which the rule may be referred to later.
            If None, a  randomly generated tag will be created.
        check : bool
            If True, call `_check_for_duplicates()` to see if a dataframe already implements this rule.
        Returns
        -------
        None.

        """
        if check:
            if self._check_for_duplicates(dataframe):
                return
        if tag is not None:
            row["TAG"] = tag
        else:
            gentag = []
            # guarantee a unique seed by
            # including relevant key=value, which will be unique
            for k, v in row.items():
                if v is not None and v != "":
                    gentag.append(f"{k}={v}")
            row["TAG"] = self._generate_tag(gentag)
        row["ID"] = self._next_id
        row["# SELECTED"] = len(dataframe)
        self._selection_rules[row["ID"]] = dataframe
        for k, v in row.items():
            row[k] = abbreviate_to(DEFAULT_COLUMN_WIDTH, v)
        self._table.add_row(row)
        # for some reason the table gets "unsorted" from its index
        # resulting in issue #457
        # so always do a sort (by primary index by default) after adding a rows
        self._table.sort(self._idtag[0])

    def _replace_time(self, **kwargs):
        """Replace astropy.Time and datetime.datetime objects in a kwargs list with numpy.datetime64 equivalent.
        This is need because UTC is a datetime64 column but we want users to be able to input Time or datetime if desired

        Parameters
        ---------
        kwargs : dict
            dictionary of keywords/values

        Returns
        ------
            dict of updated values.
        """
        kc = kwargs.copy()
        for k, v in kc.items():
            if isinstance(v, Time):
                kc[k] = v.datetime64
            elif isinstance(v, datetime.datetime):
                kc[k] = np.datetime64(v)
            # @todo could probably do this with a clever recursive call to _replace_time.
            elif isinstance(v, (Sequence, np.ndarray)) and not isinstance(v, str):
                if isinstance(v, tuple):
                    v = list(v)
                for i in range(len(v)):
                    if isinstance(v[i], Time):
                        v[i] = v[i].datetime64
                    elif isinstance(v[i], datetime.datetime):
                        v[i] = np.datetime64(v[i])
                kc[k] = v
        return kc

    def _base_select(self, tag=None, check=False, **kwargs):
        """Add one or more exact selection/flag rules, e.g., `key1 = value1, key2 = value2, ...`
        If `value` is array-like then a match to any of the array members will be selected/flagged.
        Derived classes will call this method with their own specific name, i.e. `select` or `flag`.

        Parameters
        ----------
            tag : str
                An identifying tag by which the rule may be referred to later.
                If None, a  randomly generated tag will be created.
            check : bool
                If True, check that a previous selection does not give an identical result as this one.
            key : str
                The key  (SDFITS column name or other supported key)
            value : any
                The value to select

        Returns
        -------
            True if the selection resulted in a new rule, False if not (no data selected)

        """
        # pop these before check_keys, which is intended to check SDFITS keywords
        proposed_channel_rule = kwargs.pop("proposed_channel_rule", None)
        # if called via _select_from_mixed_kwargs, then we want to merge all the
        # selections
        df = kwargs.pop("startframe", self)
        self._check_keys(kwargs.keys())
        # While not necessary for adding a row to a Table, ensuring the dict
        # has keys for all Table columns improves the performance of Table._addrow.
        row = dict.fromkeys(self._table.colnames, "")

        single_value_queries = None
        multi_value_queries = None
        for k, v in list(kwargs.items()):
            if v is None:
                continue
            ku = k.upper()
            if ku in self._aliases:
                ku = self._aliases[ku]
            v = self._sanitize_input(ku, v)
            # If a list is passed in, it must be composed of strings.
            # Numeric lists are intepreted as ranges, so must be
            # selected by user with select_range
            if isinstance(v, (Sequence, np.ndarray)) and not isinstance(v, str):
                query = None
                for vv in v:
                    if ku == "UTC":
                        self._check_type(
                            np.datetime64,
                            "Expected np.datetime64 object but got something else.",
                            **{ku: vv},
                        )
                    # if it is a string, then OR them.
                    # e.g. object = ["NGC123", "NGC234"]
                    if isinstance(vv, str):
                        thisq = f'{ku} == "{vv}"'
                    else:
                        thisq = f"{ku} == {vv}"
                    if query is None:
                        query = thisq
                    else:
                        query += f"| {thisq}"
                    # for pd.merge to give the correct answer, we would
                    # need "inner" on the first one and "outer" on subsequent
                    # df = pd.merge(df, df[df[ku] == vv], how="inner")
                if multi_value_queries is None:
                    multi_value_queries = f"({query})"
                else:
                    multi_value_queries += f"&({query})"
            else:
                if isinstance(v, str):
                    thisq = f'{ku} == "{v}"'
                else:
                    thisq = f"{ku} == {v}"
                if single_value_queries is None:
                    single_value_queries = thisq
                else:
                    single_value_queries += f"& {thisq}"
            row[ku] = v
        if multi_value_queries is not None and single_value_queries is not None:
            query = f"{multi_value_queries} & {single_value_queries}"
        elif multi_value_queries is None and single_value_queries is not None:
            query = single_value_queries
        elif multi_value_queries is not None and single_value_queries is None:
            query = multi_value_queries
        else:
            warnings.warn("There was no data selection")  # should never happen  # noqa: B028
            return False
        df = df.query(query)
        if df.empty:
            warnings.warn("Your selection rule resulted in no data being selected. Ignoring.")  # noqa: B028
            return False
        df.loc[:, "CHAN"] = proposed_channel_rule  # this column is normally None so no need to check if None first.
        self._addrow(row, df, tag, check=check)
        return True

    def _base_select_range(self, tag=None, check=False, **kwargs):
        """
        Select a range of inclusive values for a given key(s).
        e.g., `key1 = (v1,v2), key2 = (v3,v4), ...`
        will select data  `v1 <= data1 <= v2, v3 <= data2 <= v4, ... `
        Upper and lower limits may be given by setting one of the tuple values
        to None. e.g., `key1 = (None,v1)` for an upper limit `data1 <= v1` and
        `key1 = (v1,None)` for a lower limit `data >=v1`.  Lower
        limits may also be specified by a one-element tuple `key1 = (v1,)`.

        Parameters
        ----------
        tag : str, optional
            An identifying tag by which the rule may be referred to later.
            If None, a  randomly generated tag will be created.
        check : bool
            If True, check that a previous selection does not give an identical result as this one.
        key : str
            The key (SDFITS column name or other supported key)
        value : array-like
            Tuple or list giving the lower and upper limits of the range.

        Returns
        -------
        None.

        """
        self._check_keys(kwargs.keys())
        kwargs.update(self._replace_time(**kwargs))
        self._check_range(**kwargs)
        row = {}
        df = self
        for k, v in list(kwargs.items()):
            ku = k.upper()
            if ku in self._aliases:
                ku = self._aliases[ku]
            v = self._sanitize_input(ku, v)
            # deal with a tuple quantity
            if isinstance(v, Quantity):
                v = v.value
            vn = []
            # deal with quantity inside a tuple.
            for q in v:
                # ultimately will need a map of
                # desired units, so e.g. if
                # GHz used, then the value is expressed in Hz
                if isinstance(q, Quantity):
                    vn.append(q.value)
                else:
                    vn.append(q)
            v = vn
            row[ku] = v
            if len(v) == 2:
                if v[0] is not None and v[1] is not None:
                    df = pd.merge(df, df[(df[ku] <= v[1]) & (df[ku] >= v[0])], how="inner")
                elif v[0] is None:  # upper limit given
                    df = pd.merge(df, df[(df[ku] <= v[1])], how="inner")
                else:  # lower limit given (v[1] is None)
                    df = pd.merge(df, df[(df[ku] >= v[0])], how="inner")
            elif len(v) == 1:  # lower limit given
                df = pd.merge(df, df[(df[ku] >= v[0])], how="inner")
            else:
                raise Exception(f"Couldn't parse value tuple {v} for key {k} as a range.")
        if df.empty:
            warnings.warn("Your selection rule resulted in no data being selected. Ignoring.")  # noqa: B028
            return
        self._addrow(row, df, tag)

    def _base_select_within(self, tag=None, check=False, **kwargs):
        """
        Select a value within a plus or minus for a given key(s).
        e.g. `key1 = [value1,epsilon1], key2 = [value2,epsilon2], ...`
        Will select data
        `value1-epsilon1 <= data1 <= value1+epsilon1,`
        `value2-epsilon2 <= data2 <= value2+epsilon2,...`

        Parameters
        ----------
        tag : str, optional
            An identifying tag by which the rule may be referred to later.
            If None, a  randomly generated tag will be created.
        check : bool
            If True, check that a previous selection does not give an identical result as this one.
        key : str
            The key (SDFITS column name or other supported key)
        value : array-like
            Tuple or list giving the value and epsilon

        Returns
        -------
        None.

        """
        # This is just a type of range selection.
        kw = {}
        for k, v in kwargs.items():
            v1 = v[0] - v[1]
            v2 = v[0] + v[1]
            kw[k] = (v1, v2)
        self._base_select_range(tag, **kw)

    def _base_select_channel(self, channel, tag=None):
        """
        Select channels and/or channel ranges. These are NOT used in :meth:`final`
        but rather will be used to create a mask for calibration or
        flagging. Single arrays/tuples will be treated as channel lists;
        nested arrays will be treated as  *inclusive* ranges. For instance:

        ``
        # select channel 24
        select_channel(24)
        # selects channels 1 and 10
        select_channel([1,10])
        # selects channels 1 thru 10 inclusive
        select_channel([[1,10]])
        # select channel ranges 1 thru 10 and 47 thru 56 inclusive, and channel 75
        select_channel([[1,10], [47,56], 75)])
        # tuples also work, though can be harder for a human to read
        select_channel(((1,10), [47,56], 75))
        ``

        *Note* : channel numbers start at zero.

        Parameters
        ----------
        channel : number, or array-like
            The channels to select
        Returns
        -------
        None.

        """
        # We don't want to get into trying to merge
        # different, possibly exclusive, channel selections.
        # This also avoids the side effect of using self to
        # compute "# Selected" in _addrow
        if self._channel_selection is not None:
            raise Exception(
                "You can only have one channel selection rule. Remove the old rule before creating a new one."
            )
        self._check_numbers(chan=channel)
        if isinstance(channel, numbers.Number):
            channel = [int(channel)]
        self._channel_selection = channel
        # we don't care if a selection selects the same channels.  They are all pasted together in numpy later and
        # never go through a DataFrame (which is why we pass in the dummy self)
        self._addrow(
            {"CHAN": abbreviate_to(DEFAULT_COLUMN_WIDTH, channel)},
            dataframe=self,
            tag=tag,
            check=False,
        )

    # NB: using ** in doc here because `id` will make a reference to the
    # python built-in function.  Arguably we should pick a different
    # keyword but 'id' is easy for user.
    def remove(self, id=None, tag=None):
        """Remove (delete) a selection rule(s).
        You must specify either **id** or **tag** but not both. If there are
        multiple rules with the same tag, they will all be deleted.

        Parameters
        ----------
            id : int
                The ID number of the rule as displayed in `show()`
            tag : str
                An identifying tag by which the rule may be referred to later.
        """
        if id is not None and tag is not None:
            raise Exception("You can only specify one of id or tag")
        if id is None and tag is None:
            raise Exception("You must specify either id or tag")
        if id is not None:
            if id in self._selection_rules:
                # We will assume that selection_rules and table
                # have been kept in sync.  The implementation
                # should ensure this.
                del self._selection_rules[id]
                row = self._table.loc_indices["ID", id]
                # there is only one row per ID
                self._table.remove_row(row)
            else:
                raise KeyError(f"No ID = {id} found in this Selection")
        else:
            # need to find IDs of selection rules where TAG == tag.

            # This will raise keyerror if tag not matched, so no need
            # to raise our own, unless we want to change the messgae.
            matching_indices = self._table.loc_indices["TAG", tag]
            #   raise KeyError(f"No TAG = {tag} found in this Selection")
            matching = Table(self._table[matching_indices])
            for i in matching["ID"]:
                del self._selection_rules[i]
                # self._selection_rules.pop(i, None) # also works
            self._table.remove_rows(matching_indices)

    def clear(self):
        """Remove all selection rules"""
        self._selection_rules = {}
        self._flag_channel_selection = {}
        self._make_table()
        self._channel_selection = None

    def show(self):
        """
        Print the current selection rules. Only columns with a rule are shown.
        The first two columns are ID number a TAG string. Either of these may be used
        to :meth:`remove` a row.  The final column `# SELECTED` gives
        the number of rows that a given rule selects from the original.
        The :meth:`final` selection may be fewer rows because each selection rule
        is logically OR'ed to create the final selection.

        Returns
        -------
        None.

        """
        self._set_pprint_exclude_names()
        print(self._table)

    @property
    def final(self):
        """
        Create the final selection. This is done by a logical AND of each
        of the selection rules (specifically `pandas.merge(how='inner')`).

        Returns
        -------
        final : DataFrame
            The resultant selection from all the rules.
        """
        return self.merge(how="inner")

    def merge(self, how, on=None):
        """
        Merge selection rules using a specific
        type of join.

        Parameters
        ----------
        how : {‘left’, ‘right’, ‘outer’, ‘inner’, ‘cross’}, no default.
            The type of join to be performed. See :meth:`pandas.merge()`.

        on: label or list
            Column or index level names to join on. These must be found in both DataFrames.
            If on is None and not merging on indexes then this defaults to the intersection
            of the columns in both DataFrames.

        Returns
        -------
        final : DataFrame
            The resultant selection from all the rules.

        """
        if len(self._selection_rules.values()) == 0:
            # warnings.warn("Selection.merge(): upselecting now")
            return DataFrame()
        final = None
        for df in self._selection_rules.values():
            if final is None:
                # need a deepcopy here in case there
                # is only one selection rule, because
                # we don't want to return a reference to the rule
                # which the receiver might modify.
                final = deepcopy(df)
            else:
                final = pd.merge(final, df, how=how, on=on)
        return final

    def _select_from_mixed_kwargs(self, **kwargs):
        """
        Called by calibration routines which may be mixing channel selections
        and exact selections, but **not** 'within' or 'range' selections.

        Parameters
        ----------
        **kwargs : dict
            Keyword arguments.  key=value as in the public selection methods.

        Returns
        -------
        None.

        """

        # get the tag if given or generate one if not
        kwlist = list(kwargs.items())
        tag = kwargs.pop("tag", self._generate_tag(kwlist))
        if len(kwargs) == 0:
            return  # user gave no additional kwargs
        if tag is None:  # in case user did tag=None (facepalm)
            tag = self._generate_tag(kwlist)
        logger.debug(f"working TAG IS {tag}")
        # in order to pop channel we need to check case insensitively
        ukwargs = keycase(kwargs)
        chan = ukwargs.pop("CHANNEL", None)
        if chan is not None:
            self._base_select_channel(chan, tag)
        if len(ukwargs) != 0:
            logger.debug(f"selection {ukwargs}")
            self._base_select(**ukwargs, tag=tag)

    def __deepcopy__(self, memo):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=UserWarning)
            cls = self.__class__
            result = cls.__new__(cls)
            memo[id(self)] = result
            for k, v in self.__dict__.items():
                setattr(result, k, deepcopy(v, memo))
            result._table = self._table.copy()
        return result

    def get(self, key):
        """Get the selection/flag rule by its ID

        Parameters
        ----------
        key : int
            The ID value.  See :meth:`show`.

        Returns
        -------
        `~pandas.DataFrame`
            The selection/flag rule
        """
        return self._selection_rules[key]


class Selection(SelectionBase):
    """This class contains the methods for creating rules to select data from an SDFITS object.
    Data (rows) can be selected using any column name in the input SDFITS object.
    Exact selection, range selection, upper/lower limit selection, and any-of selection
    are all supported.

    Users create *selection rules* by specifying keyword (SDFITS columns) and value(s) to be selected.
    Briefly, the selection methods are:

         :meth:`select` - Select exact values

         :meth:`select_range` - Select ranges of values

         :meth:`select_within` - Select a value +/- epsilon

         :meth:`select_channel` - Select channels or ranges of channels

    The Selection object maintains a DataFrame for each selection rule created by the user. The
    :meth:`final` selection is the logical OR of these rules. Users can examine the current selections
    with :meth:`show` which will show the current
    rules and how many rows each rule selects from the unfiltered data.


    Aliases of keywords are supported. The user may add an alias for an existing SDFITS column with :meth:`alias`.   Some default :meth:`aliases` have been defined.
    """

    def select(self, tag=None, check=False, **kwargs):
        """Add one or more exact selection rules, e.g., `key1 = value1, key2 = value2, ...`
        If `value` is array-like then a match to any of the array members will be selected.
        For instance `select(object=['3C273', 'NGC1234'])` will select data for either of those
        objects and `select(ifnum=[0,2])` will select IF number 0 or IF number 2.

        Parameters
        ----------
            tag : str
                An identifying tag by which the rule may be referred to later.
                If None, a  randomly generated tag will be created.
            check : bool
                If True, check that a previous selection does not give an identical result as this one.
            key : str
                The key  (SDFITS column name or other supported key)
            value : any
                The value to select

        """
        self._base_select(tag, check=check, **kwargs)

    def select_range(self, tag=None, **kwargs):
        """
        Select a range of inclusive values for a given key(s).
        e.g., `key1 = (v1,v2), key2 = (v3,v4), ...`
        will select data  `v1 <= data1 <= v2, v3 <= data2 <= v4, ... `
        Upper and lower limits may be given by setting one of the tuple values
        to None. e.g., `key1 = (None,v1)` for an upper limit `data1 <= v1` and
        `key1 = (v1,None)` for a lower limit `data >=v1`.  Lower
        limits may also be specified by a one-element tuple `key1 = (v1,)`.

        Parameters
        ----------
        tag : str, optional
            An identifying tag by which the rule may be referred to later.
            If None, a  randomly generated tag will be created.
        key : str
            The key (SDFITS column name or other supported key)
        value : array-like
            Tuple or list giving the lower and upper limits of the range.

        Returns
        -------
        None.

        """
        self._base_select_range(tag, **kwargs)

    def select_within(self, tag=None, **kwargs):
        """
        Select a value within a plus or minus for a given key(s).
        e.g. `key1 = [value1,epsilon1], key2 = [value2,epsilon2], ...`
        Will select data
        `value1-epsilon1 <= data1 <= value1+epsilon1,`
        `value2-epsilon2 <= data2 <= value2+epsilon2,...`

        Parameters
        ----------
        tag : str, optional
            An identifying tag by which the rule may be referred to later.
            If None, a  randomly generated tag will be created.
        key : str
            The key (SDFITS column name or other supported key)
        value : array-like
            Tuple or list giving the value and epsilon

        Returns
        -------
        None.

        """
        self._base_select_within(tag, **kwargs)

    def select_channel(self, channel, tag=None):
        """
        Select channels and/or channel ranges. These are NOT used in :meth:`final`
        but rather will be used to create a mask for calibration or
        flagging. Single arrays/tuples will be treated as channel lists;
        nested arrays will be treated as  *inclusive* ranges. For instance:

        ``
        # select channel 24
        select_channel(24)
        # selects channels 1 and 10
        select_channel([1,10])
        # selects channels 1 thru 10 inclusive
        select_channel([[1,10]])
        # select channel ranges 1 thru 10 and 47 thru 56 inclusive, and channel 75
        select_channel([[1,10], [47,56], 75)])
        # tuples also work, though can be harder for a human to read
        select_channel(((1,10), [47,56], 75))
        ``

        *Note* : channel numbers start at zero.

        Parameters
        ----------
        channel : number, or array-like
            The channels to select

        Returns
        -------
        None.

        """
        self._base_select_channel(channel, tag)


class Flag(SelectionBase):
    """This class contains the methods for creating rules to flag data from an SDFITS object.
    Data (rows) can be selected for flagging using any column name in the input SDFITS object.
    Exact selection, range selection, upper/lower limit selection, and any-of selection
    are all supported.

    Users create *flag rules* by specifying keyword (SDFITS columns) and value(s) to be flagged.
    Briefly, the flag methods are:

         :meth:`flag` - Flag exact values

         :meth:`flag_range` - Flag ranges of values

         :meth:`flag_within` - Flag a value +/- epsilon

         :meth:`flag_channel` - Flag channels or ranges of channels

    The Flag object maintains a DataFrame for each flag rule created by the user. The
    :meth:`final` flag is the logical OR of these rules. Users can examine the current flags
    with :meth:`show` which will show the current
    rules and how many rows each rule selects for flagging from the unfiltered data.

    The actual flags, which are per channel, are stored in the GBTFITSLoad object,
    not in the Flag object. The Flag object just contains the flagging rules.

    Aliases of keywords are supported. The user may add an alias for an existing SDFITS column with :meth:`alias`.   Some default :meth:`aliases` have been defined.

    GBTIDL Flags can be read in with :meth:`read`.
    """

    @property
    def final(self):
        """
        Create the final flag selection. This is done by a logical OR of each
        of the flag rules (specifically `pandas.merge(how='outer')`).
        Unlike Selection which uses AND logic to progressively narrow down data,
        Flag uses OR logic to cumulatively flag any data matching any rule.

        Returns
        -------
        final : DataFrame
            The resultant flagged rows from all the rules.
        """
        return self.merge(how="outer")

    def flag(self, tag=None, check=False, **kwargs):
        """Add one or more exact flag rules, e.g., `key1 = value1, key2 = value2, ...`
        If `value` is array-like then a match to any of the array members will be flagged.
        For instance `flag(object=['3C273', 'NGC1234'])` will select data for either of those
        objects and `flag(ifnum=[0,2])` will flag IF number 0 or IF number 2.  Channels for selected data
        can be flagged using keyword `channel`, e.g., `flag(object='MBM12',channel=[0,23])`
        will flag channels 0 through 23 *inclusive* for object MBM12.

        Parameters
        ----------
        tag : str
            An identifying tag by which the rule may be referred to later.
            If None, a  randomly generated tag will be created.
        check : bool
            If True, check that a previous selection does not give an identical result as this one.
        key : str
            The key  (SDFITS column name or other supported key)
        value : any
            The value to select

        """
        chan = kwargs.pop("channel", None)
        if chan is not None:
            if isinstance(chan, numbers.Number):
                chan = [int(chan)]
            self._check_numbers(chan=chan)
        if len(kwargs) == 0:
            # The user only entered channel as a keyword, so just call flag_channel
            self.flag_channel(channel=chan, tag=tag)
        else:
            # Select on the other kwargs then add channel to it.
            # Since we are allowing the behavior that the user can select
            # identical rows with different channel flags, we must
            # use a 'proposed channel rule' because the "CHAN" column is not normally set
            # before the _check_for_duplicates call inside _base_select.
            # The selection rules dataframes are allowed to be identical if the
            # the CHAN columns will be different.
            if chan is None:
                kwargs["proposed_channel_rule"] = ALL_CHANNELS
            else:
                kwargs["proposed_channel_rule"] = str(chan)
            success = self._base_select(tag, check=check, **kwargs)  # don't do this unless chan input is good.
            if not success:
                return
            idx = len(self._table) - 1
            if chan is not None:
                cc = abbreviate_to(DEFAULT_COLUMN_WIDTH, chan)
                self._table.loc[idx]["CHAN"] = cc
                self._flag_channel_selection[idx] = chan
                self._selection_rules[idx].loc[:, "CHAN"] = str(chan)
            else:
                self._flag_channel_selection[idx] = ALL_CHANNELS
                self._selection_rules[idx].loc[:, "CHAN"] = ALL_CHANNELS

    def flag_channel(self, channel, tag=None, **kwargs):
        """
        Flag  channels and/or channel ranges for *all data*. These are NOT used in :meth:`final`
        but rather will be used to create a mask for
        flagging. Single arrays/tuples will be treated as channel lists;
        nested arrays will be treated as  *inclusive* ranges. For instance:

        ```
        # flag channel 24
        flag_channel(24)
        # flag channels 1 and 10
        flag_channel([1,10])
        # flags channels 1 thru 10 inclusive
        flag_channel([[1,10]])
        # flag channel ranges 1 thru 10 and 47 thru 56 inclusive, and channel 75
        flag_channel([[1,10], [47,56], 75)])
        # tuples also work, though can be harder for a human to read
        flag_channel(((1,10), [47,56], 75))
        ```

        *Note* : channel numbers start at zero


         Parameters
        ----------
        channel : number, or array-like
            The channels to flag

        Returns
        -------
        None.

        """
        # okay to use base method because we are flagging all rows
        self._base_select_channel(channel, tag, **kwargs)
        idx = len(self._table) - 1
        self._flag_channel_selection[idx] = channel
        self._selection_rules[idx].loc[:, "CHAN"] = str(channel)
        self._channel_selection = None  # unused for flagging

    def flag_range(self, tag=None, check=False, **kwargs):
        """Flag a range of inclusive values for a given key(s).
        e.g., `key1 = (v1,v2), key2 = (v3,v4), ...`
        will flag data  `v1 <= data1 <= v2, v3 <= data2 <= v4, ... `
        Upper and lower limits may be given by setting one of the tuple values
        to None. e.g., `key1 = (None,v1)` for an upper limit `data1 <= v1` and
        `key1 = (v1,None)` for a lower limit `data >=v1`.  Lower
        limits may also be specified by a one-element tuple `key1 = (v1,)`.

        Parameters
        ----------
        tag : str, optional
            An identifying tag by which the rule may be referred to later.
            If None, a  randomly generated tag will be created.
        check : bool
            If True, check that a previous selection does not give an identical result as this one.
        key : str
            The key (SDFITS column name or other supported key)
        value : array-like
            Tuple or list giving the lower and upper limits of the range.

        Returns
        -------
        None.
        """
        self._base_select_range(tag, check=check, **kwargs)
        idx = len(self._table) - 1
        self._flag_channel_selection[idx] = ALL_CHANNELS
        self._selection_rules[idx].loc[:, "CHAN"] = ALL_CHANNELS
        self._channel_selection = None  # unused for flagging

    def flag_within(self, tag=None, check=False, **kwargs):
        """
        Flag a value within a plus or minus for a given key(s).
        e.g. `key1 = [value1,epsilon1], key2 = [value2,epsilon2], ...`
        Will select data
        `value1-epsilon1 <= data1 <= value1+epsilon1,`
        `value2-epsilon2 <= data2 <= value2+epsilon2,...`

        Parameters
        ----------
        tag : str, optional
            An identifying tag by which the rule may be referred to later.
            If None, a  randomly generated tag will be created.
        check : bool
            If True, check that a previous selection does not give an identical result as this one.
        key : str
            The key (SDFITS column name or other supported key)
        value : array-like
            Tuple or list giving the value and epsilon

        Returns
        -------
        None.

        """
        self._base_select_within(tag, check=check, **kwargs)
        idx = len(self._table) - 1
        self._flag_channel_selection[idx] = ALL_CHANNELS
        self._selection_rules[idx].loc[:, "CHAN"] = ALL_CHANNELS
        self._channel_selection = None  # unused for flagging

    def read(self, fileobj, ignore_vegas=False, **kwargs):
        """Read a GBTIDL flag file and instantiate Flag object.

        Parameters
        ----------
        fileobj : str, file-like or `pathlib.Path`
            File to read.  If a file object, must be opened in a
            readable mode.
        ignore_vegas : bool
            If True, ignore any flag rules which contain 'VEGAS_SPUR' in the line, as these
            are usually flagged via algorithm. See :meth:`~dysh.util.core.calc_vegas_spurs`.
        **kwargs : dict
            Extra keyword arguments to apply to the flag rule.  (This is mainly for internal use.)

        Returns
        -------
        None.

        """
        # GBTIDL flag files two sections [header] and [flags]
        # In the [header] section is information about file creation.
        # The [flags] section containes the flag table
        # The table has 10 columns. Its rows have vertical bar (|) separated columns, while
        # the table header is separated by commas and begins with a #
        # The columns are:
        #
        # ID - flag ID number, same as dysh's flag rule `id`
        # RECNUM - range of the selected record numbers given as low:high inclusive
        # SCAN - range of the selected scan numbers given as low:high inclusive
        # INTNUM - range of the selected integration numbers given as low:high inclusive
        # PLNUM - range of the selected polarization numbers given as low:high inclusive
        # IFNUM - range of the selected IF numbers given as low:high inclusive
        # BCHAN - beginning channel flagged (inclusive, starting from zero)
        # ECHAN - end channel flagged (inclusive)
        # IDSTRING - Reason for flagging, same as dysh's flag rule `tag`
        #
        # Numeric alues can be a single integer or comma-separated list of integers.  If BCHAN and ECHAN
        # are a comma-separated list then they must be pair up as [bchan_i,echan+i]
        # A wildcard appears in a column if it had no selection (meaning all values were selected).
        # Example file:
        # [header]
        # created = Wed Jan  5 16:48:37 2022
        # version = 1.0
        # created_by = sdfits
        # [flags]
        # #RECNUM,SCAN,INTNUM,PLNUM,IFNUM,FDNUM,BCHAN,ECHAN,IDSTRING
        # *|6|*|*|2|0|3072|3072|VEGAS_SPUR
        #
        # It is possible there is a space after the *GBTIDL flag files can also indicate ranges with a : and can indicate upper or lower limits
        # by not including a number. For instancer here is scan range 42 to 51 and channel range with
        # lower limit of 2299
        # *|20|42:51|*|*|*|2299|*|unspecified

        # Because the table header and table row delimeters are different,
        # Table.read() can't work.  So construct it row by row.
        f = open(fileobj)
        lines = f.read().splitlines()  # gets rid of \n
        f.close()
        header = [
            "RECNUM",
            "SCAN",
            "INTNUM",
            "PLNUM",
            "IFNUM",
            "FDNUM",
            "BCHAN",
            "ECHAN",
            "IDSTRING",
        ]
        found_header = False
        for l in lines[lines.index("[flags]") + 1 :]:
            # ignore VEGAS_SPUR flags if requested
            if ignore_vegas and l.strip().lower().endswith("vegas_spur"):
                continue
            vdict = {}
            if l.startswith("#"):
                if not found_header:
                    # its the header
                    colnames = l[1:].split(",")
                    if colnames != header:
                        raise Exception(f"Column names {colnames} do not match expectated {header}")
                    found_header = True
            else:
                values = l.split("|")
                for i, v in enumerate(values):
                    if v.strip() == "*":
                        continue
                    elif header[i] == "IDSTRING":
                        vdict[header[i]] = v
                    # handle comma-separated lists
                    elif "," in v:
                        vdict[header[i]] = [int(float(x)) for x in v.split(",")]
                    # handle colon-separated ranges by expanding into a comma-separated list.
                    elif ":" in v:
                        vdict[header[i]] = [int(float(x)) for x in range(*map(int, v.split(":")))] + [
                            int(v.split(":")[-1])
                        ]
                    # handle single values
                    else:
                        vdict[header[i]] = int(float(v))

                # our tag is gbtidl's idstring
                tag = vdict.pop("IDSTRING", None)
                bchan = vdict.pop("BCHAN", None)
                echan = vdict.pop("ECHAN", None)
                if bchan is not None and echan is not None:
                    if not isinstance(bchan, list):
                        bchan = [bchan]
                    bchan = [int(float(x)) for x in bchan]
                    if not isinstance(echan, list):
                        echan = [echan]
                    echan = [int(float(x)) for x in echan]
                    # pair up echan and bchan
                    vdict["channel"] = list(zip(bchan, echan, strict=False))
                elif bchan is not None and echan is None:
                    if not isinstance(bchan, list):
                        bchan = [bchan]
                    bchan = [int(float(x)) for x in bchan]
                    echan = [2**25] * len(
                        bchan
                    )  # Set to a large number so it effectively spans the whole range from `bchan`.
                    vdict["channel"] = tuple(zip(bchan, echan, strict=False))
                elif bchan is None and echan is not None:
                    if not isinstance(echan, list):
                        echan = [echan]
                    echan = [int(float(x)) for x in echan]
                    bchan = [0] * len(echan)
                    vdict["channel"] = tuple(zip(bchan, echan, strict=False))

                if kwargs is not None:
                    vdict.update(kwargs)
                logger.debug(f"flag({tag=},{vdict})")
                self.flag(tag=tag, check=False, **vdict)
            self._table.sort(self._idtag[0])
