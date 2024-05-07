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

# from ..fits import default_sdfits_columns
from . import gbt_timestamp_to_time, generate_tag, keycase

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
    "intnum": "row",  # integration number
    "subref": "subref_state",  # subreflector state
}


# workaround to avoid circular import error in sphinx (and only sphinx)
def _default_sdfits_columns():
    from ..fits import default_sdfits_columns

    return default_sdfits_columns()


class Selection(DataFrame):
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

    def __init__(self, initobj, aliases=default_aliases, **kwargs):
        if hasattr(initobj, "_index"):  # it's an SDFITSLoad object
            super().__init__(initobj._index, copy=True)
            DEFKEYS = list(initobj._index.keys())
        else:
            super().__init__(initobj, copy=True)  # it's a Selection or DataFrame
            DEFKEYS = _default_sdfits_columns()
        # adding attributes that are not columns will result
        # in a UserWarning, which we can safely ignore.
        warnings.simplefilter("ignore", category=UserWarning)
        self._add_utc_column()
        # if we want Selection to replace _index in sdfits
        # construction this will have to change. if hasattr("_index") etc
        self._idtag = ["ID", "TAG"]
        DEFKEYS.extend(["CHAN", "UTC", "# SELECTED"])
        # add ID and TAG as the first columns
        for i in range(len(self._idtag)):
            DEFKEYS.insert(i, self._idtag[i])
        # add channel, astropy-based timestamp, and number rows selected
        DEFKEYS = np.array(DEFKEYS)
        # set up object types for the np.array
        dt = np.array([str] * (len(DEFKEYS) - 1))
        # add number selected column which is an int
        dt = np.insert(dt, len(dt), np.int32)
        # ID is also an int
        dt[0] = np.int32
        self._defkeys = DEFKEYS
        self._deftypes = dt
        self._make_table()
        self._valid_coordinates = ["RA", "DEC", "GALLON", "GALLAT", "GLON", "GLAT", "CRVAL2", "CRVAL3"]
        self._selection_rules = {}
        self._aliases = {}
        self.alias(aliases)
        self._channel_selection = None
        warnings.resetwarnings()

    def _add_utc_column(self):
        """
        Add column to the selection dataframe with a
        representation of the SDFITS UTC timestamp, which is a string,
        as an ~astropy.time.Time.

        Returns
        -------
        None.
        """
        self["UTC"] = [gbt_timestamp_to_time(q) for q in self.TIMESTAMP]

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

    def _sanitize_input(self, key, value):
        """
        Sanitize a key-value pair for. Coordinate types are checked for.

        Parameters
        ----------
        key : str
            upper case key value

        value : any
            The value for the key

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
        v = self._sanitize_coordinates(key, value)
        # deal with Time here or later?
        self._check_for_disallowed_chars(key, value)
        return v

    def _sanitize_coordinates(self, key, value):
        """
        Sanitize a coordinate selection key-value pair. Coordinates will be
        converted to floats before the final value is created.

        Parameters
        ----------
        key : str
            upper case key value

        value : any
            The value for the key.  It can be a single float,
            a single Angle (Quantity), a tuple of Angles
            (a1,a2,a3) or an Angle tuple, e.g., (n1,n2)*u.degree

        Returns
        -------
        sanitized_value : str
            The sanitized value.
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

    def _check_for_disallowed_chars(self, key, value):
        # are there any?  coordinates will already
        # be transformed to decimal degrees
        pass

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
        unrecognized = []
        ku = [k.upper() for k in keys]
        for k in ku:
            if k not in self and k not in self._aliases:
                unrecognized.append(k)
        # print("KU, K", ku, k)
        if len(unrecognized) > 0:
            raise KeyError(f"The following keywords were not recognized: {unrecognized}")

    def _check_numbers(self, **kwargs):
        self._check_type(numbers.Number, "Expected numeric value for these keywords but did not get a number", **kwargs)

    def _check_range(self, **kwargs):
        bad = []
        badtime = []
        for k, v in kwargs.items():
            ku = k.upper()
            # print(ku)
            if not isinstance(v, (tuple, list, np.ndarray)):
                raise ValueError(f"Invalid input for key {ku}={v}. Range inputs must be tuple or list.")
            for a in v:
                if a is not None:
                    if isinstance(a, Quantity):
                        a = self._sanitize_coordinates(ku, a)
                    try:
                        if ku == "UTC":
                            badtime = self._check_type(Time, "Expected Time", silent=True, **{ku: a})
                            # print("BADTIME ", badtime)
                        else:
                            self._check_numbers(**{ku: a})
                    except ValueError:
                        bad.append(ku)
        if len(bad) > 0 or len(badtime) > 0:
            msg = "Expected"
            a = " "
            if len(badtime) > 0:
                msg += f" Time object for {badtime}"
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
        df : ~pandas.DataFrame
            The selection to check


        Returns
        -------
        bool
           True if a duplicate was found, False if not.

        """
        #        Raises
        #        ------
        #        Exception
        #            If an identical rule (DataFrame) has already been added.
        for _id, s in self._selection_rules.items():
            if s.equals(df):
                # print(s, df)
                tag = self._table.loc[_id]["TAG"]
                # raise Exception(
                warnings.warn(
                    f"A rule that results in an identical selection has already been added: ID: {_id}, TAG:{tag}."
                    " Ignoring."
                )
                return True
                # )
        return False

    def _addrow(self, row, dataframe, tag=None):
        """
        Common code to add a tagged row to the internal table after the selection has been created.
        Should be called in select* methods.

        Parameters
        ----------
        row : dict
            key, value pairs of the selection
        dataframe : ~pandas.DataFrame
            The dataframe created by the selection.
        tag : str, optional
            An identifying tag by which the rule may be referred to later.
            If None, a  randomly generated tag will be created.
        Returns
        -------
        None.

        """
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
        self._table.add_row(row)

    def select(self, tag=None, **kwargs):
        """Add one or more exact selection rules, e.g., `key1 = value1, key2 = value2, ...`
        If `value` is array-like then a match to any of the array members will be selected.
        For instance `select(object=['3C273', 'NGC1234'])` will select data for either of those
        objects and `select(ifnum=[0,2])` will select IF number 0 or IF number 2.

        Parameters
        ----------
            tag : str
                An identifying tag by which the rule may be referred to later.
                If None, a  randomly generated tag will be created.
            key : str
                The key  (SDFITS column name or other supported key)
            value : any
                The value to select

        """
        # @todo ?? MAYBE allow chan(nel) in here, e.g.
        # chan = kwargs.pop(chan,None)
        # if chan is not None:
        #   self.select_channel(chan,tag=tag)
        #
        self._check_keys(kwargs.keys())
        row = {}
        # if called via _select_from_mixed_kwargs, then we want to merge all the
        # selections
        df = kwargs.pop("startframe", self)
        for k, v in list(kwargs.items()):
            ku = k.upper()
            if ku in self._aliases:
                ku = self._aliases[ku]
            v = self._sanitize_input(ku, v)
            # If a list is passed in, it must be composed of strings.
            # Numeric lists are intepreted as ranges, so must be
            # selected by user with select_range
            if isinstance(v, (Sequence, np.ndarray)) and not isinstance(v, str):
                # print(ku, v)
                query = None
                for vv in v:
                    if ku == "UTC":
                        self._check_type(
                            Time,
                            "Expected Time object but got something else.",
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
                # print("final query ", query)
                df = df.query(query)
            else:
                df = pd.merge(df, df[df[ku] == v], how="inner")
            row[ku] = str(v)
        if df.empty:
            warnings.warn("Your selection rule resulted in no data being selected. Ignoring.")
            return
        self._addrow(row, df, tag)
        # return df

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
        # @todo ?? MAYBE allow chan(nel) in here, e.g.
        # chan = kwargs.pop(chan,None)
        # if chan is not None:
        #   self.select_channel(chan,tag=tag)
        self._check_keys(kwargs.keys())
        self._check_range(**kwargs)
        row = {}
        df = self
        for k, v in list(kwargs.items()):
            ku = k.upper()
            if ku in self._aliases:
                ku = self._aliases[ku]
            v = self._sanitize_input(ku, v)
            # print(f"{ku}={v}")
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
            row[ku] = str(v)
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
            warnings.warn("Your selection rule resulted in no data being selected. Ignoring.")
            return
        self._addrow(row, df, tag)

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
        # This is just a type of range selection.
        kw = {}
        for k, v in kwargs.items():
            v1 = v[0] - v[1]
            v2 = v[0] + v[1]
            kw[k] = (v1, v2)
        self.select_range(tag, **kw)

    def select_channel(self, chan, tag=None):
        """
        Select channels and/or channel ranges. These are NOT used in :meth:`final`
        but rather will be used to create a mask for calibration or
        flagging. Single arrays/tuples will be treated as channel lists;
        nested arrays will be treated as ranges, for instance

        ``
        # selects channels 1 and 10
        select_channel([1,10])
        # selects channels 1 thru 10 inclusive
        select_channel([[1,10]])
        # select channel ranges 1 thru 10 and 47 thru 56 inclusive, and channel 75
        select_channel([[1,10], [47,56], 75)])
        # tuples also work, though can be harder for a human to read
        select_channel(((1,10), [47,56], 75))
        ``

        Parameters
        ----------
        chan : number, or array-like
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
        self._check_numbers(chan=chan)
        self._channel_selection = chan
        self._addrow({"CHAN": str(chan)}, dataframe=self, tag=tag)

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
        self._make_table()

    def show(self):
        """
        Print the current selection rules. Only columns with a rule are shown.
        The first two columns are ID number a TAG string. Either of these may be used
        to :meth:remove a row.  The final column `# SELECTED` gives
        the number of rows that a given rule selects from the original.
        The :meth:final selection may be fewer rows because each selection rule
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
        tag = kwargs.pop("tag", self._generate_tag(kwargs))
        debug = kwargs.pop("debug", False)
        if len(kwargs) == 0:
            return  # user gave no additional kwargs
        if tag is None:  # in case user did tag=None (facepalm)
            tag = self._generate_tag(kwargs)
        if debug:
            print(f"working TAG IS {tag}")
        # in order to pop channel we need to check case insensitively
        ukwargs = keycase(kwargs)
        chan = ukwargs.pop("CHANNEL", None)
        if chan is not None:
            self.select_channel(chan, tag)
        if len(ukwargs) != 0:
            if debug:
                print(f"selection {ukwargs}")
            self.select(**ukwargs, tag=tag)

    def __deepcopy__(self, memo):
        warnings.simplefilter("ignore", category=UserWarning)
        cls = self.__class__
        result = cls.__new__(cls)
        memo[id(self)] = result
        for k, v in self.__dict__.items():
            setattr(result, k, deepcopy(v, memo))
        warnings.resetwarnings()
        return result
