import numbers

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import Angle
from astropy.table import Table  # , TableAttribute
from pandas import DataFrame

# from ..fits import default_sdfits_columns
from . import generate_tag

# add ID and TAG to available keys before creating tag
# why not move this to constructor?
idtag = ["ID", "TAG"]
# DEFKEYS = np.array(default_sdfits_columns())
# DEFKEYS = np.insert(DEFKEYS, 0, idtag)


class Selection(DataFrame):
    def __init__(self, sdfits, **kwargs):
        super().__init__(sdfits._index, copy=True)
        DEFKEYS = np.array(sdfits._index.keys())
        DEFKEYS = np.insert(DEFKEYS, 0, idtag)
        dt = np.array([str] * len(DEFKEYS))
        dt[0] = int
        self._table = Table(data=None, names=DEFKEYS, dtype=dt)
        for t in idtag:
            self._table.add_index(t)
        self._valid_coordinates = ["RA", "DEC", "GALLON", "GALLAT"]
        self._selection_rules = dict()

    def _set_pprint_exclude_names(self):
        """Use pprint_exclude_names to set the list
        columns that have no entries.
        """
        if len(self._table) > 0:
            emptycols = np.array(self._table.colnames)[
                [np.all([self._table[k].data[i] == "" for i in range(len(self._table))]) for k in self._table.colnames]
            ]
            self._table.pprint_exclude_names.set(emptycols)

    # def __repr__(self):
    # when printing to screen we only want to include
    # columns that are not empty.
    # This will not affect writing to a file.
    #     self._set_pprint_exclude_names()
    #     return super().__repr__()

    def _sanitize_input(self, key, value):
        """
        Sanitize a key-value pair for. List and coordinate types are checked for.

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
        # @todo 1.  Allow minimum match str for key
        #      2.  Allow synonyms, e.g. source for object,
        #         elevation for elevatio, etc.
        if key not in self:
            raise KeyError(f"{key} is not a recognized column name.")
        # v = self._sanitize_list(value)
        v = self._sanitize_coordinates(key, value)
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
            The value for the key

        Returns
        -------
        sanitized_value : str
            The sanitized value.
        """
        if key not in self._valid_coordinates:
            return value
        if isinstance(value, float):
            a = Angle(value * u.degree)
        else:  # it should be a str or Quantity
            a = Angle(value)
        return a.degree

    def _check_for_disallowed_chars(self, key, value):
        # are there any?  coordinates will already
        # be transformed to decimal degrees
        pass

    def _parse_time(self, value):
        """


        Parameters
        ----------
        value : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
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

    def _check_keys(self, **kwargs):
        """
        Check a list for unrecognized keywords.  This method is called in any select method to check inputs.

        Parameters
        ----------
        **kwargs : dict or key=value
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
        ku = [k.upper() for k in kwargs.keys()]
        for k in ku:
            if k not in self:
                unrecognized.append(k)
        # print("KU, K", ku, k)
        if len(unrecognized) > 0:
            raise KeyError(f"The following keywords were not recognized: {unrecognized}")

    # @todo this could be made generic to check for any type
    def _check_numbers(self, **kwargs):
        # @todo allow Quantities
        """
        Check that a list of keyword arguments is all numbers.

        Parameters
        ----------
        **kwargs : dict or key=value
           Keyword arguments

        Raises
        ------
        ValueError
            If one or more of the values is not numeric.

        Returns
        -------
        None.

        """
        ku = np.ma.masked_array([k.upper() for k in kwargs.keys()])
        ku.mask = np.array([isinstance(x, numbers.Number) for x in kwargs.values()])
        if not np.all(ku.mask):
            raise ValueError(f"Expected numeric value for these keywords but did not get a number {ku[~ku.mask]}")

    def _check_for_duplicates(self, df):
        """
        Check that the user hasn't already added a rule matching this one

        Parameters
        ----------
        df : ~pandas.DataFrame
            The selection to check

        Raises
        ------
        Exception
            If an identical rule (DataFrame) has already been added.

        Returns
        -------
        None.

        """

        for s in self._selection_rules.values():
            if s.equals(df):
                raise Exception("A identical selection rule has already been added.")

    def _addrow(self, row, dataframe, tag=None):
        """
        Common code to add a tagged row to the internal table after the selection has been created.
        Should be called in select* methods.

        Parameters
        ----------
        tag : str, optional
            An identifying tag by which the rule may be referred to later.
            If None, a  randomly generated tag will be created.
        row : dict
            key, value pairs of the selection
        dataframe : ~pandas.DataFrame
            The dataframe created by the selection.

        Returns
        -------
        None.

        """
        self._check_for_duplicates(dataframe)
        if tag is not None:
            row["TAG"] = tag
        else:
            row["TAG"] = self._generate_tag(list(row.values()))
        row["ID"] = self._next_id
        self._selection_rules[row["ID"]] = dataframe
        self._table.add_row(row)

    def select(self, tag=None, **kwargs):
        """Add one or more exact selection rules, e.g., `key1 = value1, key2 = value2, ...`

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
        self._check_keys(**kwargs)
        row = dict()
        df = self
        for k, v in list(kwargs.items()):
            ku = k.upper()
            v = self._sanitize_input(ku, v)
            row[ku] = str(v)
            df = pd.merge(df, df[df[ku] == v], how="inner")
        self._addrow(row, df, tag)
        # return df

    def select_range(self, tag=None, **kwargs):
        """
        Select a range of inclusive values for a given key(s).
        e.g., `key1 = (v1,v2), key2 = (v3,v4), ...`
        Will select data  `v1 <= data1 <= v2, v3 <= data2 <= v4, ... `
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
        self._check_keys(**kwargs)
        row = dict()
        df = self
        # @todo need to handle upper and lower limits (None,v), (v,None), (v,)
        # @todo need to handle nested arrays [[v1,v2],[v3,v4]] - what's the use case?
        for k, v in list(kwargs.items()):
            ku = k.upper()
            row[ku] = str(v)
            if len(v) == 2:
                if v[0] is not None and v[1] is not None:
                    df = pd.merge(df, df[(df[ku] <= v[1]) & (df[ku] >= v[0])], how="inner")
                elif v[0] is None:  # upper limit given
                    df = pd.merge(df, df[(df[ku] <= v[1])], how="inner")
                else:  # lower limit given (v[1] is None)
                    df = pd.merge(df, df[(df[ku] >= v[0])], how="inner")
            elif len(v) == 1:  # lower limit given
                df = pd.merge(df, df[(df[ku] <= v[1])], how="inner")
            else:
                raise Exception(f"Couldn't parse value tuple {v} for key {k}")

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
        kw = dict()
        for k, v in kwargs.items():
            v1 = v[0] - v[1]
            v2 = v[0] + v[1]
            kw[k] = (v1, v2)
        self.select_range(tag, **kw)

    def remove(self, id=None, tag=None):
        """Remove (delete) a selection rule(s).
        You must specify either `id` or `tag` but not both. If there are
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
            matching = self._table[matching_indices]
            self._table.remove_rows(matching_indices)

            for i in matching["ID"]:
                del self._selection_rules[i]
                # self._selection_rules.pop(i, None) # also works

    def show(self):
        self._set_pprint_exclude_names()
        print(self._table)
