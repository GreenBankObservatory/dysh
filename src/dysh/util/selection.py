import numpy as np
import pandas as pd
from astropy.table import PprintIncludeExclude, Table, TableAttribute

from ..fits import default_sdfits_columns
from . import generate_tag

# add ID and TAG to available keys before creating tag
# why not move this to constructor?
idtag = ["ID", "TAG"]
DEFKEYS = np.array(default_sdfits_columns())
DEFKEYS = np.insert(DEFKEYS, 0, idtag)


class Selection(Table):
    foobar = TableAttribute()  # example of adding a custom attribute

    def __init__(self, *args, **kwargs):
        dt = np.array([str] * len(DEFKEYS))
        dt[0] = int  # I
        # make a table with columns and str dtype
        # self.foobar = "hello"
        # print(DEFKEYS)
        super().__init__(data=None, names=DEFKEYS, dtype=dt)
        for t in idtag:
            self.add_index(t)

    def _set_pprint_exclude_names(self):
        """Use pprint_exclude_names to set the list
        columns that have no entries.
        """
        emptycols = np.array(self.colnames)[
            [np.all([self[k].data[i] == "" for i in range(len(self))]) for k in self.colnames]
        ]
        self.pprint_exclude_names.set(emptycols)

    def __repr__(self):
        # when printing to screen we only want to include
        # columns that are not empty.
        # This will not affect writing to a file.
        self._set_pprint_exclude_names()
        return super().__repr__()

    def _parse(self, key, value):
        """


        Parameters
        ----------
        key : TYPE
            DESCRIPTION.
        value : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        pass

    def _parse_coordinates(self, value):
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

    def _parse_value(self, value):
        """
        return a string that can be stored in a Table cell

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
        return sorted(self["ID"])[-1] + 1

    # all SDFITS keywords are uppercase, but we can allow
    # for lower/mixed access this way.
    # This may not actually be helpful
    def __getitem__(self, item):
        if isinstance(item, str):
            return super().__getitem__(item.upper())
        else:
            return super().__getitem__(item)

    def select(self, tag=None, **kwargs):
        """Add a selection rule

        Parameters
        ----------
            tag : str
                An identifying tag by which the rule may be referred to later.
            key : str
                The key  (SDFITS column name)
            value : any
                The value or expression to select

        """
        row = dict()
        for k in list(kwargs.keys()):
            row[k] = self._parse_value(kwargs[k])
        if tag is not None:
            row["TAG"] = tag
        else:
            row["TAG"] = self._generate_tag(list(row.values()))
        row["ID"] = self._nextid
        self._table.add_row(list(row.values()))

    def remove(self, id=None, tag=None):
        """Remove (delete) a selection rule.
        You must specify either `id` or `tag` but not both.

        Parameters
        ----------
            id : int
                The ID number of the rule as displayed in `show()`
            tag : str
                An identifying tag by which the rule may be referred to later.
        """
        if id is not None and tag is not None:
            raise Exception("You can only specify one of id or tag")
        if id:
            try:
                row = self.loc["ID", id]
            except KeyError:
                raise KeyError(f"No ID = {id} found in this Selection")
        else:
            try:
                row = self.loc["TAG", tag]
            except KeyError:
                raise KeyError(f"No TAG = {tag} found in this Selection")
        self.remove_row(row.index)

    def show(self):
        pass

    def to_pandas(self):
        """
        Convert Selection to a pandas DataFrame
        Returns
        -------
            df : ~pandas.DataFrame
        """
        pass
