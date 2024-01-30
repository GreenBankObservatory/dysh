import hashlib

import numpy as np
import pandas as pd
from astropy.table import Table


class Selection(object):
    def __init__(self, **kwargs):
        # get these from a recent SDFITS file somehow
        self._available_keys = []
        # add ID and TAG to available keys before creating tag
        updated = np.array(...)
        dt = np.array([str] * len(updated))
        # make a table with columns and str dtype
        self._table = Table(data=none, names=updated, dtype=dt)

    def _parse(self, key, value):
        pass

    def _parse_coordinates(self, value):
        pass

    def _parse_time(self, value):
        pass

    def _parse_value(self, value):
        """return a string that can be stored in a Table cell"""
        pass

    def _generate_tag(self, values, hashlen=9):
        data = "".join(map(str, values))
        hash_object = hashlib.sha256(data.encode())
        unique_id = hash_object.hexdigest()
        return unique_id[0:hashlen]

    def select(self, tag=None, **kwargs):
        """Add a selection rule

        Parameters
        ----------
            tag : str
                An identifying tag by which the rule may be referred to later.
            key : str
                The key value (SDFITS column name)
            value : any
                The value to match

        """
        row = dict(zip(self.available_keys, [None] * len(self.available_keys)))
        for k in list(kwargs.keys()):
            row[key] = self._parse_value(kwargs[k])
        if tag is not None:
            row["tag"] = tag
        else:
            row["tag"] = self._generate_tag(list(row.values()))
        self._table.add_row(list(row.values()))

    def remove(self, id=None, tag=None):
        """Remove (delete) a selection rule
        You must specifiy either `id` or `tag` but not both.
        Parameters
        ----------
            id : int
                The ID number of the rule as displayed in `show()`
            tag : str
                An identifying tag by which the rule may be referred to later.
        """
        if id is not None and tag is not None:
            raise Exception("You can only specify one of id or tag")
        # if id is None:
        #    #remove row(s) based on tag

    #
    #        else:
    #            table.

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

    @property
    def available_keys(self):
        pass

    def add_key(self, key):
        pass
