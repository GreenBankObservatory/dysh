#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from pathlib import Path
from typing import Literal, Optional

import astropy.units as u
from astropy.table import Table
from astropy.units.quantity import Quantity
from astroquery.splatalogue import Splatalogue

from ..util import (
    append_docstr_nosections,
    docstring_parameter,
    get_project_data,
    minimum_string_match,
)

_VALID_EXCLUDE = ("potential", "atmospheric", "probable", "known", "none")

# @todo this could be a dict with keys=colname and values=description,
# then passed to fancy docstring manipulation.
_default_columns_to_return = [
    "species_id",
    "name",
    "chemical_name",
    "resolved_QNs",
    "linelist",
    "LovasASTIntensity",
    "lower_state_energy",
    "upper_state_energy",
    "sijmu2",
    "sij",
    "aij",
    "intintensity",
    "Lovas_NRAO",
    "orderedfreq",
    "lower_state_energy_K",
    "upper_state_energy_K",
    #    "orderedFreq",
    #    "measFreq",
    "upperStateDegen",
    "moleculeTag",
    "qnCode",
    "labref_Lovas_NIST",
    "rel_int_HFS_Lovas",
    "unres_quantum_numbers",
    "lineid",
    "transition_in_space",
    "transition_in_G358",
    "obsref_Lovas_NIST",
    "source_Lovas_NIST",
    "telescope_Lovas_NIST",
    #    "transitionBandColor",
    "searchErrorMessage",
    #'sqlquery',
    #'requestnumber']
]

_allowable_remote_cats = ["splatalogue"]
_allowable_local_cats = [
    "top20",
    "planetaryatmosphere",
    "hotcores",
    "darkclouds",
    "diffuseclouds",
    "comets",
    "agb/ppn/pn",
    "extragalactic",
]

_all_cats = _allowable_remote_cats + _allowable_local_cats
# Grab splatalogue keywords from its query method so we are always in sync with it.
# Remove description and first two frequency parameters since we have our own.
__splatdoc__ = Splatalogue.query_lines.__doc__
i = __splatdoc__.index("chemical_name")
__splatdoc__ = __splatdoc__.replace(__splatdoc__[0:i], "")


class SpectralLineSearchClass:
    def __init__(self):
        self._tables = {}
        self._create_recomb_lines()

    # order of these decorators matters!

    @append_docstr_nosections(__splatdoc__, sections=[])
    @docstring_parameter(str(_all_cats), str(_default_columns_to_return))
    @u.quantity_input(min_frequency=u.GHz, equivalencies=u.spectral())
    @u.quantity_input(max_frequency=u.GHz, equivalencies=u.spectral())
    def query_lines(
        self,
        min_frequency: Quantity,
        max_frequency: Quantity,
        # cat: Literal[*_all_cats] | Path = "splatalogue",  # * in index allowed in Python 3.11+
        cat: Literal[(x for x in _all_cats)] | Path = "splatalogue",
        columns: str | list = _default_columns_to_return,
        asynchronous: bool = False,
        **kwargs,
    ) -> Table:  # @todo should we return pandas DataFrame instead?
        """Query the locally or remotely for lines and return a table object. The query returns lines
        with rest frequencies in the range [`min_frequency`,`max_frequency`].

        Parameters
        ----------
        min_frequency : `~astropy.units.quantity.Quantity`
            The minimum frequency to search (or any :meth:`u.spectral` equivalent)

        max_frequency : `~astropy.units.quantity.Quantity`
            The maximum frequency to search (or any :meth:`u.spectral` equivalent)

        cat : str or Path
            The catalog to use.  One of: {0}  (minimum string match) or a valid Path to a local astropy-compatible table.  The local table
            must have all the columns listed in the `columns` parameter.

        columns: str or list
            The query result columns to include in the returned table.  Any of {1}. The default is all columns.

        asynchronous: bool
            Use asynchronous query

        cache: bool
            Search the latest cache instead of doing a new remote query.
            See https://astroquery.readthedocs.io/en/latest/index.html#caching

        """
        if isinstance(cat, Path):
            mc = cat
        else:
            mc = minimum_string_match(cat.lower(), _all_cats)
        if mc is None:
            raise ValueError(f"Unrecognized catalog {cat}. Valid catalogs are {_all_cats}.")
        # user-friendly keywords
        if False:
            if "line_lists" in kwargs:
                kwargs["line_lists"] = list(minimum_string_match(kwargs["line_lists"], Splatalogue.ALL_LINE_LISTS))
            if "line_strengths" in kwargs:
                kwargs["line_strengths"] = list(
                    minimum_string_match(kwargs["line_strengths"], Splatalogue.VALID_LINE_STRENGTHS)
                )
        if mc == "splatalogue":
            if asynchronous:
                table = Splatalogue._parse_result(Splatalogue.query_lines_async(min_frequency, max_frequency, **kwargs))
            else:
                table = Splatalogue.query_lines(min_frequency, max_frequency, **kwargs)
        else:
            # search a local table
            return self.localquery(min_frequency, max_frequency, cat=mc, columns=columns, **kwargs)

        if columns is not None:
            return table[columns]
        else:
            return table

    @u.quantity_input(min_frequency=u.GHz, equivalencies=u.spectral())
    @u.quantity_input(max_frequency=u.GHz, equivalencies=u.spectral())
    def localquery(
        self,
        min_frequency: Quantity,
        max_frequency: Quantity,
        cat: Literal[(x for x in _allowable_local_cats)] | Path = "top20",
        columns: str | list = _default_columns_to_return,
        chemical_name: Optional[str] = None,
        chem_re_flags: int = 0,
        energy_min: Optional[float] = None,
        energy_max: Optional[float] = None,
        energy_type: Literal[(x for x in Splatalogue.VALID_ENERGY_TYPES)] | None = None,
        intensity_lower_limit=None,
        intensity_type: Literal[(x for x in Splatalogue.VALID_INTENSITY_TYPES)] | None = None,
        line_lists: Literal[(x for x in Splatalogue.ALL_LINE_LISTS)] | None = None,
        line_strengths: Literal[(x for x in Splatalogue.VALID_LINE_STRENGTHS)] | None = None,
        cache: bool = False,
        **kwargs,  # ignore the rest!
    ) -> Table:
        if isinstance(cat, str) and cat not in _allowable_local_cats:
            raise ValueError(f"Unrecognized catalog {cat}. Must be one of {_allowable_local_cats}")
        if not isinstance(cat, Path):
            _cat = get_project_data / (cat + ".tab.gz")
        else:
            _cat = cat
        if _cat in self._tables:
            _table = self._tables[cat]
        else:
            _table = Table.read(_cat)
        if cache:
            self._cache_local_table(str(_cat), _table)
        # now do the work of downselecting the table.
        # The easiest way to do this through pandas; using the Table interface
        # is too cumbersome.
        species = self.get_species_ids(chemical_name)
        splist = list(species.values())
        df = _table.to_pandas()

        # Select the frequency range
        # fmt: off
        df = df[
            (
                df["orderedfreq"] >= min_frequency.to("MHz").value &
                df["orderedfreq"] <= max_frequency.to("MHz").value
            )
        ]
        # fmt: on

        # chemical name and re_flags via species_id
        df = df[df["species_id"].isin(splist)]

        # energies
        if energy_type == "el_cm1":
            k = "lower_state_energy"
            df = df[(df[k] >= energy_min & df[k] <= energy_max)]
        elif energy_type == "eu_cm1":
            k = "upper_state_energy"
            df = df[(df[k] >= energy_min & df[k] <= energy_max)]
        elif energy_type == "el_k":
            k = "lower_state_energy_k"  # helps with auto-formatting of next line
            df = df[(df[k] >= energy_min & df[k] <= energy_max)]
        elif energy_type == "eu_k":
            k = "upper_state_energy_k"
            df = df[(df[k] >= energy_min & df[k] <= energy_max)]

        # line lists
        df = df[df["linelist"].isin(list(line_lists))]
        # line strengths
        df = df[df["intensity_type"] >= intensity_lower_limit]

        table = Table.from_pandas(df)
        # @todo Should we add units to the table?
        if columns is not None:
            return table[columns]
        else:
            return table

    @docstring_parameter(Splatalogue.get_species_ids.__doc__)
    def get_species_ids(self, species_regex=None, *, reflags=0, recache=False):
        """
        Convenience call-through to :meth:`~astroquery.splatalogue.SplatalogueClass.get_species_id`.
        {0}
        """
        return Splatalogue.get_species_ids(species_regex, reflags, recache)

    def clear_cache(self):
        """
        Clear the local caches. This will clear the Splatalogue cache and any local tables
        that have been cached.Splatalogue.ALL_LINE_LISTS

        See https://astroquery.readthedocs.io/en/stable/splatalogue/splatalogue.html#troubleshooting

        Returns
        -------
        None.
        """
        Splatalogue.clear_cache()
        self._tables = {}

    def _cache_local_table(self, tablename: str, table: Table) -> None:
        """
        Cache local files in a dict so they don't have to be opened each time.

        Parameters
        ----------
        tablename : str
            table name key string
        table : Table
            the table to cache

        Returns
        -------
        None
        """
        # is this a bad idea? memory hog?
        self._tables[tablename] = table

    @property
    def colnames(self):
        """
        Returns
        -------
        colname : list
            The list of column names present in the default returned table.
            You can choose a subset of these when performing a :meth:`search`.
        """
        return _default_columns_to_return

    def _create_recomb_lines(self):
        """
        make the recombination line ascii to unicode map
        """
        # Splatalogue wants the unicode Greek characters in recombination lines.
        # Create a mapping to allow users to type in e.g. Halpha instead of H\u03B1 or Hα # noqa

        self._recomb_dict = {}
        unicode_map = {
            "alpha": "\u03b1",
            "beta": "\u03b2",
            "gamma": "\u03b3",
            "delta": "\u03b4",
            "epsilon": "\u03b6",
            "zeta": "\u0364",
        }
        for line in ["H", "He", "C"]:
            for k, v in unicode_map.items():
                self._recomb_dict[f"{line}{k}"] = f"{line}{v}"
        self._altrecomb = {
            "H": "Hydrogen",
            "C": "Carbon",
            "hydrogen": "Hydrogen",
            "carbon": "Carbon",
            "helium": "Helium",
        }

    def recomb(
        self,
        min_frequency: Quantity,
        max_frequency: Quantity,
        line: str,
        # cat: Literal[*_all_cats] | Path = "splatalogue",  # allowed in Python 3.11+
        cat: Literal[(x for x in _all_cats)] | Path = "splatalogue",
        convert_to_unicode: bool = True,
        **kwargs,
    ) -> Table:
        """
        Search for recombination lines of H, He, and C.

        Parameters
        ----------
        min_frequency : `~astropy.units.quantity.Quantity`
            The minimum frequency to search
        max_frequency : `~astropy.units.quantity.Quantity`
            The maximum frequency to search
        line : str
            A string describing the line or series to search for, e.g. "Hydrogen", "Halpha", "Hebeta", "C", "carbon"
        cat : str or Path
            The catalog to use.  One of: {0}  (minimum string match) or a valid Path to a local astropy-compatible table.  The local table
            must have all the columns listed in the `columns` parameter. Default is 'splatalogue'.
        convert_to_unicode : bool, optional
            Splatalogue stores line names using the unicode characters for Greek symbols, e.g. `\u03b1` for 'alpha'.  dysh will convert for you, if you put in e.g., 'Halpha'.
            You should only change this if a) you are inputing unicode or b) you are searching a local file that you know doesn't use unicode. The default is True.

        Returns
        -------
        Table
            An astropy table containing the results of the search

        """
        if line in self._altrecomb:
            line = self._altrecomb[line]
        else:
            for k in self._recomb_dict:
                if k in line:
                    line = line.replace(k, self._recomb_dict[k])

        return self.query_lines(
            min_frequency,
            max_frequency,
            chemical_name=line,
            cat=cat,
            line_lists=["Recombination"],
            **kwargs,
        )


SpectralLineSearch = SpectralLineSearchClass()
