#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  6 10:49:11 2025

@author: mpound
"""

from pathlib import Path
from typing import Literal, Optional

import astropy.units as u
from astropy.quantity import Quantity
from astropy.table import Table
from astroquery.splatalogue import Splatalogue

from dysh.util import get_project_data, minimum_string_match

from . import docstring_parameter, prepend_docstr_nosections

_VALID_EXCLUDE = ("potential", "atmospheric", "probable", "known", "none")

_default_columns_to_return = [
    [
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
        "orderedFreq",
        "measFreq",
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
        "transitionBandColor",
        "searchErrorMessage",
    ]  #'sqlquery',
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


class SpectralLineSearch:
    def __init__(self):
        self._tables = {}
        self._create_recomb_lines()

    @docstring_parameter(_all_cats, _default_columns_to_return)
    @prepend_docstr_nosections("\n" + Splatalogue._parse_kwargs.__doc_)
    @u.quantity_input(min_frequency=u.GHz, equivalencies=u.spectral())
    @u.quantity_input(max_frequency=u.GHz, equivalencies=u.spectral())
    def query_lines(
        self,
        min_frequency: Quantity,
        max_frequency: Quantity,
        # min_energy : Quantity = None,
        # max_energy :Quantity = None,
        cat: Literal(_all_cats) | Path = "splatalogue",
        columns: str | list = _default_columns_to_return,
        asynchronous: bool = False,
        **kwargs,
    ) -> Table:  # @todo should we return pandas DataFrame instead?
        """Query the locally or remotely for lines and return a table object.

        Parameters
        ----------
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
        kwargs["line_lists"] = list(minimum_string_match(kwargs["line_lists"]))
        kwargs["line_strengths"] = list(minimum_string_match(kwargs["line_strengths"]))
        if mc == "splatalogue":
            if asynchronous:
                table = Splatalogue._parse_result(Splatalogue.query_lines_async(min_frequency, max_frequency, **kwargs))
            else:
                table = Splatalogue.query(min_frequency, max_frequency, **kwargs)
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
        cat: str | Path = "top20lines",
        columns: str | list = _default_columns_to_return,
        chemical_name: Optional[str] = None,
        chem_re_flags: int = 0,
        energy_min: Optional[float] = None,
        energy_max: Optional[float] = None,
        energy_type: Literal(Splatalogue.VALID_ENERGY_TYPES) = None,
        intensity_lower_limit=None,
        intensity_type: Literal(Splatalogue.VALID_INTENSITY_TYPES) = None,
        exclude: Literal(_VALID_EXCLUDE) = None,
        line_lists: Literal(Splatalogue.ALL_LINE_LISTS) = None,
        line_strengths: Literal(Splatalogue.VALID_LINE_STRENGTHS) = None,
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
        species = Splatalogue.get_species_ids(chemical_name)
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

        table = Table.from_pandas(df)
        if columns is not None:
            return table[columns]
        else:
            return table

    def clear_cache(self):
        Splatalogue.clear_cache()
        self._tables = {}

    def _cache_local_table(self, tablename: str, table: Table) -> None:
        # is this a bad idea? memory hog?
        self._tables[tablename] = table

    @property
    def colnames(self):
        return self._tier1table.colnames

    def _create_recomb_lines(self, lines):
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

    def recomb(
        self,
        min_frequency: Quantity,
        max_frequency: Quantity,
        line: str,
        cat: str | Path = "splatalogue",
    ) -> Table:
        pass
