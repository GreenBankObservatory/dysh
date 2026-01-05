#!/usr/bin/env python3
import re
from pathlib import Path
from typing import Literal

import astropy.units as u
from astropy.table import Table
from astropy.units.quantity import Quantity
from astroquery.splatalogue import Splatalogue

from ..util import (
    append_docstr_nosections,
    docstring_parameter,
    get_project_data,
    minimum_list_match,
    minimum_string_match,
    replace_col_astype,
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
    "gbtlines",  # note gbtlines does not include recombination lines
    "gbtrecomb",
    # we don't have any of the below until splatalogue is fixed or Tony R. provides.
    #    "top20",
    #    "planetaryatmosphere",
    #    "hotcores",
    #    "darkclouds",
    #    "diffuseclouds",
    #    "comets",
    #    "agb/ppn/pn",
    #    "extragalactic",
]

_all_cats = _allowable_remote_cats + _allowable_local_cats
# Grab splatalogue keywords from its query method so we are always in sync with it.
# Remove description and first two frequency parameters since we have our own.
__splatdoc__ = Splatalogue.query_lines.__doc__


def all_cats():
    # needed to access dunder variable from outside this module
    return _all_cats


# Remove the first part of splatalogue doc in favor of our own
i = __splatdoc__.index("chemical_name")
__splatdoc__ = __splatdoc__.replace(__splatdoc__[0:i], "")
# Replace how astroquery displays the Return section with our typical format.
ours = "Returns\n-------\nTable\n\tAn `astropy.table.Table` containing the results of the search."
i = __splatdoc__.index("Returns")
__splatdoc__ = __splatdoc__.replace(__splatdoc__[i:], ours)


class SpectralLineSearchClass:
    def __init__(self):
        self._tables = {}
        self._create_recomb_lines()

    def _process_cat(self, cat: str | Path) -> str:
        """The input catalog can be a string or a Path.  If it is a string, first check
        if is one of the special strings indicating a local dysh catalog. Otherwise check
        that it is a valid path to a file.

        Returns
        -------
            string to a valid path or a special strings
        """
        if isinstance(cat, Path):
            return str(cat)
        if (mc := minimum_string_match(cat.lower(), _all_cats)) is not None:
            return mc
        cp = Path(cat)
        if cp.is_file():
            return str(cp)
        else:
            raise ValueError(f"Unrecognized catalog {cat}. Valid catalogs are {_all_cats} or a valid path name.")

    def _patch_line_lists(self, line_lists: list) -> list:
        # This is to fix an inconsistency in splatalogue data.
        # Although the keyword to trigger recombination line search is 'Recombination', the
        # returned value in the Table is 'Recomb'.  So you can't search the database and the
        # result table with the same keyword! Ugh!
        # We are guaranteed by minimum_list_match that line_lists will contain the full word
        # if it refers to recombination.
        lowerlist = list(map(str.lower, line_lists))
        if "recombination" in lowerlist:
            lowerlist[lowerlist.index("recombination")] = "Recomb"
        return lowerlist

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
        cache: bool = True,
        format: str = "ascii.ecsv",
        only_NRAO_recommended=True,
        **kwargs,
    ) -> Table:  # @todo should we return pandas DataFrame instead?
        """Query locally or remotely for lines and return a table object. The query returns lines
        with rest frequencies in the range [`min_frequency`,`max_frequency`].

        **Note:** If the search parameters result in no matches, a zero-length Table will be returned.

        Parameters
        ----------
        min_frequency : `~astropy.units.quantity.Quantity`
            The minimum frequency to search (or any :meth:`u.spectral` equivalent)

        max_frequency : `~astropy.units.quantity.Quantity`
            The maximum frequency to search (or any :meth:`u.spectral` equivalent)

        cat : str or Path
            The catalog to use.  One of: {0}  (minimum string match) or a valid Path to a local astropy-compatible table.  The local table
            must have all the columns listed in the `columns` parameter.

                - `'gbtlines'` is a local catalog of spectral lines between 300 MHz and 120 GHz with CDMS/JP log(intensity) > -9.

                - `'gbtrecomb'` is a local catalog of H, He, and C recombination lnes between 300 MHz and 120 GHz.
        columns: str or list
            The query result columns to include in the returned table.  Any of {1}. The default is all columns.
        cache: bool
            For a local file query, make an in-memory copy of the input catalog to be used in subsequent queries to this catalog.
        asynchronous: bool
            Use asynchronous query
        format: str
            Stringe describing the format of a local input table. Must be a valid `astropy.io.ascii <https://docs.astropy.org/en/latest/io/ascii/index.html>`_ format string.  Default is 'ascii.ecsv'
        """
        mc = self._process_cat(cat)
        # we are overlaying this kwarg with a parameter to expose that we are changing the default.
        kwargs.update({"only_NRAO_recommended": only_NRAO_recommended})
        # user-friendly keywords
        if kwargs.get("line_lists", None) is not None:
            kwargs["line_lists"] = minimum_list_match(kwargs["line_lists"], Splatalogue.ALL_LINE_LISTS, casefold=True)
        if kwargs.get("line_strengths", None) is not None:
            kwargs["line_strengths"] = minimum_list_match(
                kwargs["line_strengths"], Splatalogue.VALID_LINE_STRENGTHS, casefold=True
            )
        if kwargs.get("intensity_type", None) is not None:
            kwargs["intensity_type"] = minimum_string_match(
                kwargs["intensity_type"], Splatalogue.VALID_INTENSITY_TYPES, casefold=True
            )
        if mc == "splatalogue":
            if asynchronous:
                table = Splatalogue._parse_result(
                    Splatalogue.query_lines_async(
                        min_frequency,
                        max_frequency,
                        **kwargs,
                    )
                )
            else:
                table = Splatalogue.query_lines(min_frequency, max_frequency, **kwargs)
        else:
            # search a local table
            return self.localquery(min_frequency, max_frequency, cat=mc, columns=columns, cache=cache, **kwargs)
        replace_col_astype(table, "intintensity", float, -1e20)
        if columns is not None and len(table) != 0:
            return table[columns]
        else:
            return table

    @docstring_parameter(str(_all_cats), str(_default_columns_to_return))
    @u.quantity_input(min_frequency=u.GHz, equivalencies=u.spectral())
    @u.quantity_input(max_frequency=u.GHz, equivalencies=u.spectral())
    def localquery(
        self,
        min_frequency: Quantity,
        max_frequency: Quantity,
        cat: Literal[(x for x in _allowable_local_cats)] | Path = "gbtlines",
        columns: str | list = _default_columns_to_return,
        chemical_name: str | None = None,
        chem_re_flags: int = re.I,
        energy_min: float | None = None,
        energy_max: float | None = None,
        energy_type: Literal[(x for x in Splatalogue.VALID_ENERGY_TYPES)] | None = None,
        intensity_lower_limit=None,
        intensity_type: Literal[(x for x in Splatalogue.VALID_INTENSITY_TYPES)] | None = None,
        line_lists: Literal[(x for x in Splatalogue.ALL_LINE_LISTS)] | None = None,
        cache: bool = False,
        format: str = "ascii.ecsv",
        **kwargs,  # ignore the rest!
    ) -> Table:
        """Query a local file for lines and return a table object. The query returns lines
        with rest frequencies in the range [`min_frequency`,`max_frequency`].

        **Note:**
         - If the search parameters result no matches, a zero-length Table will be returned.
         - Many of the keywords are only supported if `cat='splatalogue'` because local tables do not have all the columns that the Splatalogue database has.

        Parameters
        ----------
        min_frequency : `~astropy.units.quantity.Quantity`
            The minimum frequency to search (or any :meth:`u.spectral` equivalent). No default.
        max_frequency : `~astropy.units.quantity.Quantity`
            The maximum frequency to search (or any :meth:`u.spectral` equivalent). No default.
        cat : str
            The catalog to use.  One of: {0}  (minimum string match) or a valid local astropy-compatible table.  The local table
            must have all the columns listed in the `columns` parameter.  The default is a GBT-specific line catalog, 'gbtlines'.
        columns: str or list
            The query result columns to include in the returned table.  Any of {1}. The default is all columns.
        chemical_name : str
            Name of the chemical to search for. Treated as a regular
            expression.  An empty set will match *any*
            species. Examples:

            ``'H2CO'`` - 13 species have H2CO somewhere in their formula.

            ``'Formaldehyde'`` - There are 8 isotopologues of Formaldehyde
                                 (e.g., H213CO).

            ``'formaldehyde'`` - Thioformaldehyde,Cyanoformaldehyde.

            ``'formaldehyde',chem_re_flags=re.I`` - Formaldehyde,thioformaldehyde,
                                                    and Cyanoformaldehyde.

            ``' H2CO '`` - Just 1 species, H2CO. The spaces prevent including
                           others.
        chem_re_flags : int
            See the `~re` module
        energy_min : `None` or float
            Energy range to include.  See `energy_type`
        energy_max : `None` or float
            Energy range to include.  See `energy_type`
        energy_type : ``'el_cm1'``, ``'eu_cm1'``, ``'eu_k'``, ``'el_k'``
            Type of energy to restrict.  L/U for lower/upper state energy,
            cm/K for *inverse* cm, i.e. wavenumber, or K for Kelvin
        intensity_lower_limit : `None` or float
            Lower limit on the intensity.  See `intensity_type`
        intensity_type : `None`, ``'CDMS/JPL (log)'``, ``'Sij-mu2'``, ``'Aij (log)'``
            The type of intensity on which to place a lower limit
        line_lists : list
            Options:
            Lovas, SLAIM, JPL, CDMS, ToyaMA, OSU, Recombination, RFI
        cache: bool
            Make an in-memory copy of the input table to be used in subsequent queries to this catalog.
        format: str
            The astropy IO format string for the input table.  Default is ECSV format.

        Returns
        -------
        Table
            An astropy table containing the results of the search
        """
        _cat = self._process_cat(cat)
        if _cat in _allowable_local_cats:
            _cat = str(get_project_data() / (_cat + ".csv.gz"))
        # using the cache is about 10x faster for gbtlines
        if _cat in self._tables:
            _table = self._tables[_cat]
        else:
            _table = Table.read(_cat, format=format)
        if cache:
            self._cache_local_table(_cat, _table)
        # now do the work of downselecting the table.
        # The easiest way to do this through pandas; using the Table interface
        # is too cumbersome.
        if chemical_name is not None:
            species = self.get_species_ids(species_regex=chemical_name, reflags=chem_re_flags)
            if len(species) == 0:
                raise ValueError(f"Unable to find species matching {chemical_name}")
            # get species id returns string but 'species_id' column in tables returned by splatalogue is int!
            splist = list(map(int, species.values()))
        df = _table.to_pandas()

        # Select the frequency range
        # fmt: off
        df = df[
            (
                (df["orderedfreq"] >= min_frequency.to("MHz").value) &
                (df["orderedfreq"] <= max_frequency.to("MHz").value)
            )
        ]
        # fmt: on

        # chemical name and re_flags via species_id
        if chemical_name is not None:
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
        if line_lists is not None:
            if (line_lists := minimum_list_match(line_lists, Splatalogue.ALL_LINE_LISTS, casefold=True)) is None:
                raise ValueError(
                    f"list_lists must be one or more of {Splatalogue.ALL_LINE_LISTS} (case insensitive, minimum match)."
                )
            line_lists = self._patch_line_lists(line_lists)
            df = df[df["linelist"].isin(line_lists)]
        # line strengths
        if intensity_lower_limit is not None:
            if intensity_type is None:
                raise ValueError(
                    f"If you specify an intensity lower limit, you must also specify its intensity_type. One of  {Splatalogue.VALID_INTENSITY_TYPES} (case insensitive, minimum_match)."
                )
            elif (
                intensity_type := minimum_string_match(intensity_type, Splatalogue.VALID_INTENSITY_TYPES, casefold=True)
            ) is None:
                raise ValueError(
                    f"intensity_type must be one of {Splatalogue.VALID_INTENSITY_TYPES} (case insensitive, minimum match ."
                )
            else:
                df = df[df["intintensity"] >= intensity_lower_limit]
        table = Table.from_pandas(df)
        # @todo Should we add units to the table?
        if columns is not None:
            return table[columns]
        else:
            return table

    @docstring_parameter(Splatalogue.get_species_ids.__doc__)
    def get_species_ids(self, species_regex, reflags=0, recache=False):
        """
        Convenience call-through to :meth:`~astroquery.splatalogue.SplatalogueClass.get_species_id`.
        {0}
        """
        return Splatalogue.get_species_ids(species_regex, reflags=reflags, recache=recache)

    def clear_cache(self):
        """
        Clear the local caches. This will clear the Splatalogue cache and any local tables
        that have been cached.

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
            You can choose a subset of these when performing a search.
        """
        return _default_columns_to_return

    def _create_recomb_lines(self):
        """
        Make the recombination line ascii to unicode map.  Also add some convenient aliases for users.
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

        # aliases since Splatalogue is picky about Case.
        self._altrecomb = {
            "H": "Hydrogen",
            "C": "Carbon",
            # "He" already works
            "hydrogen": "Hydrogen",
            "carbon": "Carbon",
            "helium": "Helium",
        }

    def recomb(
        self,
        min_frequency: Quantity,
        max_frequency: Quantity,
        line: str,  # @todo let this be a list? simillar to 'recomball"
        # cat: Literal[*_all_cats] | Path = "splatalogue",  # allowed in Python 3.11+
        cat: Literal[(x for x in _all_cats)] | Path = "splatalogue",
        columns: str | list = _default_columns_to_return,
        convert_to_unicode: bool = True,
        only_NRAO_recommended: bool = True,
        **kwargs,
    ) -> Table:
        """
        Search for recombination lines of H, He, and C.

        Parameters
        ----------
        min_frequency : `~astropy.units.quantity.Quantity`
            The minimum frequency to search (or any :meth:`u.spectral` equivalent). No default.
        max_frequency : `~astropy.units.quantity.Quantity`
            The maximum frequency to search (or any :meth:`u.spectral` equivalent). No default.
        line : str
            A string describing the line or series to search for, e.g. "Hydrogen", "Halpha", "Hebeta", "C", "carbon". No default.
        cat : str or Path
            The catalog to use.  One of: {0}  (minimum string match) or a valid Path to a local astropy-compatible table.  The local table
            must have all the columns listed in the `columns` parameter. Default is 'splatalogue'.
        columns: str or list
            The query result columns to include in the returned table.  Any of {1}. The default is all columns.
        convert_to_unicode : bool, optional
            Splatalogue stores line names using the unicode characters for Greek symbols, e.g. `\u03b1` for 'alpha'.  dysh will convert for you, if you put in e.g., 'Halpha'.
            You should only change this if a) you are inputing unicode or b) you are searching a local file that you know doesn't use unicode. The default is True.
        \\*\\*kwargs : dict
            Other keyword arguments supported by :meth:`query_lines` if `cat` is 'splatalogue'.

        Returns
        -------
        Table
            An astropy table containing the results of the search

        """
        if line in self._altrecomb:
            line = self._altrecomb[line]
        elif line is not None:
            for k in self._recomb_dict:
                if k in line:
                    line = line.replace(k, self._recomb_dict[k])
        return self.query_lines(
            min_frequency,
            max_frequency,
            chemical_name=line,
            cat=cat,
            line_lists=["Recomb"],
            columns=columns,
            only_NRAO_recommended=only_NRAO_recommended,
            chem_re_flags=re.I,
            **kwargs,
        )

    def recomball(
        self,
        min_frequency: Quantity,
        max_frequency: Quantity,
        cat: Literal[(x for x in _all_cats)] | Path = "splatalogue",
        columns: str | list = _default_columns_to_return,
        cache: bool = False,
        only_NRAO_recommended: bool = True,
        **kwargs,
    ) -> Table:
        """
        Fetch all recombination lines of H, He, C in the given frequency range from the catalog.

        Parameters
        ----------
        min_frequency : `~astropy.units.quantity.Quantity`
            The minimum frequency to search (or any :meth:`u.spectral` equivalent). No default.
        max_frequency : `~astropy.units.quantity.Quantity`
            The maximum frequency to search (or any :meth:`u.spectral` equivalent). No default.
        cat : str or Path
            The catalog to use.  One of: {0}  (minimum string match) or a valid Path to a local astropy-compatible table.  The local table
            must have all the columns listed in the `columns` parameter. Default is 'splatalogue'.
        columns: str or list
            The query result columns to include in the returned table.  Any of {1}. The default is all columns.
        cache: bool
            For a local file query, make an in-memory copy of the input table to be used in subsequent queries to this catalog.
        \\*\\*kwargs : dict
            Other keyword arguments supported by :meth:`query_lines` if `cat` is 'splatalogue'.

        Returns
        -------
        Table
            An astropy table containing the results of the search.

        """
        return self.recomb(
            min_frequency=min_frequency,
            max_frequency=max_frequency,
            cat=cat,
            line=None,
            cache=cache,
            columns=columns,
            only_NRAO_recommended=only_NRAO_recommended,
            **kwargs,
        )
