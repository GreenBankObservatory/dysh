"""Load SDFITS files produced by the Green Bank Telescope"""

import copy
import inspect
import itertools
import numbers
import os
import platform
import time
import warnings
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.io import fits
from astropy.units.quantity import Quantity

from dysh.log import logger

from ..coordinates import Observatory, decode_veldef, eq2hor, hor2eq
from ..log import HistoricalBase, log_call_to_history, log_call_to_result
from ..spectra.core import mean_data
from ..spectra.scan import (
    FSScan,
    NodScan,
    PSScan,
    ScanBase,
    ScanBlock,
    Spectrum,
    SubBeamNodScan,
    TPScan,
)
from ..spectra.tcal import TCal
from ..spectra.vane import VaneSpectrum
from ..util import (
    Flag,
    Selection,
    calc_vegas_spurs,
    consecutive,
    convert_array_to_mask,
    eliminate_flagged_rows,
    get_valid_channel_range,
    keycase,
    select_from,
    show_dataframe,
    uniq,
)
from ..util.calibrator import Calibrator
from ..util.files import dysh_data
from ..util.gaincorrection import GBTGainCorrection
from ..util.selection import Flag, Selection  # noqa: F811
from ..util.weatherforecast import GBTWeatherForecast
from . import conf, core
from .sdfitsload import SDFITSLoad

# from GBT IDL users guide Table 6.7
# @todo what about the Track/OnOffOn in e.g. AGBT15B_287_33.raw.vegas  (EDGE HI data)
# _PROCEDURES = ["Track", "OnOff", "OffOn", "OffOnSameHA", "Nod", "SubBeamNod"]


class GBTFITSLoad(SDFITSLoad, HistoricalBase):
    """
    GBT-specific container to represent one or more SDFITS files

    Parameters
    ----------
    fileobj : str or `pathlib.Path`
        File to read or directory path.  If a directory, all
        FITS files within will be read in.
    source  : str
        target source to select from input file(s). Default: all sources
    hdu : int or list
        Header Data Unit to select from input file. Default: all HDUs

    skipflags: bool
        If True, do not read any flag files associated with these data. Default:False

    flag_vegas: bool
        If True, flag VEGAS spurs using the algorithm described in :meth:`~dysh.util.core.calc_vegas_spurs`
        and ignore VEGAS_SPUR flag rules in flag files. Note this parameter is independent of 'skip_flags', which
        controls only the reading of the flag file.  If you want no flags at all, use `skipflags=True, flag_vegas=False`.


        +---------+-----------+--------------------------------------------------------------------------------------------+
        |skipflags|flag_vegas | behavior                                                                                   |
        +=========+===========+============================================================================================+
        |False    | False     | VEGAS and other flags are created based on the flags file                                  |
        +---------+-----------+--------------------------------------------------------------------------------------------+
        |True     | False     | No flags are created                                                                       |
        +---------+-----------+--------------------------------------------------------------------------------------------+
        |True     | True      | VEGAS flags are created based on the FITS header                                           |
        +---------+-----------+--------------------------------------------------------------------------------------------+
        |False    | True      | VEGAS flags are created based on the FITS header.  Other flags are read from the flags file|
        +---------+-----------+--------------------------------------------------------------------------------------------+
    """

    @log_call_to_history
    def __init__(self, fileobj, source=None, hdu=None, skipflags=False, flag_vegas=True, **kwargs):
        kwargs_opts = {
            "index": True,
            "verbose": False,
            "fix_ka": True,
        }  # only set index to False for performance testing.
        HistoricalBase.__init__(self)
        kwargs_opts.update(kwargs)
        path = Path(fileobj)
        self._sdf = []
        self._selection = None
        self._tpnocal = None  # should become True or False once known
        self._flag = None

        self.GBT = Observatory["GBT"]
        if path.is_file():
            logger.debug(f"Treating given path {path} as a file")
            self._sdf.append(SDFITSLoad(path, source, hdu, **kwargs_opts))
            if not hasattr(self, "_filename"):
                self._filename = self._sdf[0].filename
        elif path.is_dir():
            logger.debug(f"Treating given path {path} as a directory")
            # Find all the FITS files in the directory and sort alphabetically
            # because e.g., VEGAS does A,B,C,D,E
            nf = 0  # performance testing
            for f in sorted(path.glob("*.fits")):
                logger.debug(f"Selecting {f} to load")
                if kwargs.get("verbose", None):
                    print(f"Loading {f}")
                if nf < kwargs.get("nfiles", 99999):  # performance testing limit number of files loaded
                    self._sdf.append(SDFITSLoad(f, source, hdu, **kwargs_opts))
                    nf += 1
            if len(self._sdf) == 0:  # fixes issue 381
                raise Exception(f"No FITS files found in {fileobj}.")
            if not hasattr(self, "_filename"):
                self._filename = self._sdf[0].filename.parent
            self.add_history(f"This GBTFITSLoad encapsulates the files: {self.filenames()}", add_time=True)
        else:
            raise Exception(f"{fileobj} is not a file or directory path")
        # Add in any history/comment that were in the previous file(s)
        for sdf in self._sdf:
            for h in sdf._hdu:
                self.add_history(h.header.get("HISTORY", []))
                self.add_comment(h.header.get("COMMENT", []))
        self._remove_duplicates()
        if kwargs_opts["index"]:
            self._create_index_if_needed(skipflags, flag_vegas)
            self._update_radesys()
            # This only works if the index was created.
            if kwargs_opts["fix_ka"]:
                self._fix_ka_rx_if_needed()
        # We cannot use this to get mmHg as it will disable all default astropy units!
        # https://docs.astropy.org/en/stable/api/astropy.units.cds.enable.html#astropy.units.cds.enable
        # cds.enable()  # to get mmHg

        # ushow/udata depend on the index being present, so check that index is created.
        if kwargs.get("verbose", None) and kwargs_opts["index"]:
            print(f"==GBTLoad {fileobj}")
            self.ushow("OBJECT", 0)
            self.ushow("SCAN", 0)
            self.ushow("SAMPLER", 0)
            self.ushow("PLNUM")
            self.ushow("IFNUM")
            self.ushow("FDNUM")
            self.ushow("SIG", 0)
            self.ushow("CAL", 0)
            self.ushow("PROCSEQN", 0)
            self.ushow("PROCSIZE", 0)
            self.ushow("OBSMODE", 0)
            self.ushow("SIDEBAND", 0)

        lsdf = len(self._sdf)
        if lsdf > 1:
            print(f"Loaded {lsdf} FITS files")
        if kwargs_opts["index"]:
            self.add_history(f"Project ID: {self.projectID}", add_time=True)
        else:
            print("Reminder: No index created; many functions won't work.")

        self._qd_corrected = False

    def __repr__(self):
        return str(self.files)

    def __str__(self):
        return str(self.filenames)

    @property
    def _index(self):
        # for backwards compatibility after removing _index
        # as a separate object
        return self._selection

    @property
    def projectID(self):
        """
        The project identification

        Returns
        -------
        str
            The project ID string
        """
        return uniq(self["PROJID"])[0]

    @property
    def total_rows(self):
        """Returns the total number of rows summed over all files and binary table HDUs"""
        return sum([s.total_rows for s in self._sdf])

    @property
    def columns(self):
        """The column names in the binary table, minus the DATA column

        Returns
        -------
        `~pandas.Index`
            The column names as a DataFrame Index
        """
        # return a list instead?
        return self._selection.columns

    @property
    def selection(self):
        """
        The data selection object

        Returns
        -------
        `~dysh.util.Selection`
            The Selection object

        """
        return self._selection

    @property
    def final_selection(self):
        """
        The merged selection rules in the Selection object.
        See :meth:`~dysh.util.Selection.final`

        Returns
        -------
        `~pandas.DataFrame`
            The final merged selection

        """
        return self._selection.final

    @property
    def files(self):
        """
        The list of SDFITS file(s) that make up this GBTFITSLoad object

        Returns
        -------
        files : list
            list of `~PosixPath` objects

        """
        files = []
        for sdf in self._sdf:
            files.append(sdf.filename)
        return files

    @property
    def flags(self):
        """
        The data flag object

        Returns
        -------
        `~dysh.util.Flag`
            The Flag object

        """
        return self._flag

    @property
    def final_flags(self):
        # this method is not particularly useful. consider removing it
        """
        The merged flag rules in the Flag object.
        See :meth:`~dysh.util.SelectionBase.final`

        Returns
        -------
        `~pandas.DataFrame`
            The final merged flags

        """
        # all_channels_flagged = np.where(self._table["CHAN"] == "")j
        return self._flag.final

    def filenames(self):
        """
        The list of SDFITS filenames(s) that make up this GBTFITSLoad object

        Returns
        -------
        filenames : list
            list of str filenames

        """
        return [p.as_posix() for p in self.files]

    def index(self, hdu=None, bintable: int = None, fitsindex=None):  # noqa: RUF013
        """
        Return The index table

        Parameters
        ----------
        hdu : int or list
            Header Data Unit to select from the index. Default: all HDUs
        bintable :  int
            The index of the `bintable` attribute, None means all bintables
        fitsindex: int
            The index of the FITS file contained in this GBTFITSLoad.
            Default:None meaning return one index over all files.

        Returns
        -------
        index : `~pandas.DataFrame`
            The index of this GBTFITSLoad

        """
        if fitsindex is None:
            df = self._selection
        else:
            df = self._sdf[fitsindex]._index

        if hdu is None and bintable is None:
            return df
        if hdu is not None:
            df = df[df["HDU"] == hdu]
        if bintable is not None:
            df = df[df["BINTABLE"] == bintable]
        return df

    def stats(self, bintable=0):
        """
        Return some basic statistics of the GBTFITSLoad.
        Useful for performance testing.  A dictionary with
        the following keys and values is returned:

            nfiles : number of FITS files
            nrows  : number of data rows
            fdnum  : number of unique feeds
            ifnum  : number of unique IFs
            plnum  : number of unique polarizations
            sig    : number of unique SIG integrations
            cal    : number of unique CAL integrations

        Parameters
        ----------
        bintable :  int
            The index of the `bintable` attribute to probe.

        Returns
        -------
        stats : dict
            A dictionary with keys
        """

        s = {}
        df = self.index(bintable=bintable)
        s["nrows"] = len(df)
        s["nfiles"] = len(self.files)
        for k in ["fdnum", "ifnum", "plnum", "intnum", "sig", "cal"]:
            s[k] = len(uniq(df[k.upper()]))
        s["nchan"] = self._sdf[0].nchan(0)
        return s

    # override sdfits version
    def rawspectra(self, bintable, fitsindex, setmask=False):
        """
        Get the raw (unprocessed) spectra from the input bintable.

        Parameters
        ----------
        bintable :  int
            The index of the `bintable` attribute
        fitsindex: int
            the index of the FITS file contained in this GBTFITSLoad.  Default:0
        setmask : boolean
            If True, set the mask according to the current flags. Default:False

        Returns
        -------
        rawspectra : `~numpy.ndarray`
            The DATA column of the input bintable, masked according to `setmask`

        """
        return self._sdf[fitsindex].rawspectra(bintable, setmask=setmask)

    def rawspectrum(self, i, bintable=0, fitsindex=0, setmask=False):
        """
        Get a single raw (unprocessed) spectrum from the input bintable.

        Parameters
        ----------
        i :  int
            The row index to retrieve.
        bintable :  int or None
            The index of the `bintable` attribute. If None, the underlying bintable is computed from i
        fitsindex: int
            the index of the FITS file contained in this GBTFITSLoad.  Default:0
        setmask : bool
            If True, set the data mask according to the current flags. Default:False
            Note: if :meth:`apply_flags` has not been called, flags will not yet be set.
        Returns
        -------
        rawspectrum : `~numpy.ma.MaskedArray`
            The i-th row of DATA column of the input bintable, masked according to `setmask`

        """
        return self._sdf[fitsindex].rawspectrum(i, bintable, setmask=setmask)

    def getspec(self, i, bintable=0, observer_location=Observatory["GBT"], fitsindex=0, setmask=False):
        """
        Get a row (record) as a Spectrum

        Parameters
        ----------
        i : int
            The record (row) index to retrieve
        bintable : int, optional
             The index of the `bintable` attribute. default is 0.
        observer_location : `~astropy.coordinates.EarthLocation`
            Location of the observatory. See `~dysh.coordinates.Observatory`.
            This will be transformed to `~astropy.coordinates.ITRS` using the time of
            observation DATE-OBS or MJD-OBS in
            the SDFITS header.  The default is the location of the GBT.
        fitsindex: int
            the index of the FITS file contained in this GBTFITSLoad.  Default:0
        setmask : bool
            If True, set the data mask according to the current flags. Default:False
            Note: if :meth:`apply_flags` has not been called, flags will not yet be set.
        Returns
        -------
        s : `~dysh.spectra.spectrum.Spectrum`
            The Spectrum object representing the data row.

        """
        return self._sdf[fitsindex].getspec(i, bintable, observer_location, setmask=setmask)

    def _validate_summary_columns(self, columns, col_defs, needed=None, verbose=False):
        """
        Sanitize `columns` for `~dysh.fits.gbtfitsload.GBTFITSLoad.get_summary`.

        Parameters
        ----------
        columns : list or str
            Columns to sanitize. If a string, multiple column names must be comma separated.
        col_defs : dict
            Dictionary with column definitions. See `~dysh.fits.core.summary_column_definitions` for the expected format.
        needed : list
            List of columns needed to build the index.
        verbose : bool
            `~dysh.fits.gbtfitsload.GBTFITSLoad.get_summary` verbose mode.
        """
        if isinstance(columns, str):
            # Remove spaces and split by commas.
            columns = "".join(columns.split()).split(",")
        # Check for any kind of list, and rule out str which is a type of Sequence.
        if isinstance(columns, (Sequence, np.ndarray)) and not isinstance(columns, str):
            cols_set = set(columns)
            if len(cols_set) == 0:
                raise ValueError("Empty 'columns'.")
        else:
            raise TypeError(f"columns must be list-like, got a {type(columns)} instead.")
        # Selected columns must be defined in col_defs.
        col_defs_set = set(col_defs.keys())
        diff = cols_set - col_defs_set
        if len(diff) > 0 and not verbose:
            raise ValueError(f"Column(s) {diff} are not handled yet. Known columns are: {', '.join(col_defs_set)}")
        # No duplicate columns.
        if len(cols_set) < len(columns):
            logger.warning("columns contains duplicated values. Removing them.")
        # Sort the columns back to their input order.
        columns = sorted(cols_set, key=columns.index)
        # Can't deal with only the columns used to group the index.
        if needed is not None:
            if set(needed) >= set(columns):
                raise ValueError(f"Can't show only {' and/or '.join(needed)} columns. Add another column.")

        return columns

    def get_summary(self, scan=None, verbose=False, columns=None, add_columns=None, col_defs=None):
        """
        Create a summary of the input dataset as a `~pandas.DataFrame`.

        Parameters
        ----------
        scan : int or 2-tuple
            The scan(s) to use. A 2-tuple represents (beginning, ending) scans. Default: show all scans
        verbose : bool
            If verbose=False (default), the records are grouped by scan number and project id and aggregated
            according to the column. For example, the records for columns RESTFREQ, AZIMUTH and ELEVATIO are
            averaged for every scan. For columns IFNUM, PLNUM and FDNUM it counts the unique number of records.
            For column OBJECT it shows the value of the first record for the scan. For more details and a full
            list of the supported columns see `~dysh.fits.core.summary_column_definitions`.
            If True, list every record.
        columns : list or str
            List of columns for the output summary. If not set and `verbose=False`, the default list will contain SCAN, OBJECT,
            VELOCITY, PROC, PROCSEQN, RESTFREQ, DOPFREQ, IFNUM (# IF), PLNUM (# POL), INTNUM (# INT), FDNUM (# FEED), AZIMUTH,
            and ELEVATIO (ELEVATION).
            If not set and `verbose=True`, it will contain SCAN, OBJECT, VELOCITY, PROC, PROCSEQN, PROCSIZE, RESTFREQ,
            DOPFREQ, IFNUM, FEED, AZIMUTH, ELEVATIO, FDNUM, INTNUM, PLNUM, SIG, CAL, and DATE-OBS.
            If a string, multiple column names must be comma separated.
        add_columns : list
            List of columns to be added to the default `columns`.
            If `columns` is not None, then this will be ignored.
            If a string, multiple column names must be comma separated.
        col_defs : dict
            Dictionary with column definitions. See `~dysh.fits.core.summary_column_definitions` for the expected format.

        Returns
        -------
        summary : `~pandas.DataFrame`
            Summary of the data as a DataFrame.

        Raises
        ------
        TypeError
            If `column` is not a list.
        ValueError
            If one of the column names in `column` is not defined.
        KeyError
            If one of the column names in `column` is not part of the index.
        """

        # @todo set individual format options on output by
        # changing these to dicts(?)

        if col_defs is None:
            col_defs = core.summary_column_definitions()

        needed = ["PROJID", "BINTABLE", "SCAN"]

        # Initial handling of `add_columns` keyword.
        if add_columns is not None and columns is not None:
            logger.warning("Both 'columns' and 'add_columns' set. Will ignore 'add_columns'.")
            add_columns = []
        elif add_columns is None:
            add_columns = []
        else:
            add_columns = self._validate_summary_columns(add_columns, col_defs, verbose=verbose)

        # Deafult columns to show.
        if columns is None:
            if verbose:
                columns = [
                    "SCAN",
                    "OBJECT",
                    "VELOCITY",
                    "PROC",
                    "PROCSEQN",
                    "PROCSIZE",
                    "RESTFREQ",
                    "DOPFREQ",
                    "IFNUM",
                    "PLNUM",
                    "FDNUM",
                    "FEED",
                    "SIG",
                    "CAL",
                    "INTNUM",
                    "AZIMUTH",
                    "ELEVATIO",
                    "DATE-OBS",
                ]
            else:
                columns = [
                    "SCAN",
                    "OBJECT",
                    "VELOCITY",
                    "PROC",
                    "PROCSEQN",
                    "RESTFREQ",
                    "DOPFREQ",
                    "IFNUM",
                    "PLNUM",
                    "INTNUM",
                    "FDNUM",
                    "AZIMUTH",
                    "ELEVATIO",
                ]
        else:
            # Check that the user input won't break anything.
            columns = self._validate_summary_columns(columns, col_defs, needed, verbose)
        ocols = columns + add_columns  # Output columns.
        ocols = sorted(set(ocols), key=ocols.index)  # Remove duplicates preserving order.
        _columns = ocols.copy()
        for n in needed:
            try:
                _columns.remove(n)
            except ValueError:
                continue
        cols = _columns + needed  # All columns to fetch.
        # skipflags and flag_vegas from the constructor have been lost at this
        # point, but if the user skipped the index with index=False, then the flags
        # including vegas flags would have been skipped as well. So we can set them
        # here reasonably.  If index=True in the constructor, then the index was already
        # created and this method is a no-op.
        self._create_index_if_needed(skipflags=True, flag_vegas=False)

        # Define column types.
        col_dtypes = {k: v.type for k, v in col_defs.items() if k in cols}

        # make a copy here because we can't guarantee if this is a
        # view or a copy without it. See https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#returning-a-view-versus-a-copy
        df = self[cols].copy().astype(col_dtypes)

        # Scale columns.
        for cn in columns:
            try:
                if col_defs[cn].scale != 1:
                    df[cn] *= col_defs[cn].scale
            except KeyError:
                logger.warning(f"Column {cn} undefined. Please submit an issue.")
                continue

        if scan is not None:
            if isinstance(scan, int):
                scan = [scan]
            if len(scan) == 1:
                scan = [scan[0], scan[0]]
            df = self._select_scans(scan, df)

        if not verbose:
            # Short summary version.
            # Set column operations for aggregation.
            col_ops = {k: v.operation for k, v in col_defs.items() if k in _columns}
            # We have to reset the index and column types.
            df = df.groupby(needed).agg(col_ops).reset_index().astype(col_dtypes)
            # Post operations.
            col_post_ops = {k: v.post for k, v in col_defs.items() if k in _columns and v.post is not None}
            if len(col_post_ops) > 0:
                df[list(col_post_ops.keys())] = df.apply(col_post_ops)
            # Sort rows.
            df = df.sort_values(by=needed)
            # Keep only the columns to be shown.
            df = df[ocols]
            # Set column names.
            col_names = {k: v.name if k in _columns and v.name is not None else k for k, v in col_defs.items()}
            new_columns = [col_defs[c].name if col_defs[c].name is not None else c for c in ocols]
            df = df.rename(columns=col_names)
            df = df[new_columns]
        else:
            # Ensure column order is preserved.
            df = df[ocols]

        return df

    def summary(self, scan=None, verbose=False, max_rows=-1, show_index=False, columns=None, add_columns=None):
        """
        Show a summary of the `~dysh.fits.GBTFITSLoad` object.
        To retrieve the underlying `~pandas.DataFrame` use
        `~dysh.fits.GBTFITSLoad.get_summary()`.

        Parameters
        ----------
        scan : int or 2-tuple
            The scan(s) to use. A 2-tuple represents (beginning, ending) scans.
            Default: show all scans
        verbose : bool
            If True, list every record, otherwise return a compact summary.
            The compact summary averages some of the columns over scan number (e.g.,
            RESTFREQ, AZIMUTH, ELEVATIO), and lists the number of spectral windows (IFs),
            polarizations (# POL), feeds (# FEED), and integrations (# INT).
        max_rows : int or None
            Maximum number of rows to display. If less than the total number of rows, then
            the first `max_rows/2` and last `max_rows/2` rows will be shown, separated
            by ellipsis. If set to -1 (Default), the value found in the dysh configuration
            file for `summary_max_rows` will be used. Set to `None` for unlimited rows.
        show_index : bool
            Show index of the `~pandas.DataFrame`.
        columns : list or str
            List of columns for the output summary. If not set and `verbose=False`, the default list will contain SCAN,
            OBJECT, VELOCITY, PROC, PROCSEQN, RESTFREQ, DOPFREQ, IFNUM (# IF), PLNUM (# POL), INTNUM (# INT), FDNUM (# FEED),
            AZIMUTH, and ELEVATIO (ELEVATION).
            If not set and `verbose=True`, it will contain SCAN, OBJECT, VELOCITY, PROC, PROCSEQN, PROCSIZE, RESTFREQ,
            DOPFREQ, IFNUM, FEED, AZIMUTH, ELEVATIO, FDNUM, INTNUM, PLNUM, SIG, CAL, and DATE-OBS.
            If a string, multiple column names must be comma separated.
        add_columns : list or str
            List of columns to be added to the default `columns`.
            If `columns` is not None, then this will be ignored.
            If a string, multiple column names must be comma separated.
        """

        df = self.get_summary(scan=scan, verbose=verbose, columns=columns, add_columns=add_columns)

        if max_rows == -1:
            max_rows = conf.summary_max_rows
        max_cols = 1500

        show_dataframe(df, show_index=show_index, max_rows=max_rows, max_cols=max_cols)

    def velocity_convention(self, veldef):
        """Given the GBT VELDEF FITS string return the specutils
        velocity convention, e.g., "doppler_radio"

        Parameters
        ----------
            veldef : str
                The FITS header VELDEF string

        Returns
        -------
            convention : str
                The velocity convention
        """
        (convention, _frame) = decode_veldef(veldef)
        return convention

    def velocity_frame(self, veldef):
        """Given the GBT VELDEF FITS string return the
        velocity frame, e.g., "heliocentric".

        Parameters
        ----------
            veldef : str
                The FITS header VELDEF string

        Returns
        -------
            frame: str
                The velocity frame
        """
        (_convention, frame) = decode_veldef(veldef)
        return frame

    def _select_scans(self, scans, df):
        return df[(df["SCAN"] >= scans[0]) & (df["SCAN"] <= scans[1])]

    # @todo maybe move all selection/flag methods to sdfitsload after adding Selection/Flag
    # to sdfitsload
    # @todo maybe write a Delegator class to autopass to Selection.
    # See, e.g., https://michaelcho.me/article/method-delegation-in-python/
    @log_call_to_history
    def select(self, tag=None, check=False, **kwargs):
        """Add one or more exact selection rules, e.g., `key1 = value1, key2 = value2, ...`
        If `value` is array-like then a match to any of the array members will be selected.
        For instance `select(object=['3C273', 'NGC1234'])` will select data for either of those
        objects and `select(ifnum=[0,2])` will select IF number 0 or IF number 2.
        See `~dysh.util.selection.Selection`.

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
        self._selection.select(tag=tag, check=check, **kwargs)

    @log_call_to_history
    def select_range(self, tag=None, check=False, **kwargs):
        """
        Select a range of inclusive values for a given key(s).
        e.g., `key1 = (v1,v2), key2 = (v3,v4), ...`
        will select data  `v1 <= data1 <= v2, v3 <= data2 <= v4, ... `
        Upper and lower limits may be given by setting one of the tuple values
        to None. e.g., `key1 = (None,v1)` for an upper limit `data1 <= v1` and
        `key1 = (v1,None)` for a lower limit `data >=v1`.  Lower
        limits may also be specified by a one-element tuple `key1 = (v1,)`.

        For time values, :class:`~astropy.time.Time`, :class:`~np.datetime64` and :class:`~datetime.datetime` are supported.
        See `~dysh.util.selection.Selection`.

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
        self._selection.select_range(tag=tag, check=check, **kwargs)

    @log_call_to_history
    def select_within(self, tag=None, check=False, **kwargs):
        """
        Select a value within a plus or minus for a given key(s).
        e.g. `key1 = [value1,epsilon1], key2 = [value2,epsilon2], ...`
        Will select data

        `value1-epsilon1 <= data1 <= value1+epsilon1,`
        `value2-epsilon2 <= data2 <= value2+epsilon2,...`

        For time values, :class:`~astropy.time.Time`, :class:`~np.datetime64` and :class:`~datetime.datetime` are supported.
        See `~dysh.util.selection.Selection`.

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
        self._selection.select_within(tag=tag, check=check, **kwargs)

    @log_call_to_history
    def select_channel(self, channel, tag=None):
        """
        Select channels and/or channel ranges. These are NOT used in :meth:`final`
        but rather will be used to create a mask for calibration or
        flagging. Single arrays/tuples will be treated as channel lists;
        nested arrays will be treated as ranges, for instance

        .. code::

            # selects channels 1 and 10
            select_channel([1,10])
            # selects channels 1 thru 10 inclusive
            select_channel([[1,10]])
            # select channel ranges 1 thru 10 and 47 thru 56 inclusive, and channel 75
            select_channel([[1,10], [47,56], 75)])
            # tuples also work
            select_channel(((1,10), [47,56], 75))


        See `~dysh.util.selection.Selection`.

        Parameters
        ----------
        channel : number, or array-like
            The channels to select

        Returns
        -------
        None.
        """

        self._selection.select_channel(tag=tag, channel=channel)

    @log_call_to_history
    def clear_selection(self):
        """Clear all selections for these data"""
        self._selection.clear()

    @log_call_to_history
    def flag(self, tag=None, check=False, **kwargs):
        """Add one or more exact flag rules, e.g., `key1 = value1, key2 = value2, ...`
        If `value` is array-like then a match to any of the array members will be flagged.
        For instance `flag(object=['3C273', 'NGC1234'])` will select data for either of those
        objects and `flag(ifnum=[0,2])` will flag IF number 0 or IF number 2.  Channels for selected data
        can be flagged using keyword `channel`, e.g., `flag(object='MBM12',channel=[0,23])`
        will flag channels 0 through 23 *inclusive* for object MBM12.
        See `~dysh.util.selection.Flag`.

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
        self._flag.flag(tag=tag, check=check, **kwargs)

    @log_call_to_history
    def flag_range(self, tag=None, check=False, **kwargs):
        """
        Flag a range of inclusive values for a given key(s).
        e.g., `key1 = (v1,v2), key2 = (v3,v4), ...`
        will select data  `v1 <= data1 <= v2, v3 <= data2 <= v4, ...`

        Upper and lower limits may be given by setting one of the tuple values
        to None. e.g., `key1 = (None,v1)` for an upper limit `data1 <= v1` and
        `key1 = (v1,None)` for a lower limit `data >=v1`.  Lower
        limits may also be specified by a one-element tuple `key1 = (v1,)`.

        For time values, :class:`~astropy.time.Time`, :class:`~np.datetime64` and :class:`~datetime.datetime` are supported.

        See `~dysh.util.selection.Flag`.

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
        self._flag.flag_range(tag=tag, check=check, **kwargs)

    @log_call_to_history
    def flag_within(self, tag=None, check=False, **kwargs):
        """
        Flag a value within a plus or minus for a given key(s).
        e.g. `key1 = [value1,epsilon1], key2 = [value2,epsilon2], ...`
        Will select data

        `value1-epsilon1 <= data1 <= value1+epsilon1,`
        `value2-epsilon2 <= data2 <= value2+epsilon2,...`

        For time values, :class:`~astropy.time.Time`, :class:`~np.datetime64` and :class:`~datetime.datetime` are supported.

        See `~dysh.util.selection.Flag`.

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
        self._flag.flag_within(tag=tag, check=check, **kwargs)

    @log_call_to_history
    def flag_channel(self, channel, tag=None):
        """
        Flag channels and/or channel ranges. These will be used to create a mask for
        data when calibrating (see :meth:`apply_flags`). Single arrays/tuples will be treated as channel lists;
        nested arrays will be treated as ranges, for instance

        .. code::

            # flag channel 128
            flag_channel(128)
            # flags channels 1 and 10
            flag_channel([1,10])
            # flags channels 1 thru 10 inclusive
            flag_channel([[1,10]])
            # flags channel ranges 1 thru 10 and 47 thru 56 inclusive, and channel 75
            flag_channel([[1,10], [47,56], 75)])
            # tuples also work
            flag_channel(((1,10), [47,56], 75))


        See `~dysh.util.selection.Flag`.

        Parameters
        ----------
        channel : number, or array-like
            The channels to flag

        Returns
        -------
        None.
        """
        self._flag.flag_channel(tag=tag, channel=channel)

    def is_vegas(self):
        """Check if these data appear to use the VEGAS backend

        Returns
        -------
            True if FITS HEADER Keyword INSTRUME or BACKEND is present and equals 'VEGAS', False otherwise
        """
        if "INSTRUME" in self._selection:
            instrument = str(next(iter(set(self["INSTRUME"])))).upper()
        else:
            instrument = ""
        if "BACKEND" in self._selection:
            backend = str(next(iter(set(self["BACKEND"])))).upper()
        else:
            backend = ""
        if instrument == "VEGAS" or backend == "VEGAS":
            return True
        else:
            return False

    @log_call_to_history
    def flag_vegas_spurs(self, flag_central=False):
        """
        Flag VEGAS SPUR channels.

        Parameters
        ----------
        flag_central : bool, optional
            Whether to flag the central VEGAS spur location or not.
            The GBO SDFITS writer by default replaces the value at the central SPUR with the average of the
            two adjacent channels, and hence the central channel is not typically flagged.

        Returns
        -------
        None.

        """
        if not self.is_vegas():
            logger.warning(
                "This does not appear to be VEGAS data. Check if FITS Header keywords 'INSTRUME' or 'BACKEND' are present and equal 'VEGAS'. No channels will be flagged."
            )
            return
        try:
            df = self._selection.groupby(["FITSINDEX", "BINTABLE"])
            for _i, ((fi, bi), g) in enumerate(df):
                vsprval = g["VSPRVAL"].to_numpy()
                vspdelt = g["VSPDELT"].to_numpy()
                vsprpix = g["VSPRPIX"].to_numpy()
                rows = g["ROW"].to_numpy()
                maxnchan = self._sdf[fi].nchan(bi) - 1
                spurs = calc_vegas_spurs(vsprval, vspdelt, vsprpix, maxnchan, flag_central)
                if spurs.shape[0] != len(rows):
                    raise ValueError(f"spurs array length {spurs.shape[0]} != selected number of rows {len(rows)}")
                else:
                    p = []
                    # there should be a way to do this without a for loop.
                    for a in spurs:
                        mask = np.full(maxnchan + 1, False)
                        mask[a[~a.mask]] = True
                        p.append(mask)
                    self._sdf[fi]._additional_channel_mask[bi][rows] |= np.array(p)
        except KeyError as k:
            logger.warning(
                f"Can't determine VEGAS spur locations because one or more VSP keywords are missing from the FITS header {k}"
            )

    @log_call_to_history
    def apply_flags(self):
        """
        Set the channel flags according to the rules specified in the `flags` attribute.
        This sets numpy masks in the underlying `SDFITSLoad` objects.

        Returns
        -------
        None.

        """
        # Loop over the dict of flagged channels, which
        # have the same key as the flag rules.
        # For all SDFs in each flag rule, set the flag mask(s)
        # for their rows.  The index of the sdf._flagmask array is the bintable index
        for key, chan in self._flag._flag_channel_selection.items():
            selection = self._flag.get(key)
            # chan will be a list or a list of lists
            # If it is a single list, it is just a list of channels
            # if it is list of lists, then it is upper lower inclusive
            dfs = selection.groupby(["FITSINDEX", "BINTABLE"])
            # the dict key for the groups is a tuple (fitsindex,bintable)
            for _i, ((fi, bi), g) in enumerate(dfs):
                chan_mask = convert_array_to_mask(chan, self._sdf[fi].nchan(bi))
                rows = g["ROW"].to_numpy()
                logger.debug(f"Applying {chan} to {rows=}")
                logger.debug(f"{np.where(chan_mask)}")
                self._sdf[fi]._flagmask[bi][rows] |= chan_mask
        # now any additional channel flags, i.e. VEGAS flags
        self._apply_additional_flags()

    def _apply_additional_flags(self):
        """apply the additional channel flags created by, e.g., flag_vegas"""
        for k in self._sdf:
            if k._additional_channel_mask is not None and k._flagmask is not None:
                k._flagmask |= k._additional_channel_mask

    @log_call_to_history
    def clear_flags(self):
        """Clear all flags for these data"""
        for sdf in self._sdf:
            sdf._init_flags()
        self._flag.clear()

    def _create_index_if_needed(self, skipflags=False, flag_vegas=True):
        """
        Parameters
        ----------
        skipflags : bool, optional
            If True, do not read any flag files associated with these data.  The default is False.

        Returns
        -------
        None.

        """

        if self._selection is not None:
            return
        i = 0
        df = None
        if self._selection is None:
            for s in self._sdf:
                if s._index is None:
                    s.create_index()
                # add a FITSINDEX column
                s._index["FITSINDEX"] = i * np.ones(len(s._index), dtype=int)
                if df is None:
                    df = s._index
                else:
                    df = pd.concat([df, s._index], axis=0, ignore_index=True)
                i = i + 1
        self._selection = Selection(df)
        self._flag = Flag(df)
        self._construct_procedure()
        self._construct_integration_number()

        if flag_vegas and self.is_vegas():
            self.flag_vegas_spurs()

        if skipflags:
            return

        # for directories with multiple FITS files and possibly multiple FLAG files
        # we have to ensure the right flag file goes with the right FITS tile.
        # The GBT convention is the same filename with '.flag' instead of '.fits'.
        # We construct the flagfile and also pass in FITSINDEX column to ensure
        # only the data associated with that file are flagged.
        found_flags = False
        for s in self._sdf:
            p = Path(s.filename)
            flagfile = p.with_suffix(".flag")
            if flagfile.exists():
                fi = uniq(s["FITSINDEX"])[0]
                self.flags.read(flagfile, fitsindex=fi, ignore_vegas=flag_vegas)
                found_flags = True
        if found_flags and len(self.flags._table) != 0:
            logger.info("Flags were created from existing flag files. Use GBTFITSLoad.flags.show() to see them.")

    def _construct_procedure(self):
        """
        Construct the procedure string (PROC) from OBSMODE and add it to the index (i.e., a new SDFITS column).
        OBSTYPE and SUBOBSMODE are also created here.  OBSMODE has the form like 'PROC:OBSTYPE:SUBOBSMODE', e.g.
        OnOff:PSWITCHON:TPWCAL.

        """
        if self._selection is None:
            warnings.warn("Couldn't construct procedure string: index is not yet created.")  # noqa: B028
            return
        if "OBSMODE" not in self._index:
            warnings.warn("Couldn't construct procedure string: OBSMODE is not in index.")  # noqa: B028
            return
        df = self["OBSMODE"].str.split(":", expand=True)
        for obj in [self._index, self._flag]:
            obj["PROC"] = df[0]
            # Assign these to something that might be useful later,
            # since we have them
            obj["OBSTYPE"] = df[1]
            obj["SUBOBSMODE"] = df[2]
        for sdf in self._sdf:  # Note: sdf._index is a Dataframe, not a Selection
            df = sdf._index["OBSMODE"].str.split(":", expand=True)
            sdf._index["PROC"] = df[0]
            sdf._index["OBSTYPE"] = df[1]
            sdf._index["SUBOBSMODE"] = df[2]

    def _construct_integration_number(self):
        """Construct the integration number (INTNUM) for all scans and add it to the index (i.e., a new SDFITS column)
        Integration number starts at zero and is incremented when the DATE-OBS changes. It resets to
        zero when the scan number changes.
        """
        if self._index is None:
            warnings.warn("Couldn't construct integration number: index is not yet created.", stacklevel=2)
            return

        # check it hasn't been constructed before.
        if "INTNUM" in self._index:
            return
        # check that GBTIDL didn't write it out at some point.
        if "INT" in self._index:
            # This is faster than using
            # self._selection = Selection(self._selection.rename(columns={"INT": "INTNUM"}))
            # Keep, unless we find a faster way.
            self._index.rename(columns={"INT": "INTNUM"}, inplace=True)  # noqa: PD002
            for s in self._sdf:
                s._rename_binary_table_column("int", "intnum")
            return

        intnumarray = np.empty(len(self._selection), dtype=int)
        # Leverage pandas to group things by scan and observing time.
        dfs = self._selection.groupby(["SCAN"])
        for _, group in dfs:
            # Group by FITSINDEX since different banks can have different time stamps.
            dfsf = group.groupby("FITSINDEX")
            for _, fg in dfsf:
                dfsft = fg.groupby("DATE-OBS")
                intnums = np.arange(0, len(dfsft.groups))
                for i, (_, g) in enumerate(dfsft):
                    idx = g.index
                    intnumarray[idx] = intnums[i]
        self._selection["INTNUM"] = intnumarray
        self._flag["INTNUM"] = intnumarray

    def _normalize_channel_range(self, channel: list) -> list | None:
        """Ensure channel range is [first,last] for calibration"""
        if channel is None and self._selection._channel_selection is None:
            return None
        elif channel is not None and self._selection._channel_selection is not None:
            raise ValueError(
                f"A channel selection was previously made: {self._selection._channel_selection}.  Clear that selection before attempting to calibrate with a different channel range."
            )
        elif channel is None and self._selection._channel_selection is not None:
            return get_valid_channel_range(self._selection._channel_selection)
        else:
            return get_valid_channel_range(channel)

    def _common_selection(self, ifnum, plnum, fdnum, **kwargs):
        """Do selection and flag application common to all calibration methods.
        Flags are not applied unless selection results in non-zero length data selection.

        Parameters
        ----------
        fdnum: int
            The feed number
        ifnum : int
            The intermediate frequency (IF) number
        plnum : int
            The polarization number
        kwargs : dict
            Additional selections

        Returns
        -------
        scan_df : tuple
            A tuple consisting of a list of scan numbers selected and a `~pandas.DataFrame` of the selection.
        """
        kwargs = keycase(kwargs)
        apply_flags = kwargs.pop("APPLY_FLAGS", True)

        if len(self._selection._selection_rules) > 0:
            _final = self._selection.final
        else:
            _final = self._index
        scans = kwargs.get("SCAN", None)
        if scans is None:
            scans = uniq(_final["SCAN"])
        elif isinstance(scans, int):
            scans = list([scans])
        if "REF" in kwargs:
            scans.append(kwargs.pop("REF"))
            scans = uniq(scans)
        preselected = {}
        preselected["SCAN"] = scans
        preselected["FDNUM"] = fdnum
        preselected["IFNUM"] = ifnum
        preselected["PLNUM"] = plnum
        for k, v in preselected.items():
            if k not in kwargs:
                kwargs[k] = v
        # For PS and Nod scans we must find the full pairs of ON/OFF scans since the
        # user may have input only the ONs or OFFs.
        # _common_scan_list_selection returns a dict of scan numbers that contain all the ON/OFF pairs
        # This will replace the input list of scans.
        scans_to_add = []
        if "PROCKEY" in kwargs:
            missing = self._common_scan_list_selection(scans, _final, kwargs["PROCKEY"], kwargs["PROCVALS"], check=True)
            scans_to_add = set(missing["ON"]).union(missing["OFF"])
            logger.debug(f"after removing preselected {preselected['SCAN']}, scans_to_add={scans_to_add}")
            kwargs.pop("PROCKEY")
            kwargs.pop("PROCVALS")
        if len(scans_to_add) != 0:
            # add a rule selecting the missing scans :-)
            logger.debug(f"adding rule scan={scans_to_add}")
            kwargs["SCAN"] = list(scans_to_add)
        if len(_final[_final["SCAN"].isin(scans)]) == 0:
            raise ValueError(f"Scans {scans} not found in selected data")
        ps_selection = copy.deepcopy(self._selection)
        # now downselect with any additional kwargs
        ps_selection._select_from_mixed_kwargs(**kwargs)
        _sf = ps_selection.final
        # now remove rows that have been entirely flagged
        if apply_flags:
            _sf = eliminate_flagged_rows(_sf, self.flags.final)
        if len(_sf) == 0:
            raise Exception("Didn't find any unflagged scans matching the input selection criteria.")
        # Don't apply flags until we are sure that selection succeeded
        if apply_flags:
            self.apply_flags()
        return (scans, _sf)

    def info(self):
        """Return information on the HDUs contained in this object. See :meth:`~astropy.HDUList/info()`"""
        for s in self._sdf:
            s.info()

    @log_call_to_result
    def gettp(
        self,
        fdnum: int,
        ifnum: int,
        plnum: int,
        sig: bool | None = None,
        cal: bool | None = None,
        calibrate: bool = True,
        apply_flags: bool = True,
        t_sys=None,
        t_cal=None,
        channel: list | None = None,
        vane=None,
        **kwargs,
    ):
        """
        Get a total power scan, optionally calibrating it.

        Parameters
        ----------
        fdnum: int
            The feed number
        ifnum : int
            The intermediate frequency (IF) number
        plnum : int
            The polarization number
        sig : bool or None
            True to use only integrations where signal state is True, False to use reference state (signal state is False). None to use all integrations.
        cal: bool or None
            True to use only integrations where calibration (diode) is on, False if off. None to use all integrations regardless calibration state.
            The system temperature will be calculated from both states regardless of the value of this variable.
        calibrate: bool
            whether or not to calibrate the data.  If `True`, the data will be (calon + caloff)*0.5, otherwise it will be SDFITS row data.
            Default:True
        apply_flags : boolean, optional.  If True, apply flags before calibration.
            See :meth:`apply_flags`. Default: True
        t_sys : float
            System temperature. If provided, it overrides the value computed using the noise diode.
            If no noise diode is fired, and `t_sys=None`, then the column "TSYS" will be used instead.
        t_cal : None or float
            Noise diode temperature. If provided, this value is used instead of the value found in the
            TCAL column of the SDFITS file. If no value is provided, default, then the TCAL column is
            used.
        channel: list or None
            An inclusive list of `[firstchan, lastchan]` to use in the calibration. The channel list is zero-based. If provided,
            only data channels in the inclusive range `[firstchan,lastchan]` will be used. If a reference spectrum has been given, it will also be
            trimmed to `[firstchan,lastchan]` before any smoothing. If channels have already been selected through
            :meth:`GBTFITSLoad.select_channel`, a ValueError will be raised.
        vane : None
            Used to suppress info message about use of TSYS column in case this is being used to make a `~dysh.spectra.vane.VaneSpectrum`.
        **kwargs : dict
            Optional additional selection  keyword arguments, typically
            given as key=value, though a dictionary works too.
            e.g., `source="NGC132", intnum=range(20)` etc.

        Returns
        -------
        data : `~dysh.spectra.scan.ScanBlock`
            A ScanBlock containing one or more `~dysh.spectra.scan.TPScan`

        """
        _channel = self._normalize_channel_range(channel)
        (scans, _sf) = self._common_selection(fdnum=fdnum, ifnum=ifnum, plnum=plnum, apply_flags=apply_flags, **kwargs)
        tsys = _parse_tsys(t_sys, scans)
        _tsys = None
        _tcal = t_cal
        _bintable = kwargs.get("bintable", None)
        TF = {True: "T", False: "F"}
        scanblock = ScanBlock()
        calrows = {}
        for scan in scans:
            _sifdf = select_from("SCAN", scan, _sf)
            if len(_sifdf) == 0:
                continue
            dfcalT = select_from("CAL", "T", _sifdf)
            dfcalF = select_from("CAL", "F", _sifdf)
            calrows["ON"] = dfcalT["ROW"].to_numpy()
            calrows["OFF"] = dfcalF["ROW"].to_numpy()
            if len(calrows["ON"]) != len(calrows["OFF"]):
                if len(calrows["ON"]) > 0:
                    raise Exception(f"unbalanced calrows {len(calrows['ON'])} != {len(calrows['OFF'])}")
            # sig and cal are treated specially since
            # they are not in kwargs and in SDFITS header
            # they are not booleans but chars
            if sig is not None:
                _sifdf = select_from("SIG", TF[sig], _sifdf)
            if _bintable is None:
                _bintable = self._get_bintable(_sifdf)
            if t_cal is not None:
                _tcal = t_cal
            else:
                _tcal = self._get_tcal(dfcalF["TCAL"])
            if len(calrows["ON"]) == 0:
                if tsys is None:
                    _tsys = dfcalF["TSYS"].to_numpy()
                    if vane is None:
                        logger.info("Using TSYS column")
                    logger.debug(f"Scan: {scan}")
            # Use user provided system temperature.
            if tsys is not None:
                _tsys = tsys[scan][0]
            # The rows with the selected sig state and all cal states.
            tprows = _sifdf["ROW"].to_numpy()
            logger.debug(f"TPROWS len={len(tprows)}")
            logger.debug(f"CALROWS on len={len(calrows['ON'])}")
            logger.debug(f"fitsindex={_sifdf['FITSINDEX'].iloc[0]}")
            if len(tprows) == 0:
                continue
            if "TSCALE" in _sifdf:
                tscale = _sifdf["TSCALE"].unique()
            else:
                tscale = ["Raw"]
            if len(tscale) > 1:
                raise ValueError(
                    f"More than one TSCALE value in the previously-calibrated input file {tscale}; can't create a TPScan."
                )
            g = TPScan(
                self._sdf[_sifdf["FITSINDEX"].iloc[0]],
                scan,
                sig,
                cal,
                tprows,
                calrows,
                fdnum=fdnum,
                ifnum=ifnum,
                plnum=plnum,
                bintable=_bintable,
                calibrate=calibrate,
                apply_flags=apply_flags,
                tsys=_tsys,
                tcal=_tcal,
                tscale=tscale[0],
                channel=_channel,
            )
            tscalefac = _sifdf.get("TSCALFAC", None)
            if tscalefac is not None:
                # the data were previously calibrated, preserve the scale factor
                g._tscale_fac = np.array(tscalefac)
            g.merge_commentary(self)
            scanblock.append(g)
            # Reset these variables for the next scan.
            _tsys = tsys
            _tcal = t_cal
            _bintable = kwargs.get("bintable", None)
        if len(scanblock) == 0:
            raise Exception("Didn't find any scans matching the input selection criteria.")
        scanblock.merge_commentary(self)
        return scanblock
        # end of gettp()

    @log_call_to_result
    def getsigref(
        self,
        scan: int | list | np.ndarray,
        ref: int | Spectrum,
        fdnum: int,
        ifnum: int,
        plnum: int,
        calibrate: bool = True,
        smoothref: int = 1,
        apply_flags: str = True,
        units: str = "ta",
        zenith_opacity: float | None = None,
        weights="tsys",
        t_sys=None,
        t_cal=None,
        nocal: bool = False,
        ap_eff: float | None = None,
        surface_error: Quantity | None = None,
        channel: list | None = None,
        vane: int | VaneSpectrum | None = None,
        t_atm: float | None = None,
        t_bkg: float | None = None,
        t_warm: float | None = None,
        **kwargs,
    ) -> ScanBlock:
        r"""
        Retrieve and calibrate position-switched data using a custom reference scan.  Also known as `Flexible Off.`

        Parameters
        ----------
        scan : int or list or `numpy.array`
            The signal scan numbers to calibrate
        ref : int or Spectrum
            The reference scan number or a `~dysh.spectra.spectrum.Spectrum` object.  If an integer is given,
            the reference spectrum will be the total power time-averaged spectrum using the weights given.
            If `channel` is given, the reference spectrum will be trimmed to the `channel` range before calibration.
        fdnum : int
            The feed number.
        ifnum : int
            The intermediate frequency (IF) number.
        plnum : int
            The polarization number.
        calibrate : boolean, optional
            Calibrate the scans. The default is True.
        smooth_ref : int, optional
            If >1 smooth the reference with a boxcar kernel with a width of `smooth_ref` channels. The default is to not smooth the reference.
        apply_flags : boolean, optional
            If True, apply flags before calibration.
            See :meth:`apply_flags`. Default: True
        units : str, optional
            The brightness scale unit for the output scan, must be one of (case-insensitive)
                    - 'ta'   : Antenna Temperature
                    - 'ta*'  : Antenna temperature corrected to above the atmosphere
                    - 'flux' : flux density in Jansky
            If 'ta*' or 'flux' the zenith opacity must also be given. Default: 'ta'
        zenith_opacity : float, optional
            The zenith opacity to use in calculating the scale factors for the integrations.  Default: None
        t_sys : float, optional
            If given, this is the system temperature in Kelvin. It overrides the values calculated using the noise diodes.
            If not given, and signal and reference are scan numbers, the system temperature will be calculated from the reference
            scan and the noise diode. If not given, and the reference is a `Spectrum`, the reference system temperature as given
            in the metadata header will be used. The default is to use the noise diode or the metadata, as appropriate.
            If `vane` is provided, then `t_sys` will be ignored and `vane` will be used to derive the system temperature.
        t_cal : None or float
            Noise diode temperature. If provided, this value is used instead of the value found in the
            TCAL column of the SDFITS file. If no value is provided, default, then the TCAL column is
            used. If `t_sys` is provided, `t_cal` will be ignored.
        nocal : bool, optional
            Is the noise diode being fired? False means the noise diode was firing.
            By default it will figure this out by looking at the "CAL" column.
            It can be set to True to override this.
        ap_eff : float or None
            Aperture efficiency o be used when scaling data to brightness temperature of flux. The provided aperture
            efficiency must be a number between 0 and 1.  If None, `dysh` will calculate it as described in
            :meth:`~GBTGainCorrection.aperture_efficiency`. Only one of `ap_eff` or `surface_error`
            can be provided.
        surface_error : `~astropy.units.Quantity` or None
            Surface rms error, in units of length (typically microns), to be used in the Ruze formula when calculating the
            aperture efficiency.  If None, `dysh` will use the known GBT surface error model.  Only one of `ap_eff` or `surface_error`
            can be provided.
        channel: list or None
            An inclusive list of `[firstchan, lastchan]` to use in the calibration. The channel list is zero-based. If provided,
            only data channels in the inclusive range `[firstchan,lastchan]` will be used. If a reference spectrum has been given, it will also be
            trimmed to `[firstchan,lastchan]`.  System temperature calculation will use 80% of the trimmed channel range.  If channels have already been selected through
            :meth:`GBTFITSLoad.select_channel`, a ValueError will be raised.
        vane : int or `~dysh.spectra.vane.VaneSpectrum` or None
            Vane scalibration scan. This will be used to derive the system temperature.
            If provided, `t_sys` will be ignored.
        t_atm : float or None
            Atmospheric temperature in K. If `vane` is a `~dysh.spectra.vane.VaneSpectrum` it won't be used.
            If `vane` is an `int`, then the resulting `~dysh.spectra.vane.VaneSpectrum` will use this value for
            the atmospheric temperature. If not provided and `vane` is an `int`, `~dysh.spectra.vane.VaneSpectrum` will try to fetch a
            value from the GBO weather forecast script (only available at GBO).
        t_bkg : float or None
            Background temperature in K. If `vane` is a `~dysh.spectra.vane.VaneSpectrum` it won't be used.
            If `vane` is an `int`, then the resulting `~dysh.spectra.vane.VaneSpectrum` will use this value for
            the background temperature. If not provided, it will take a default value of 2.725 K, i.e., the CMB at 3 mm.
        t_warm : float or None
            Vane temperature in K. If `vane` is a `~dysh.spectra.vane.VaneSpectrum` it won't be used.
            If `vane` is an `int`, then the resulting `~dysh.spectra.vane.VaneSpectrum` will use this value for the vane temperature.
            If not provided and `vane` is an `int`, it will take the value found in the "TWARM" column of the SDFITS.
        **kwargs : dict
            Optional additional selection keyword arguments, typically
            given as key=value, though a dictionary works too.
            e.g., `source='NGC123', ` etc.

        Raises
        ------
        Exception
            If scans matching the selection criteria are not found.

        Returns
        -------
        scanblock : `~dysh.spectra.scan.ScanBlock`
            ScanBlock containing one or more `~dysh.spectra.scan.PSScan`.

        """
        ScanBase._check_tscale(units)
        ScanBase._check_gain_factors(ap_eff, surface_error)
        self._check_vane_and_t_sys_args(vane, t_sys)

        if vane is not None:
            vane, units, requested_units, zenith_opacity = self._vane_setup(
                vane, fdnum, ifnum, plnum, units, zenith_opacity, t_warm, t_atm, t_bkg, t_cal, apply_flags
            )

        if units.lower() != "ta" and zenith_opacity is None and ap_eff is None:
            raise ValueError("Can't scale the data without a valid zenith opacity")
        if not isinstance(ref, int) and not isinstance(ref, Spectrum):
            raise TypeError("Reference scan ('ref') must be either an integer scan number or a Spectrum object")
        if isinstance(scan, Spectrum):
            raise TypeError(
                "Spectrum object not allowed for 'scan'.  You can use Spectrum arithmetic if both 'scan' and 'ref' are Spectrum objects"
            )
        _channel = self._normalize_channel_range(channel)
        if t_sys is not None and t_cal is not None:
            warnings.warn("Both t_cal and t_sys were set. Only t_sys will be used.", stacklevel=2)

        scanlist = {}
        if isinstance(scan, int):
            scan = [scan]
        elif isinstance(scan, np.ndarray):
            scan = list(scan)
        (scans, _sf) = self._common_selection(
            fdnum=fdnum,
            ifnum=ifnum,
            plnum=plnum,
            apply_flags=apply_flags,
            scan=scan,
            **kwargs,
        )
        tsys = _parse_tsys(t_sys, scans)
        _tsys = None
        _tcal = t_cal
        _nocal = nocal
        _bintable = kwargs.get("bintable", None)
        scanlist["ON"] = scans
        scanlist["OFF"] = [None] * len(scans)
        if isinstance(ref, int):
            # make an average reference spectrum
            refspec = self.gettp(
                scan=ref,
                fdnum=fdnum,
                ifnum=ifnum,
                plnum=plnum,
                calibrate=calibrate,
                apply_flags=apply_flags,
                t_cal=t_cal,
                channel=_channel,
                vane=vane,
                **kwargs,
            ).timeaverage(weights=weights)
        else:
            refspec = ref._copy()  # Needs to be a copy since we will change it afterwards.

        if t_cal is not None:
            _tcal = t_cal
            # Scale the system temperature.
            refspec.meta["TSYS"] *= _tcal / refspec.meta["TCAL"]
            refspec.meta["TCAL"] = _tcal
        else:
            _tcal = refspec.meta["TCAL"]

        # Check if `refspec` has a system temperature.
        if tsys is None:
            tsys = self._get_refspec_tsys(refspec)
            tsys = _parse_tsys(tsys, scans)  # Put it in a known format.

        scanblock = ScanBlock()
        for i in range(len(self._sdf)):
            _df = select_from("FITSINDEX", i, _sf)
            if len(_df) == 0:
                continue
            if len(scanlist["ON"]) == 0 or len(scanlist["OFF"]) == 0:
                logger.debug(f"scans {scans} not found, continuing")
                continue
            rows = {}

            for on, off in zip(scanlist["ON"], scanlist["OFF"], strict=False):
                _ondf = select_from("SCAN", on, _df)
                _offdf = select_from("SCAN", off, _df)
                if _bintable is None:
                    _bintable = self._get_bintable(_ondf)
                rows["ON"] = list(_ondf["ROW"])
                rows["OFF"] = list(_offdf["ROW"])
                for key in rows:
                    if len(rows[key]) == 0 and off is not None:
                        raise Exception(f"{key} scans not found in scan list {scans}")
                # do not pass scan list here. We need all the cal rows. They will
                # be intersected with scan rows in PSScan
                calrows = {}
                dfcalT = select_from("CAL", "T", _df)
                dfcalF = select_from("CAL", "F", _df)
                calrows["ON"] = list(dfcalT["ROW"])
                calrows["OFF"] = list(dfcalF["ROW"])
                d = {"ON": on, "OFF": off}
                if len(calrows["ON"]) == 0 or nocal:
                    _nocal = True
                    if tsys is None:
                        dfoncalF = select_from("CAL", "F", _ondf)
                        _tsys = dfoncalF["TSYS"].to_numpy()
                        if vane is None:
                            logger.info("Using TSYS column")
                # Use user provided system temperature.
                if tsys is not None:
                    try:
                        _tsys = tsys[on][0]
                    except KeyError:
                        _tsys = tsys[off][0]
                g = PSScan(
                    self._sdf[i],
                    scan=d,
                    scanrows=rows,
                    calrows=calrows,
                    fdnum=fdnum,
                    ifnum=ifnum,
                    plnum=plnum,
                    bintable=_bintable,
                    calibrate=calibrate,
                    smoothref=smoothref,
                    apply_flags=apply_flags,
                    tscale=units,
                    zenith_opacity=zenith_opacity,
                    refspec=refspec,
                    nocal=_nocal,
                    tsys=_tsys,
                    tcal=_tcal,
                    ap_eff=ap_eff,
                    surface_error=surface_error,
                    channel=_channel,
                    vane=vane,
                )
                g._refscan = ref
                # If calibrated with a vane change the units (so ugly >.<).
                if vane is not None:
                    self._set_scale_vane(g, requested_units, zenith_opacity)
                g.merge_commentary(self)
                scanblock.append(g)
                # Reset these variables in case they change for the next scan.
                _nocal = nocal
                _tsys = None
                _bintable = kwargs.get("bintable", None)

        if len(scanblock) == 0:
            raise Exception("Didn't find any scans matching the input selection criteria.")
        scanblock.merge_commentary(self)
        return scanblock
        # end of getsigref()

    @log_call_to_result
    def getps(
        self,
        fdnum: int,
        ifnum: int,
        plnum: int,
        calibrate: bool = True,
        smoothref: int = 1,
        apply_flags: str = True,
        units: str = "ta",
        zenith_opacity: float | None = None,
        t_sys=None,
        nocal=False,
        t_cal=None,
        ap_eff: float | None = None,
        surface_error: Quantity | None = None,
        channel: list | None = None,
        vane: int | VaneSpectrum | None = None,
        t_atm: float | None = None,
        t_bkg: float | None = None,
        t_warm: float | None = None,
        **kwargs,
    ) -> ScanBlock:
        """
        Retrieve and calibrate position-switched data.

        Parameters
        ----------
        fdnum : int
            The feed number.
        ifnum : int
            The intermediate frequency (IF) number.
        plnum : int
            The polarization number.
        calibrate : boolean, optional
            Calibrate the scans. The default is True.
        smooth_ref : int, optional
            The number of channels in the reference to boxcar smooth prior to calibration.
        apply_flags : boolean, optional.  If True, apply flags before calibration.
            See :meth:`apply_flags`. Default: True
        units : str, optional
            The brightness scale unit for the output scan, must be one of (case-insensitive)
                    - 'ta'   : Antenna Temperature
                    - 'ta*'  : Antenna temperature corrected to above the atmosphere
                    - 'flux' : flux density in Jansky
            If 'ta*' or 'flux' the zenith opacity must also be given. Default: 'ta'
        zenith_opacity: float, optional
            The zenith opacity to use in calculating the scale factors for the integrations.  Default:None
        t_sys : float
            System temperature. If provided, it overrides the value computed using the noise diode.
            If no noise diode is fired, and `t_sys=None`, then the column "TSYS" will be used instead.
            If `vane` is provided, then `t_sys` will be ignored and `vane` will be used to derive the system temperature.
        t_cal : None or float
            Noise diode temperature. If provided, this value is used instead of the value found in the
            TCAL column of the SDFITS file. If no value is provided, default, then the TCAL column is
            used.
        nocal : bool, optional
            Is the noise diode being fired? False means the noise diode was firing.
            By default it will figure this out by looking at the "CAL" column.
            It can be set to True to override this. Default: False
        ap_eff : float or None
            Aperture efficiency o be used when scaling data to brightness temperature of flux. The provided aperture
            efficiency must be a number between 0 and 1.  If None, `dysh` will calculate it as described in
            :meth:`~GBTGainCorrection.aperture_efficiency`. Only one of `ap_eff` or `surface_error`
            can be provided.
        surface_error : `~astropy.units.Quantity` or None
            Surface rms error, in units of length (typically microns), to be used in the Ruze formula when calculating the
            aperture efficiency.  If None, `dysh` will use the known GBT surface error model.  Only one of `ap_eff` or `surface_error`
            can be provided.
        channel: list or None
            An inclusive list of `[firstchan, lastchan]` to use in the calibration. The channel list is zero-based. If provided,
            only data channels in the inclusive range `[firstchan,lastchan]` will be used. If a reference spectrum has been given, it will also be
            trimmed to `[firstchan,lastchan]`.  System temperature calculation will use 80% of the trimmed channel range.  If channels have already been selected through
            :meth:`GBTFITSLoad.select_channel`, a ValueError will be raised.
        vane : int or `~dysh.spectra.vane.VaneSpectrum` or None
            Vane scalibration scan. This will be used to derive the system temperature.
            If provided, `t_sys` will be ignored.
        t_atm : float or None
            Atmospheric temperature in K. If `vane` is a `~dysh.spectra.vane.VaneSpectrum` it won't be used.
            If `vane` is an `int`, then the resulting `~dysh.spectra.vane.VaneSpectrum` will use this value for
            the atmospheric temperature. If not provided and `vane` is an `int`, `~dysh.spectra.vane.VaneSpectrum` will try to fetch a
            value from the GBO weather forecast script (only available at GBO).
        t_bkg : float or None
            Background temperature in K. If `vane` is a `~dysh.spectra.vane.VaneSpectrum` it won't be used.
            If `vane` is an `int`, then the resulting `~dysh.spectra.vane.VaneSpectrum` will use this value for
            the background temperature. If not provided, it will take a default value of 2.725 K, i.e., the CMB at 3 mm.
        t_warm : float or None
            Vane temperature in K. If `vane` is a `~dysh.spectra.vane.VaneSpectrum` it won't be used.
            If `vane` is an `int`, then the resulting `~dysh.spectra.vane.VaneSpectrum` will use this value for the vane temperature.
            If not provided and `vane` is an `int`, it will take the value found in the "TWARM" column of the SDFITS.
        **kwargs : dict
            Optional additional selection keyword arguments, typically
            given as key=value, though a dictionary works too.
            e.g., `scan=[27,30], source='NGC123', ` etc.

        Raises
        ------
        Exception
            If scans matching the selection criteria are not found.

        Returns
        -------
        scanblock : `~dysh.spectra.scan.ScanBlock`
            ScanBlock containing one or more `~dysh.spectra.scan.PSScan`.
        """
        ScanBase._check_tscale(units)
        ScanBase._check_gain_factors(ap_eff, surface_error)
        self._check_vane_and_t_sys_args(vane, t_sys)

        if vane is not None:
            vane, units, requested_units, zenith_opacity = self._vane_setup(
                vane, fdnum, ifnum, plnum, units, zenith_opacity, t_warm, t_atm, t_bkg, t_cal, apply_flags
            )

        if units.lower() != "ta" and zenith_opacity is None and ap_eff is None:
            raise ValueError("Can't scale the data without a valid zenith opacity")
        _channel = self._normalize_channel_range(channel)

        prockey = "OBSTYPE"
        procvals = {"ON": "PSWITCHON", "OFF": "PSWITCHOFF"}
        (scans, _sf) = self._common_selection(
            fdnum=fdnum,
            ifnum=ifnum,
            plnum=plnum,
            apply_flags=apply_flags,
            prockey=prockey,
            procvals=procvals,
            **kwargs,
        )
        tsys = _parse_tsys(t_sys, scans)
        _tsys = tsys
        _tcal = t_cal
        _nocal = nocal
        _bintable = kwargs.get("bintable", None)

        scanblock = ScanBlock()
        for i in range(len(self._sdf)):
            _df = select_from("FITSINDEX", i, _sf)
            if len(_df) == 0:  # If nothing was selected go to next file.
                continue
            scanlist = self._common_scan_list_selection(scans, _df, prockey=prockey, procvals=procvals, check=False)
            if len(scanlist["ON"]) == 0 or len(scanlist["OFF"]) == 0:
                logger.debug(f"scans {scans} not found, continuing")
                continue
            rows = {}
            # loop over scan pairs
            for on, off in zip(scanlist["ON"], scanlist["OFF"], strict=False):
                _ondf = select_from("SCAN", on, _df)
                _offdf = select_from("SCAN", off, _df)
                rows["ON"] = list(_ondf["ROW"])
                rows["OFF"] = list(_offdf["ROW"])
                for key in rows:
                    if len(rows[key]) == 0:
                        raise Exception(f"{key} scans not found in scan list {scans}")
                # do not pass scan list here. We need all the cal rows. They will
                # be intersected with scan rows in PSScan
                calrows = {}
                calrows["ON"] = list(select_from("CAL", "T", _ondf)["ROW"]) + list(
                    select_from("CAL", "T", _offdf)["ROW"]
                )
                calrows["OFF"] = list(select_from("CAL", "F", _ondf)["ROW"]) + list(
                    select_from("CAL", "F", _offdf)["ROW"]
                )
                if t_cal is not None:
                    _tcal = t_cal
                else:
                    _tcal = self._get_tcal(_offdf["TCAL"])
                if len(calrows["ON"]) == 0 or nocal:
                    _nocal = True
                    if tsys is None:
                        dfoncalF = select_from("CAL", "F", _ondf)
                        _tsys = dfoncalF["TSYS"].to_numpy()
                        if vane is None:
                            logger.info("Using TSYS column")
                # Use user provided system temperature.
                if tsys is not None:
                    try:
                        _tsys = tsys[on][0]
                    except KeyError:
                        _tsys = tsys[off][0]
                d = {"ON": on, "OFF": off}
                if _bintable is None:
                    _bintable = self._get_bintable(_ondf)
                g = PSScan(
                    self._sdf[i],
                    scan=d,
                    scanrows=rows,
                    calrows=calrows,
                    fdnum=fdnum,
                    ifnum=ifnum,
                    plnum=plnum,
                    bintable=_bintable,
                    calibrate=calibrate,
                    smoothref=smoothref,
                    apply_flags=apply_flags,
                    tscale=units,
                    zenith_opacity=zenith_opacity,
                    nocal=_nocal,
                    tsys=_tsys,
                    tcal=_tcal,
                    ap_eff=ap_eff,
                    surface_error=surface_error,
                    channel=_channel,
                    vane=vane,
                )
                # If calibrated with a vane change the units (so ugly >.<).
                if vane is not None:
                    self._set_scale_vane(g, requested_units, zenith_opacity)
                g.merge_commentary(self)
                scanblock.append(g)
                # Reset these variables for the next scan.
                _tsys = tsys
                _tcal = t_cal
                _nocal = nocal
                _bintable = kwargs.get("bintable", None)
        if len(scanblock) == 0:
            raise Exception("Didn't find any scans matching the input selection criteria.")
        scanblock.merge_commentary(self)
        return scanblock
        # end of getps()

    @log_call_to_result
    def getnod(
        self,
        ifnum: int,
        plnum: int,
        fdnum: None | int = None,
        calibrate: bool = True,
        smoothref: int = 1,
        apply_flags: bool = True,
        t_sys=None,
        t_cal=None,
        nocal=False,
        units="ta",
        zenith_opacity=None,
        ap_eff: float | None = None,
        surface_error: Quantity | None = None,
        channel: list | None = None,
        vane: int | VaneSpectrum | None = None,
        t_atm: float | None = None,
        t_bkg: float | None = None,
        t_warm: float | None = None,
        **kwargs,
    ):
        """
        Retrieve and calibrate nodding data.

        Parameters
        ----------
        ifnum : int
            The intermediate frequency (IF) number
        plnum : int
            The polarization number
        fdnum :  2-tuple, optional
            The feed numbers. A pair of feed numbers may be given to choose different nodding beams than were used to obtain the observations.  Default: None which means use the beams found in the data.
        calibrate : boolean, optional
            Calibrate the scans.
            The default is True.
        smooth_ref : int, optional
            Smooth the reference spectra using a boxcar kernel with a width of `smooth_ref` channels.
            The default is to not smooth the reference spectra.
        apply_flags : boolean, optional.
            If True, apply flags before calibration.
            See :meth:`apply_flags`. Default: True
        units : str, optional
            The brightness scale unit for the output scan, must be one of (case-insensitive)
                    - 'ta'   : Antenna Temperature
                    - 'ta*'  : Antenna temperature corrected to above the atmosphere
                    - 'flux' : flux density in Jansky
            If 'ta*' or 'flux' the zenith opacity must also be given. Default: 'ta'
        t_sys : float or list or list of lists or dict, optional
            System temperature. If provided, it overrides the value computed using the noise diode.
            If no noise diode is fired, and `t_sys=None`, then the column "TSYS" will be used instead.
            For example, `t_sys = np.array([[30], [50]])` would use a system temperature of 30 K for
            the first feed and 50 K for the second feed. Another example, `t_sys = {1: [[50, 60]], 2: [[45],[65]], 3: [[60],[70]]}`
            would use a system temperature of 50 K for the first feed in scan 1, 60 K for the second feed in scan 1, 45 K for the first feed in scan 2, 65 K for the second feed in scan 2, 60 K for the first feed in scan 3, and 70 K for the second feed in scan 3. If passing a dict it should contain an item for every scan.
            If `vane` is provided, then `t_sys` will be ignored and `vane` will be used to derive the system temperature.
        t_cal : None or float
            Noise diode temperature. If provided, this value is used instead of the value found in the
            TCAL column of the SDFITS file. If no value is provided, default, then the TCAL column is
            used.
        nocal : bool, optional
            Is the noise diode being fired? False means the noise diode was firing.
            By default it will figure this out by looking at the "CAL" column.
            It can be set to True to override this. Default: False
        ap_eff : float or None
            Aperture efficiency o be used when scaling data to brightness temperature of flux. The provided aperture
            efficiency must be a number between 0 and 1.  If None, `dysh` will calculate it as described in
            :meth:`~GBTGainCorrection.aperture_efficiency`. Only one of `ap_eff` or `surface_error`
            can be provided.
        surface_error : `~astropy.units.Quantity` or None
            Surface rms error, in units of length (typically microns), to be used in the Ruze formula when calculating the
            aperture efficiency.  If None, `dysh` will use the known GBT surface error model.  Only one of `ap_eff` or `surface_error`
            can be provided.
        channel: list or None
            An inclusive list of `[firstchan, lastchan]` to use in the calibration. The channel list is zero-based. If provided,
            only data channels in the inclusive range `[firstchan,lastchan]` will be used. If a reference spectrum has been given, it will also be
            trimmed to `[firstchan,lastchan]`.  System temperature calculation will use 80% of the trimmed channel range.  If channels have already been selected through
            :meth:`GBTFITSLoad.select_channel`, a ValueError will be raised.
        vane : int or `~dysh.spectra.vane.VaneSpectrum` or None
            Vane calibration scan. This will be used to derive the system temperature.
            If provided, `t_sys` will be ignored.
        t_atm : float or None
            Atmospheric temperature in K. If `vane` is a `~dysh.spectra.vane.VaneSpectrum` it won't be used.
            If `vane` is an `int`, then the resulting `~dysh.spectra.vane.VaneSpectrum` will use this value for
            the atmospheric temperature. If not provided and `vane` is an `int`, `~dysh.spectra.vane.VaneSpectrum` will try to fetch a
            value from the GBO weather forecast script (only available at GBO).
        t_bkg : float or None
            Background temperature in K. If `vane` is a `~dysh.spectra.vane.VaneSpectrum` it won't be used.
            If `vane` is an `int`, then the resulting `~dysh.spectra.vane.VaneSpectrum` will use this value for
            the background temperature. If not provided, it will take a default value of 2.725 K, i.e., the CMB at 3 mm.
        t_warm : float or None
            Vane temperature in K. If `vane` is a `~dysh.spectra.vane.VaneSpectrum` it won't be used.
            If `vane` is an `int`, then the resulting `~dysh.spectra.vane.VaneSpectrum` will use this value for the vane temperature.
            If not provided and `vane` is an `int`, it will take the value found in the "TWARM" column of the SDFITS.
        **kwargs : dict
            Optional additional selection keyword arguments, typically
            given as key=value, though a dictionary works too.
            e.g., `intnum=1, plnum=2` etc.
            For multi-beam with more than 2 beams, fdnum=[BEAM1,BEAM2] must be selected,
            unless the data have been properly tagged using PROCSCAN which BEAM1 and BEAM2 are.

        Raises
        ------
        Exception
            If scans matching the selection criteria are not found.

        Returns
        -------
        scanblock : `~dysh.spectra.scan.ScanBlock`
            ScanBlock containing one or more `~dysh.spectra.scan.NodScan`.

        """

        ScanBase._check_tscale(units)
        ScanBase._check_gain_factors(ap_eff, surface_error)
        if units.lower() != "ta" and zenith_opacity is None and ap_eff is None:
            raise ValueError("Can't scale the data without a valid zenith opacity")
        if fdnum is not None and (type(fdnum) is int or len(fdnum) != 2):
            raise TypeError(f"fdnum={fdnum} not valid, need a list with two feeds")
        _channel = self._normalize_channel_range(channel)
        self._check_vane_and_t_sys_args(vane, t_sys)

        prockey = "PROCSEQN"
        procvals = {"ON": 1, "OFF": 2}
        (scans, _sf) = self._common_selection(
            fdnum=fdnum,
            ifnum=ifnum,
            plnum=plnum,
            apply_flags=apply_flags,
            prockey=prockey,
            procvals=procvals,
            **kwargs,
        )

        tsys = _parse_tsys(t_sys, scans)
        _tsys = None
        _tcal = t_cal
        _nocal = nocal
        _bintable = kwargs.get("bintable", None)

        beam1_selected = True
        scanblock = ScanBlock()

        for i in range(len(self._sdf)):
            df0 = select_from("FITSINDEX", i, _sf)
            if len(df0) == 0:
                continue
            scanlist = self._common_scan_list_selection(scans, df0, prockey=prockey, procvals=procvals, check=False)
            if len(scanlist["ON"]) == 0 or len(scanlist["OFF"]) == 0:
                logger.debug(f"Some of scans {scans} not found, continuing")
                continue
            # Loop over scan pairs.
            for on, off in zip(scanlist["ON"], scanlist["OFF"], strict=False):
                # Each scan could use a different pair of fdnums.
                if fdnum is None:
                    _fdnum = self.get_nod_beams(scan=on)
                else:
                    _fdnum = fdnum
                logger.debug(f"getnod using fdnum={_fdnum}")
                for j, f in enumerate(_fdnum):
                    _df = select_from("FDNUM", f, df0)
                    if len(_df) == 0:  # skip IF's and beams not part of the nodding pair.
                        continue

                    if vane is not None:
                        # Each beam needs its own vane, and we might not know the fdnums before this point.
                        _vane, units, requested_units, zenith_opacity = self._vane_setup(
                            vane, f, ifnum, plnum, units, zenith_opacity, t_warm, t_atm, t_bkg, t_cal, apply_flags
                        )
                    else:
                        _vane = None

                    beam1_selected = f == _fdnum[0]
                    logger.debug(f"SCANLIST {scanlist}")
                    logger.debug(f"FEED {f} {beam1_selected} {_fdnum[0]}")
                    logger.debug(f"PROCSEQN {set(_df['PROCSEQN'])}")
                    logger.debug(f"Sending dataframe with scans {set(_df['SCAN'])}")
                    logger.debug(f"and PROC {set(_df['PROC'])}")

                    rows = {}
                    _ondf = select_from("SCAN", on, _df)
                    _offdf = select_from("SCAN", off, _df)
                    rows["ON"] = list(_ondf["ROW"])
                    rows["OFF"] = list(_offdf["ROW"])
                    for key in rows:
                        if len(rows[key]) == 0:
                            raise Exception(f"{key} scans not found in scan list {scans}")
                    if _bintable is None:
                        _bintable = self._get_bintable(_ondf)
                    # Do not pass scan list here. We need all the cal rows. They will
                    # be intersected with scan rows in NodScan.
                    calrows = {}
                    dfcalT = select_from("CAL", "T", _df)
                    dfcalF = select_from("CAL", "F", _df)
                    calrows["ON"] = list(dfcalT["ROW"])
                    calrows["OFF"] = list(dfcalF["ROW"])
                    d = {"ON": on, "OFF": off}

                    if t_cal is not None:
                        _tcal = t_cal
                    else:
                        _tcal = self._get_tcal(_offdf["TCAL"])
                    # Check if there is a noise diode.
                    if len(calrows["ON"]) == 0 or nocal:
                        _nocal = True
                        if tsys is None:
                            dfoncalF = select_from("CAL", "F", _ondf)
                            _tsys = dfoncalF["TSYS"].to_numpy()
                            if vane is None:
                                logger.info("Using TSYS column")
                    # Use user provided system temperature.
                    if tsys is not None:
                        _tsys = tsys[on][j]

                    logger.debug(f"{i, f} SCANROWS {rows}")
                    logger.debug(f"BEAM1 {beam1_selected}")
                    g = NodScan(
                        self._sdf[i],
                        scan=d,
                        beam1=beam1_selected,
                        scanrows=rows,
                        calrows=calrows,
                        fdnum=f,
                        ifnum=ifnum,
                        plnum=plnum,
                        bintable=_bintable,
                        calibrate=calibrate,
                        smoothref=smoothref,
                        apply_flags=apply_flags,
                        nocal=_nocal,
                        tsys=_tsys,
                        tcal=_tcal,
                        tscale=units,
                        zenith_opacity=zenith_opacity,
                        ap_eff=ap_eff,
                        surface_error=surface_error,
                        channel=_channel,
                        vane=_vane,
                    )
                    # If calibrated with a vane change the units (so ugly >.<).
                    if vane is not None:
                        self._set_scale_vane(g, requested_units, zenith_opacity)
                    g.merge_commentary(self)
                    scanblock.append(g)
                    # Reset these variables for the next scan.
                    _tsys = tsys
                    _tcal = t_cal
                    _nocal = nocal
                    _bintable = kwargs.get("bintable", None)
        if len(scanblock) == 0:
            raise Exception("Didn't find any unflagged scans matching the input selection criteria.")
        if len(scanblock) % 2 == 1:
            raise Exception("Odd number of scans for getnod, check that your feeds are valid")
        # note the two nods are not merged, but added to the pool as two "independant" PS scans
        scanblock.merge_commentary(self)
        return scanblock
        # end of getnod()

    @log_call_to_result
    def getfs(
        self,
        fdnum: int,
        ifnum: int,
        plnum: int,
        calibrate: bool = True,
        fold: bool = True,
        shift_method: str = "fft",
        use_sig: bool = True,
        smoothref: int = 1,
        apply_flags: bool = True,
        units: str = "ta",
        zenith_opacity: float | None = None,
        t_sys=None,
        t_cal=None,
        nocal: bool = False,
        ap_eff: float | None = None,
        surface_error: Quantity | None = None,
        channel: list | None = None,
        vane: int | VaneSpectrum | None = None,
        t_atm: float | None = None,
        t_bkg: float | None = None,
        t_warm: float | None = None,
        **kwargs,
    ):
        """
        Retrieve and calibrate frequency-switched data.

        Parameters
        ----------
        fdnum: int
            The feed number
        ifnum : int
            The intermediate frequency (IF) number
        plnum : int
            The polarization number
        calibrate : boolean, optional
            Calibrate the scans. The default is True.
        fold : boolean, optional
            Fold the sig and ref scans.  The default is True.
        shift_method : str
            Method to use when shifting the spectra for folding. One of 'fft' or 'interpolate'.
            'fft' uses a phase shift in the time domain. 'interpolate' interpolates the signal. Default: 'fft'.
        use_sig : boolean, optional
            Return the sig or ref based spectrum. This applies to both the folded
            and unfolded option.  The default is True.
            NOT IMPLEMENTED YET
        smooth_ref: int, optional
            the number of channels in the reference to boxcar smooth prior to calibration
        apply_flags : boolean, optional.  If True, apply flags before calibration.
            See :meth:`apply_flags`. Default: True
        units : str, optional
            The brightness scale unit for the output scan, must be one of (case-insensitive)
                    - 'ta'   : Antenna Temperature
                    - 'ta*'  : Antenna temperature corrected to above the atmosphere
                    - 'flux' : flux density in Jansky
            If 'ta*' or 'flux' the zenith opacity must also be given. Default: 'ta'
        zenith_opacity: float, optional
                The zenith opacity to use in calculating the scale factors for the integrations.  Default:None
        t_sys : float, optional
            System temperature. If provided, it overrides the value computed using the noise diode.
            If no noise diode is fired, and `t_sys=None`, then the column "TSYS" will be used instead.
            If `vane` is provided, then `t_sys` will be ignored and `vane` will be used to derive the system temperature.
        t_cal : None or float
            Noise diode temperature. If provided, this value is used instead of the value found in the
            TCAL column of the SDFITS file. If no value is provided, default, then the TCAL column is
            used.
        nocal : bool, optional
            Is the noise diode being fired? False means the noise diode was firing.
            By default it will figure this out by looking at the "CAL" column.
            It can be set to True to override this.
        ap_eff : float or None
            Aperture efficiency o be used when scaling data to brightness temperature of flux. The provided aperture
            efficiency must be a number between 0 and 1.  If None, `dysh` will calculate it as described in
            :meth:`~GBTGainCorrection.aperture_efficiency`. Only one of `ap_eff` or `surface_error`
            can be provided.
        surface_error : `~astropy.units.Quantity` or None
            Surface rms error, in units of length (typically microns), to be used in the Ruze formula when calculating the
            aperture efficiency.  If None, `dysh` will use the known GBT surface error model.  Only one of `ap_eff` or `surface_error`
            can be provided.
        channel: list or None
            An inclusive list of `[firstchan, lastchan]` to use in the calibration. The channel list is zero-based. If provided,
            only data channels in the inclusive range `[firstchan,lastchan]` will be used. If a reference spectrum has been given, it will also be
            trimmed to `[firstchan,lastchan]`.  System temperature calculation will use 80% of the trimmed channel range.  If channels have already been selected through
            :meth:`GBTFITSLoad.select_channel`, a ValueError will be raised.

            **Note**: With certain choices of `channel`, folding the data with `shift_method='fft'` can result in
            a numpy array broadcast exception. If this occurs, either change `shift_method` to 'interpolate' or change
            the channel range by one channel to avoid the error.
        vane : int or `~dysh.spectra.vane.VaneSpectrum` or None
            Vane scalibration scan. This will be used to derive the system temperature.
            If provided, `t_sys` will be ignored.
        t_atm : float or None
            Atmospheric temperature in K. If `vane` is a `~dysh.spectra.vane.VaneSpectrum` it won't be used.
            If `vane` is an `int`, then the resulting `~dysh.spectra.vane.VaneSpectrum` will use this value for
            the atmospheric temperature. If not provided and `vane` is an `int`, `~dysh.spectra.vane.VaneSpectrum` will try to fetch a
            value from the GBO weather forecast script (only available at GBO).
        t_bkg : float or None
            Background temperature in K. If `vane` is a `~dysh.spectra.vane.VaneSpectrum` it won't be used.
            If `vane` is an `int`, then the resulting `~dysh.spectra.vane.VaneSpectrum` will use this value for
            the background temperature. If not provided, it will take a default value of 2.725 K, i.e., the CMB at 3 mm.
        t_warm : float or None
            Vane temperature in K. If `vane` is a `~dysh.spectra.vane.VaneSpectrum` it won't be used.
            If `vane` is an `int`, then the resulting `~dysh.spectra.vane.VaneSpectrum` will use this value for the vane temperature.
            If not provided and `vane` is an `int`, it will take the value found in the "TWARM" column of the SDFITS.
        **kwargs : dict
            Optional additional selection keyword arguments, typically
            given as key=value, though a dictionary works too.
            e.g., `ifnum=1, plnum=[2,3]` etc.

        Raises
        ------
        Exception
            If no scans matching the selection criteria are found.

        Returns
        -------
        scanblock : `~dysh.spectra.scan.ScanBlock`
            ScanBlock containing one or more`~dysh.spectra.scan.FSScan`.

        """

        ScanBase._check_tscale(units)
        ScanBase._check_gain_factors(ap_eff, surface_error)
        self._check_vane_and_t_sys_args(vane, t_sys)

        if vane is not None:
            vane, units, requested_units, zenith_opacity = self._vane_setup(
                vane, fdnum, ifnum, plnum, units, zenith_opacity, t_warm, t_atm, t_bkg, t_cal, apply_flags
            )

        if units.lower() != "ta" and zenith_opacity is None and ap_eff is None:
            raise ValueError("Can't scale the data without a valid zenith opacity")
        _channel = self._normalize_channel_range(channel)
        (scans, _sf) = self._common_selection(ifnum=ifnum, plnum=plnum, fdnum=fdnum, apply_flags=apply_flags, **kwargs)

        tsys = _parse_tsys(t_sys, scans)
        _tsys = None
        _tcal = t_cal
        _nocal = nocal
        _bintable = kwargs.get("bintable", None)

        scanblock = ScanBlock()

        for i in range(len(self._sdf)):
            logger.debug(f"Processing file {i}: {self._sdf[i].filename}")

            df = select_from("FITSINDEX", i, _sf)
            if len(df) == 0:
                continue

            # loop over scans:
            for scan in scans:
                logger.debug(f"doing scan {scan}")
                calrows = {}
                _df = select_from("SCAN", scan, df)
                if len(_df) == 0:
                    continue
                if _bintable is None:
                    _bintable = self._get_bintable(_df)
                dfcalT = select_from("CAL", "T", _df)
                dfcalF = select_from("CAL", "F", _df)
                sigrows = {}
                dfsigT = select_from("SIG", "T", _df)
                dfsigF = select_from("SIG", "F", _df)

                calrows["ON"] = list(dfcalT["ROW"])
                calrows["OFF"] = list(dfcalF["ROW"])
                sigrows["ON"] = list(dfsigT["ROW"])
                sigrows["OFF"] = list(dfsigF["ROW"])

                if t_cal is not None:
                    _tcal = t_cal
                else:
                    _tcal = self._get_tcal(dfcalF["TCAL"])

                # Is there a noise diode?
                if len(calrows["ON"]) == 0 or nocal:
                    _nocal = True
                    # User did not provide a system temperature.
                    if tsys is None:
                        _tsys = dfcalF["TSYS"].to_numpy()
                        if vane is None:
                            logger.info("Using TSYS column")
                # User provided a system temperature.
                if tsys is not None:
                    _tsys = tsys[scan][0]

                g = FSScan(
                    self._sdf[i],
                    scan=scan,
                    sigrows=sigrows,
                    calrows=calrows,
                    fdnum=fdnum,
                    ifnum=ifnum,
                    plnum=plnum,
                    bintable=_bintable,
                    calibrate=calibrate,
                    fold=fold,
                    shift_method=shift_method,
                    use_sig=use_sig,
                    smoothref=smoothref,
                    apply_flags=apply_flags,
                    tscale=units,
                    zenith_opacity=zenith_opacity,
                    tsys=_tsys,
                    tcal=_tcal,
                    nocal=_nocal,
                    ap_eff=ap_eff,
                    surface_error=surface_error,
                    channel=_channel,
                    vane=vane,
                )
                # If calibrated with a vane change the units (so ugly >.<).
                if vane is not None:
                    self._set_scale_vane(g, requested_units, zenith_opacity)
                g.merge_commentary(self)
                scanblock.append(g)
                # Reset these variables for the next scan.
                _tsys = tsys
                _tcal = t_cal
                _nocal = nocal
                _bintable = kwargs.get("bintable", None)
        if len(scanblock) == 0:
            raise Exception("Didn't find any unflagged scans matching the input selection criteria.")
        scanblock.merge_commentary(self)
        return scanblock
        # end of getfs()

    def _fix_ka_rx_if_needed(self):
        """The Ka band receiver had mislabeled FDNUM for a period of time.
        This corrects FDNUM for the given polarization (the Ka beams measure orthogonal polarizations).

        **Note**:   I don't know what date range this was effective, so this
        method may do the wrong thing for some data!

        See issue #160 https://github.com/GreenBankObservatory/dysh/issues/160
        """
        # Check if we are dealing with Ka data before the beam switch.
        rx = self["FRONTEND"].unique()
        if "Rcvr26_40" not in rx:
            return

        self._fix_column("FDNUM", 1, {"FRONTEND": "Rcvr26_40", "PLNUM": 1})
        logger.info("Fixing FDNUM mislabel for Rcvr26_40. FDNUM 0 changed to 1")
        self._fix_column("FDNUM", 0, {"FRONTEND": "Rcvr26_40", "PLNUM": 0})
        logger.info("Fixing FDNUM mislabel for Rcvr26_40. FDNUM 1 changed to 0")

    @log_call_to_result
    def gettcal(
        self,
        scan: int,
        ifnum: int,
        plnum: int,
        zenith_opacity: float,
        ref: None | int | Spectrum = None,
        fdnum: None | int = None,
        apply_flags: bool = True,
        method=None,
        name=None,
        fluxscale=None,
        method_kwargs: None | dict = None,
        ap_eff=None,
        surface_error=None,
        **kwargs,
    ):
        """
        Derive the noise diode temperature from observations of a flux calibrator.

        Parameters
        ----------
        scan : int
            The scan number.
        fdnum : int
            The feed number.
        ifnum : int
            The intermediate frequency (IF) number.
        plnum : int
            The polarization number.
        ref : int or `~dysh.spectra.spectrum.Spectrum`
            The reference scan number or a `~dysh.spectra.spectrum.Spectrum` object.  If an integer is given,
            the reference spectrum will be the total power time-averaged spectrum using the weights given.
            This is only used if `method=GBTFITSLoad.getsigref`.
        apply_flags : boolean, optional
            If True, apply flags before calibration.
            See :meth:`apply_flags`. Default: True
        zenith_opacity : float
            The zenith opacity to use in calculating the scale factors for the integrations. Default: None
        method : callable
            Method to use for calibrating the data.
            It can be one of `GBTFITSLoad.getsigref`, `GBTFITSLoad.getps`, `GBTFITSLoad.getnod` or `GBTFITSLoad.subbeamnod`.
            If None, the default, it will use `GBTFITSLoad.getsigref` for Track observations,
            `GBTFITSLoad.getps` for OnOff or OffOn observations,
            `GBTFITSLoad.getnod` for Nod observations, and
            `GBTFITSLoad.subbeamnod` for SubBeamNod observations.
        name : str
            Alternative name for the calibrator source.
            This will override the value found in the "OBJECT" column of the SDFITS.
            Useful when the "OBJECT" column contains a value not present in the calibrator catalog.
        fluxscale : str
            Name of the flux scale to use to compute the flux of the calibrator.
            "Perley-Butler 2017" and "Ott 1994" are known to dysh, although the user can provide other scales.
        method_kwargs : dict
            Dictionary with additional keywords to pass to the calibration `method`.
        ap_eff : float or None
            Aperture efficiency o be used when scaling data to brightness temperature of flux. The provided aperture
            efficiency must be a number between 0 and 1.  If None, `dysh` will calculate it as described in
            :meth:`~GBTGainCorrection.aperture_efficiency`. Only one of `ap_eff` or `surface_error`
            can be provided.
        surface_error: Quantity or None
            Surface rms error, in units of length (typically microns), to be used in the Ruze formula when calculating the
            aperture efficiency.  If None, `dysh` will use the known GBT surface error model.  Only one of `ap_eff` or `surface_error`
            can be provided.
        **kwargs : dict
            Optional additional selection keyword arguments, typically
            given as key=value, though a dictionary works too.
            e.g., `object='NGC123'`.

        Returns
        -------
        tcal : `~dysh.spectra.tcal.TCal`
            Object that contains the noise diode temperature in its flux attribute.

        Raises
        ------
        TypeError
            If more than one scan is provided.
        TypeError
            If `method` is not recognized.
        """

        valid_procs = {
            "Track": self.getsigref,
            "OnOff": self.getps,
            "OffOn": self.getps,
            "Nod": self.getnod,
            "SubBeamNod": self.subbeamnod,
        }

        if method_kwargs is None:
            method_kwargs = {}

        if not isinstance(scan, int):
            raise TypeError(f"Only a single integer value allowed for `scan`. Got {scan}")

        if method is not None:
            if method not in valid_procs.values():
                valid_methods = [m.__qualname__ for m in valid_procs.values()]
                raise TypeError(f"Unrecognized method ({method}). It should be one of {valid_methods}")

        (scans, _sf) = self._common_selection(
            fdnum=fdnum,
            ifnum=ifnum,
            plnum=plnum,
            apply_flags=apply_flags,
            scan=scan,
            **kwargs,
        )

        if name is None:
            name = _sf["OBJECT"].unique()[0]
        target = Calibrator.from_name(name, scale=fluxscale)

        proc = _sf["PROC"].unique()[0]

        if method is None:
            method = valid_procs[proc]
        logger.info(f"Will use {method.__name__} to calibrate the data.")

        method_args = {
            "scan": scans,
            "fdnum": fdnum,
            "ifnum": ifnum,
            "plnum": plnum,
            "apply_flags": apply_flags,
            "zenith_opacity": zenith_opacity,
            "t_cal": 1.0,
            "units": "flux",
            "ap_eff": ap_eff,
            "surface_error": surface_error,
        }
        if ref is not None:
            method_args["ref"] = ref

        # Merge kwargs.
        # Merging in this order implies that kwargs will supersede method_kwargs.
        method_kwargs.update(kwargs)

        obs_ta = method(**method_args, **method_kwargs).timeaverage()

        nu = obs_ta.spectral_axis
        snu = target.compute_sed(nu.quantity)

        tcal_values = (snu / obs_ta.flux).value * u.K

        tcal = TCal.from_spectrum(obs_ta, data=tcal_values, snu=snu, name=name)

        return tcal

    @log_call_to_result
    def subbeamnod(
        self,
        fdnum: int,
        ifnum: int,
        plnum: int,
        method="cycle",
        calibrate: bool = True,
        weights="tsys",
        smoothref: int = 1,
        apply_flags: bool = True,
        units="ta",
        zenith_opacity=None,
        t_sys=None,
        t_cal=None,
        ap_eff: float | None = None,
        surface_error: Quantity | None = None,
        channel: list | None = None,
        nocal: bool = False,
        vane: int | VaneSpectrum | None = None,
        t_atm: float | None = None,
        t_bkg: float | None = None,
        t_warm: float | None = None,
        **kwargs,
    ):
        """Calibrate a SubBeamNod scan.

        Parameters
        ----------
        fdnum : int
            The feed number.
        ifnum : int
            The intermediate frequency (IF) number.
        plnum : int
            The polarization number.
        method : str
            Method to use when processing. One of 'cycle' or 'scan'.
            'cycle' (default) treats each SUBREF_STATE independently, resulting in multiple signal and reference states per scan..
            'scan' averages the SUBREF_STATE rows resulting in one signal and reference state per scan.
        calibrate : bool
            Whether or not to calibrate the data.
        weights : str or None
            Weights to use for the time averaging of the sub reflector states.
            None to indicate equal weighting or 'tsys' to indicate inverse variance weights.
        smooth_ref : int, optional
            The boxcar kernel width to smooth the reference spectra prior to calibration.
        apply_flags : boolean, optional.
            If True, apply flags before calibration.
            See :meth:`apply_flags`.
        units : str, optional
            The brightness scale unit for the output scan, must be one of (case-insensitive)
                    - 'ta'   : Antenna Temperature
                    - 'ta*'  : Antenna temperature corrected to above the atmosphere
                    - 'flux' : flux density in Jansky
            If 'ta*' or 'flux' the zenith opacity must also be given. Default: 'ta'
        zenith_opacity : float, optional
            The zenith opacity to use to correct the data for atmospheric opacity.
        t_sys : float, optional
            System temperature. If provided, it overrides the value computed using the noise diode.
            If no noise diode is fired, and `t_sys=None`, then the column "TSYS" will be used instead.
            If `vane` is provided, then `t_sys` will be ignored and `vane` will be used to derive the system temperature.
        t_cal : None or float
            Noise diode temperature. If provided, this value is used instead of the value found in the
            TCAL column of the SDFITS file. If no value is provided, default, then the TCAL column is
            used.
        channel: list or None
            An inclusive list of `[firstchan, lastchan]` to use in the calibration. The channel list is zero-based. If provided,
            only data channels in the inclusive range `[firstchan,lastchan]` will be used. If a reference spectrum has been given, it will also be
            trimmed to `[firstchan,lastchan]`.  System temperature calculation will use 80% of the trimmed channel range.  If channels have already been selected through
            :meth:`GBTFITSLoad.select_channel`, a ValueError will be raised.
        ap_eff : float or None
            Aperture efficiency o be used when scaling data to brightness temperature of flux. The provided aperture
            efficiency must be a number between 0 and 1.  If None, `dysh` will calculate it as described in
            :meth:`~GBTGainCorrection.aperture_efficiency`. Only one of `ap_eff` or `surface_error`
            can be provided.
        surface_error : `~astropy.units.Quantity` or None
            Surface rms error, in units of length (typically microns), to be used in the Ruze formula when calculating the
            aperture efficiency.  If None, `dysh` will use the known GBT surface error model.  Only one of `ap_eff` or `surface_error`
            can be provided.
        vane : int or `~dysh.spectra.vane.VaneSpectrum` or None
            Vane scalibration scan. This will be used to derive the system temperature.
            If provided, `t_sys` will be ignored.
        t_atm : float or None
            Atmospheric temperature in K. If `vane` is a `~dysh.spectra.vane.VaneSpectrum` it won't be used.
            If `vane` is an `int`, then the resulting `~dysh.spectra.vane.VaneSpectrum` will use this value for
            the atmospheric temperature. If not provided and `vane` is an `int`, `~dysh.spectra.vane.VaneSpectrum` will try to fetch a
            value from the GBO weather forecast script (only available at GBO).
        t_bkg : float or None
            Background temperature in K. If `vane` is a `~dysh.spectra.vane.VaneSpectrum` it won't be used.
            If `vane` is an `int`, then the resulting `~dysh.spectra.vane.VaneSpectrum` will use this value for
            the background temperature. If not provided, it will take a default value of 2.725 K, i.e., the CMB at 3 mm.
        t_warm : float or None
            Vane temperature in K. If `vane` is a `~dysh.spectra.vane.VaneSpectrum` it won't be used.
            If `vane` is an `int`, then the resulting `~dysh.spectra.vane.VaneSpectrum` will use this value for the vane temperature.
            If not provided and `vane` is an `int`, it will take the value found in the "TWARM" column of the SDFITS.
        **kwargs : dict
            Optional additional selection keyword arguments, typically
            given as key=value, though a dictionary works too.
            e.g., `ifnum=1, plnum=[2,3]` etc.

        Returns
        -------
        data : `~dysh.spectra.scan.ScanBlock`
            A ScanBlock containing one or more `~dysh.spectra.scan.SubBeamNodScan`
        """

        ScanBase._check_tscale(units)
        ScanBase._check_gain_factors(ap_eff, surface_error)
        self._check_vane_and_t_sys_args(vane, t_sys)

        if vane is not None:
            vane, units, requested_units, zenith_opacity = self._vane_setup(
                vane, fdnum, ifnum, plnum, units, zenith_opacity, t_warm, t_atm, t_bkg, t_cal, apply_flags
            )

        if units.lower() != "ta" and zenith_opacity is None and ap_eff is None:
            raise ValueError("Can't scale the data without a valid zenith opacity")
        _channel = self._normalize_channel_range(channel)
        _bintable = kwargs.get("bintable", None)

        (scans, _sf) = self._common_selection(ifnum=ifnum, plnum=plnum, fdnum=fdnum, apply_flags=apply_flags, **kwargs)

        tsys = _parse_tsys(t_sys, scans)
        _tcal = t_cal

        scanblock = ScanBlock()

        if method == "cycle":
            # Calibrate each cycle individually and then
            # average the calibrated data.
            for sdfi in range(len(self._sdf)):
                _df = select_from("FITSINDEX", sdfi, _sf)
                for scan in scans:
                    reftp = []
                    sigtp = []
                    logger.debug(f"doing scan {scan}")
                    df = select_from("SCAN", scan, _df)
                    if len(df) == 0:
                        continue
                    if _bintable is None:
                        _bintable = self._get_bintable(df)
                    df_on = df[df["CAL"] == "T"]
                    df_off = df[df["CAL"] == "F"]
                    df_on_sig = df_on[df_on["SUBREF_STATE"] == -1]
                    df_on_ref = df_on[df_on["SUBREF_STATE"] == 1]
                    df_off_sig = df_off[df_off["SUBREF_STATE"] == -1]
                    df_off_ref = df_off[df_off["SUBREF_STATE"] == 1]
                    logger.debug(f"SCANs in df_on_sig {set(df_on_sig['SCAN'])}")
                    logger.debug(f"SCANs in df_on_ref {set(df_on_ref['SCAN'])}")
                    logger.debug(f"SCANs in df_off_sig {set(df_off_sig['SCAN'])}")
                    logger.debug(f"SCANs in df_off_ref {set(df_off_ref['SCAN'])}")
                    sig_on_rows = df_on_sig["ROW"].to_numpy()
                    ref_on_rows = df_on_ref["ROW"].to_numpy()
                    sig_off_rows = df_off_sig["ROW"].to_numpy()
                    ref_off_rows = df_off_ref["ROW"].to_numpy()

                    # Define how large of a gap between rows we will tolerate to consider
                    # a row as part of a cycle.
                    # Thinking about it, we should use the SUBREF_STATE=0 as delimiter rather
                    # than this.
                    # stepsize = len(ifnum) * len(plnum) * 2 + 1
                    stepsize = len(self.udata("IFNUM", 0)) * len(self.udata("PLNUM", 0)) * 2 + 1
                    ref_on_groups = consecutive(ref_on_rows, stepsize=stepsize)
                    sig_on_groups = consecutive(sig_on_rows, stepsize=stepsize)
                    ref_off_groups = consecutive(ref_off_rows, stepsize=stepsize)
                    sig_off_groups = consecutive(sig_off_rows, stepsize=stepsize)
                    # Make sure we have enough signal and reference pairs.
                    # Same number of cycles or less signal cycles.
                    if len(sig_on_groups) <= len(ref_on_groups):
                        pairs = {i: i for i in range(len(sig_on_groups))}
                    # One more signal cycle. Re-use one reference cycle.
                    elif len(sig_on_groups) - 1 == len(ref_on_groups):
                        pairs = {i: i for i in range(len(sig_on_groups))}
                        pairs[len(sig_on_groups) - 1] = len(ref_on_groups) - 1
                    else:
                        e = f"""There are {len(sig_on_groups)} and {len(ref_on_groups)} signal and reference cycles.
                                Try using method='scan'."""
                        raise ValueError(e)

                    if t_cal is not None:
                        _tcal = t_cal
                    else:
                        _tcal = self._get_tcal(df_off["TCAL"])

                    # Loop over cycles, calibrating each independently.
                    groups_zip = zip(ref_on_groups, sig_on_groups, ref_off_groups, sig_off_groups, strict=False)
                    for rgon, sgon, rgoff, sgoff in groups_zip:
                        # Do it the dysh way.
                        # TODO: use gettp instead of TPScan.
                        calrows = {"ON": rgon, "OFF": rgoff}
                        tprows = np.sort(np.hstack((rgon, rgoff)))
                        _reftp = TPScan(
                            self._sdf[sdfi],
                            scan,
                            None,
                            None,
                            tprows,
                            calrows,
                            fdnum=fdnum,
                            ifnum=ifnum,
                            plnum=plnum,
                            bintable=_bintable,
                            calibrate=calibrate,
                            apply_flags=apply_flags,
                            tcal=_tcal,
                            channel=_channel,
                        )
                        if tsys is not None:
                            _reftp._tsys[:] = tsys[scan][0]
                        reftp.append(_reftp)
                        calrows = {"ON": sgon, "OFF": sgoff}
                        tprows = np.sort(np.hstack((sgon, sgoff)))
                        sigtp.append(
                            TPScan(
                                self._sdf[sdfi],
                                scan,
                                None,
                                None,
                                tprows,
                                calrows,
                                fdnum=fdnum,
                                ifnum=ifnum,
                                plnum=plnum,
                                bintable=_bintable,
                                calibrate=calibrate,
                                apply_flags=apply_flags,
                                tcal=_tcal,
                                channel=_channel,
                            )
                        )
                    sb = SubBeamNodScan(
                        sigtp,
                        reftp,
                        fdnum=fdnum,
                        ifnum=ifnum,
                        plnum=plnum,
                        calibrate=calibrate,
                        smoothref=smoothref,
                        apply_flags=apply_flags,
                        tscale=units,
                        zenith_opacity=zenith_opacity,
                        weights=weights,
                        tcal=_tcal,
                        ap_eff=ap_eff,
                        surface_error=surface_error,
                        vane=vane,
                    )
                    # If calibrated with a vane change the units (so ugly >.<).
                    if vane is not None:
                        self._set_scale_vane(sb, requested_units, zenith_opacity)
                    sb.merge_commentary(self)
                    scanblock.append(sb)
                    _bintable = kwargs.get("bintable", None)
        elif method == "scan":
            # Process the whole scan as a single block.
            # This allows calibrating the data if the scan
            # was aborted and there are not enough sig/ref
            # cycles to do a per cycle calibration.
            for scan in scans:
                kwargs["scan"] = scan
                reftp = []
                sigtp = []
                tpon = self.gettp(
                    fdnum=fdnum,
                    ifnum=ifnum,
                    plnum=plnum,
                    sig=None,
                    cal=None,
                    subref=-1,
                    calibrate=calibrate,
                    apply_flags=apply_flags,
                    t_cal=t_cal,
                    channel=_channel,
                    **kwargs,
                )
                sigtp.append(tpon[0])
                tpoff = self.gettp(
                    fdnum=fdnum,
                    ifnum=ifnum,
                    plnum=plnum,
                    sig=None,
                    cal=None,
                    subref=1,
                    calibrate=calibrate,
                    apply_flags=apply_flags,
                    t_sys=t_sys,
                    t_cal=t_cal,
                    channel=_channel,
                    **kwargs,
                )
                reftp.append(tpoff[0])
                sb = SubBeamNodScan(
                    sigtp,
                    reftp,
                    fdnum=fdnum,
                    ifnum=ifnum,
                    plnum=plnum,
                    calibrate=calibrate,
                    smoothref=smoothref,
                    apply_flags=apply_flags,
                    tscale=units,
                    zenith_opacity=zenith_opacity,
                    weights=weights,
                    ap_eff=ap_eff,
                    surface_error=surface_error,
                    tcal=tpoff[0].getspec(0).meta["TCAL"],
                    vane=vane,
                )
                # If calibrated with a vane change the units (so ugly >.<).
                if vane is not None:
                    self._set_scale_vane(sb, requested_units, zenith_opacity)
                sb.merge_commentary(self)
                scanblock.append(sb)
        if len(scanblock) == 0:
            raise Exception("Didn't find any unflagged scans matching the input selection criteria.")
        scanblock.merge_commentary(self)
        return scanblock

    def _common_scan_list_selection(self, scans, selection, prockey, procvals, check=False):
        # First list item is the PROCSEQN that defines the ON.
        # Second list item is the delta between ON and OFF, delta=OFF-ON.
        proc_dict = {
            "OnOff": [1, 1],
            "OffOn": [2, -1],
            "Nod": [1, 1],
        }
        scan_selection = {"ON": [], "OFF": []}
        df = selection[selection["SCAN"].isin(scans)]
        procset = df["PROC"].unique()
        lenprocset = len(procset)
        if lenprocset == 0:
            # This is ok since not all files in a set have all the polarizations, feeds, or IFs
            return scan_selection

        for proc in procset:
            # This method should only be used by these observing procedures.
            if proc not in proc_dict.keys():
                continue

            _proc_on = df.loc[(df["PROC"] == proc) & (df["PROCSEQN"] == proc_dict[proc][0])]["SCAN"]
            _proc_off = df.loc[(df["PROC"] == proc) & (df["PROCSEQN"] == sum(proc_dict[proc]))]["SCAN"]
            proc_on = list(set(_proc_on)) + list(set(_proc_off - proc_dict[proc][1]))
            proc_off = list(set(_proc_off)) + list(set(_proc_on + proc_dict[proc][1]))

            # Check that no bogus scan numbers were added.
            df_on = selection[selection["SCAN"].isin(proc_on)]
            df_off = selection[selection["SCAN"].isin(proc_off)]
            bogus_on = df_on["PROCSEQN"] != proc_dict[proc][0]
            bogus_off = df_off["PROCSEQN"] != sum(proc_dict[proc])
            for s in set(df_on[bogus_on]["SCAN"]):
                proc_on.remove(s)
            for s in set(df_off[bogus_off]["SCAN"]):
                proc_off.remove(s)

            # Remove any scans with no pair.
            if len(proc_on) != len(proc_off):
                if len(proc_on) < len(proc_off):
                    for o in proc_off:
                        if o - proc_dict[proc][1] not in proc_on:
                            logger.warning(f"Scan {o} has no matching ON scan. Will not calibrate.")
                            proc_off.remove(o)
                else:
                    for o in proc_on:
                        if o + proc_dict[proc][1] not in proc_off:
                            logger.warning(f"Scan {o} has no matching OFF scan. Will not calibrate.")
                            proc_on.remove(o)

            # Add the remaining to the list of scans.
            scan_selection["ON"] += sorted(list(set(proc_on)))
            scan_selection["OFF"] += sorted(list(set(proc_off)))

        # Make sure the elements are unique.
        scan_selection["ON"] = sorted(list(set(scan_selection["ON"])))
        scan_selection["OFF"] = sorted(list(set(scan_selection["OFF"])))

        # Check again that they have the same number of elements.
        if len(scan_selection["ON"]) != len(scan_selection["OFF"]):
            raise Exception(
                f"ON and OFF scan list lengths differ {len(scan_selection['ON'])} != {len(scan_selection['OFF'])}"
            )

        return scan_selection

    def write(
        self,
        fileobj,
        multifile=True,
        flags=True,
        verbose=False,
        output_verify="exception",
        overwrite=False,
        checksum=False,
        **kwargs,
    ):
        """
        Write all or a subset of the `GBTFITSLoad` data to a new SDFITS file(s).

        Parameters
        ----------
        fileobj : str, file-like or `pathlib.Path`
            File to write to.  If a file object, must be opened in a
            writeable mode.
        multifile: bool, optional
            If True, write to multiple files if and only if there are multiple SDFITS files in this GBTFITSLoad.
            Otherwise, write to a single SDFITS file.
        flags: bool, optional
            If True, write the applied flags to a `FLAGS` column in the binary table.
        verbose: bool, optional
            If True, print out some information about number of rows written per file
        output_verify : str
            Output verification option.  Must be one of ``"fix"``,
            ``"silentfix"``, ``"ignore"``, ``"warn"``, or
            ``"exception"``.  May also be any combination of ``"fix"`` or
            ``"silentfix"`` with ``"+ignore"``, ``+warn``, or ``+exception"
            (e.g. ``"fix+warn"``).  See https://docs.astropy.org/en/latest/io/fits/api/verification.html for more info
        overwrite : bool, optional
            If ``True``, overwrite the output file if it exists. Raises an
            ``OSError`` if ``False`` and the output file exists. Default is
            ``False``.
        checksum : bool
            When `True` adds both ``DATASUM`` and ``CHECKSUM`` cards
            to the headers of all HDU's written to the file.
        **kwargs : dict
            Optional additional selection keyword arguments, typically
            given as key=value, though a dictionary works too.
            e.g., `ifnum=1, plnum=[2,3]` etc.
        """
        # debug = kwargs.pop("debug", False)
        logger.debug(kwargs)
        selection = Selection(self._index)
        if len(kwargs) > 0:
            selection._select_from_mixed_kwargs(**kwargs)
            logger.debug(selection.show())
            _final = selection.final
        else:
            _final = selection
        if len(_final) == 0:
            raise Exception("Your selection resulted in no rows to be written")
        fi = _final["FITSINDEX"].unique()
        logger.debug(f"fitsindex {fi} ")
        total_rows_written = 0
        if multifile:
            count = 0
            for k in fi:
                this_rows_written = 0
                # copy the primary HDU
                hdu = self._sdf[k]._hdu[0].copy()
                outhdu = fits.HDUList(hdu)
                # get the bintables rows as new bintables.
                df = select_from("FITSINDEX", k, _final)
                bintables = df.BINTABLE.unique()
                for b in bintables:  # loop over the bintables in this fitsfile
                    rows = df.ROW[df.BINTABLE == b].unique()
                    rows.sort()
                    lr = len(rows)
                    if lr > 0:
                        if flags:  # update the flags before we select rows
                            flagval = self._sdf[k]._flagmask[b].astype(np.uint8)
                            dim1 = np.shape(flagval)[1]
                            form = f"{dim1}B"
                            c = fits.Column(name="FLAGS", format=form, array=flagval)
                            self._sdf[k]._update_column({"FLAGS": c}, b)
                        ob = self._sdf[k]._bintable_from_rows(rows, b)
                        if len(ob.data) > 0:
                            outhdu.append(ob)
                        total_rows_written += lr
                        this_rows_written += lr
                if len(fi) > 1:
                    p = Path(fileobj)
                    # Note this will not preserve "A","B" etc suffixes in original FITS files.
                    outfile = p.parent / (p.stem + str(count) + p.suffix)
                    count += 1
                else:
                    outfile = fileobj
                # add comment and history cards to the primary HDU if applicable.
                # All files get all cards.
                for h in self.history:
                    outhdu[0].header["HISTORY"] = h
                for c in self.comments:
                    outhdu[0].header["COMMENT"] = c
                if verbose:
                    logger.info(f"Writing {this_rows_written} rows to {outfile}.")
                outhdu.writeto(outfile, output_verify=output_verify, overwrite=overwrite, checksum=checksum)
            if verbose:
                logger.info(f"Total of {total_rows_written} rows written to files.")
        else:
            hdu = self._sdf[fi[0]]._hdu[0].copy()
            outhdu = fits.HDUList(hdu)
            for k in fi:
                df = select_from("FITSINDEX", k, _final)
                bintables = df.BINTABLE.unique()
                for b in bintables:
                    rows = df.ROW[df.BINTABLE == b].unique()
                    rows.sort()
                    lr = len(rows)
                    if lr > 0:
                        if flags:  # update the flags before we select rows
                            flagval = self._sdf[k]._flagmask[b].astype(np.uint8)
                            dim1 = np.shape(flagval)[1]
                            form = f"{dim1}B"
                            # tdim = f"({dim1}, 1, 1, 1)" # let fitsio do this
                            c = fits.Column(name="FLAGS", format=form, array=flagval)
                            self._sdf[k]._update_column({"FLAGS": c}, b)
                        ob = self._sdf[k]._bintable_from_rows(rows, b)
                        if len(ob.data) > 0:
                            outhdu.append(ob)
                        total_rows_written += lr
            # add history and comment cards to primary header if applicable
            for h in self.history:
                outhdu[0].header["HISTORY"] = h
            for c in self.comments:
                outhdu[0].header["COMMENT"] = c
            if total_rows_written == 0:  # shouldn't happen, caught earlier
                raise Exception("Your selection resulted in no rows to be written")
            elif verbose:
                logger.info(f"Writing {total_rows_written} to {fileobj}")
            # outhdu.update_extend()  # possibly unneeded
            outhdu.writeto(fileobj, output_verify=output_verify, overwrite=overwrite, checksum=checksum)
            outhdu.close()

    def _update_radesys(self):
        """
        Updates the 'RADESYS' column of the index for cases when it is empty.
        """

        radesys = {"AzEl": "AltAz", "HADec": "hadec", "Galactic": "galactic"}

        warning_msg = (  # noqa: E731
            lambda scans,
            a,
            coord,
            limit: f"""Scan(s) {scans} have {a} {coord} below {limit}. The GBT does not go that low. Any operations that rely on the sky coordinates are likely to be inaccurate (e.g., switching velocity frames)."""
        )

        # Elevation below the GBT elevation limit (5 degrees) warning.
        low_el_mask = self["ELEVATIO"] < 5
        if low_el_mask.sum() > 0:
            low_el_scans = map(str, set(self._index.loc[low_el_mask, "SCAN"]))
            warnings.warn(warning_msg(",".join(low_el_scans), "an", "elevation", "5 degrees"))  # noqa: B028

        # Azimuth and elevation case.
        self._fix_column("RADESYS", radesys["AzEl"], {"CTYPE2": "AZ", "CTYPE3": "EL"})

        # Hour angle and declination case.
        self._fix_column("RADESYS", radesys["HADec"], {"CTYPE2": "HA"})

        # Galactic coordinates.
        self._fix_column("RADESYS", radesys["Galactic"], {"CTYPE2": "GLON"})

    def _fix_column(self, column, new_val, mask_dict):
        """
        Update the values of an existing SDFITS `column` with `new_val` where `mask_dict` is true.
        This updates `GBTFITSLoad._index` and `GBTFITSLoad._sdf.index`.
        This is mainly used to "fix" values.

        Parameters
        ----------
        column : str
            SDFITS column to update.
        new_val : str or float
            New value for `column`.
        mask_dict : dict
            Dictionary with column names and column values as keys and values.
            This will be used to determine where `GBTFITSLoad[key] == value`.
            Multiple keys and values will be combined using `numpy.logical_and`.
        """
        _mask = self._column_mask(mask_dict)
        if _mask.sum() == 0:
            return
        # Update self._index.
        self._index.loc[_mask, column] = new_val
        # Update SDFITSLoad.index.
        sdf_idx = set(self["FITSINDEX"][_mask])
        for i in sdf_idx:
            sdfi = self._sdf[i].index()
            _mask = self._sdf[i]._column_mask(mask_dict)
            sdfi.loc[_mask, column] = new_val

    def __getitem__(self, items):
        # items can be a single string or a list of strings.
        # Want case insensitivity
        # @todo deal with "DATA"
        if isinstance(items, str):
            items = items.upper()
        elif isinstance(items, (Sequence, np.ndarray)):
            items = [i.upper() for i in items]
        else:
            raise KeyError(f"Invalid key {items}. Keys must be str or list of str")
        if "DATA" in items:
            return np.vstack([s["DATA"] for s in self._sdf])
        return self._selection[items]

    @log_call_to_history
    def __setitem__(self, items, values):
        # @todo deal with "DATA"
        if isinstance(items, str):
            items = items.upper()
        # we won't support multiple keys for setting right now.
        # ultimately it could be done with recursive call to __setitem__
        # for each key/val pair
        # elif isinstance(items, (Sequence, np.ndarray)):
        #    items = [i.upper() for i in items]
        else:
            raise KeyError(f"Invalid key {items}. Keys must be str")
        if isinstance(items, str):
            iset = set([items])
        else:
            iset = set(items)
        col_exists = len(set(self.columns).intersection(iset)) > 0
        # col_in_selection =
        if col_exists:
            warnings.warn(f"Changing an existing SDFITS column {items}")  # noqa: B028
        # now deal with values as arrays
        is_array = False
        if isinstance(values, (Sequence, np.ndarray)) and not isinstance(values, str):
            if len(values) != self.total_rows:
                raise ValueError(
                    f"Length of values array ({len(values)}) for column {items} and total number of rows"
                    f" ({self.total_rows}) aren't equal."
                )
            is_array = True
        if "DATA" not in items:  # DATA is not a column in the selection
            self._selection[items] = values
        start = 0
        # loop over the individual files
        for s in self._sdf:
            if not is_array:
                s[items] = values
            else:
                s[items] = values[start : start + s.total_rows]
                start = start + s.total_rows
        selected_cols = self.selection.columns_selected()
        if items in selected_cols:
            warnings.warn(  # noqa: B028
                f"You have changed the metadata for a column that was previously used in a data selection [{items}]."
                " You may wish to update the selection. "
            )

    @log_call_to_history
    def qd_correct(self, ignore_jump: bool = False) -> None:
        """
        Apply quadrant detector (QD) corrections to sky coordinates.
        During an observation the QD records the motion of the GBT feed arm in
        the elevation and cross elevation directions. This movement results in
        pointing errors which are not automatically corrected for. This method
        allows users to correct for this movement, when possible. Typically,
        these corrections are of the order of a few arcseconds.
        More details about the QD and its use can be found in this reference:
        `<https://ui.adsabs.harvard.edu/abs/2011PASP..123..682R/abstract>`_

        Parameters
        ----------
        ignore_jump : bool
            Whether to ignore the prescence of jumps in the QD data.
            If set to `True` the corrections will be applied even if jumps are found.
            These jumps are also referred to as hysteresis events.
        """

        if self._qd_corrected:
            logger.debug("GBTFITSLoad already corrected.")
            return

        # Check for jumps.
        if self._check_qd_jump() and not ignore_jump:
            logger.warning(
                "There is a jump in the quadrant detector data. We do not recommend applying these corrections in this case. If you wish to proceed call with `ignore_jump=True`."
            )
            return

        # Check that there is only one frame type.
        # This assumes that self._update_radesys has already been called
        # so that "RADESYS" is populated.
        frame = self["RADESYS"].apply(str.lower).to_numpy().astype(str)
        if len(set(frame)) > 1:
            raise TypeError("Only a single coordinate system per observation is supported for now.")
        frame = frame[0]

        lon = self["CRVAL2"].to_numpy()
        lat = self["CRVAL3"].to_numpy()
        time = self["DATE-OBS"].to_numpy().astype(str)
        location = self.GBT

        # Transform to AzEl.
        # Only those that are not already in that frame.
        altaz_mask = frame != "AltAz"
        altaz = eq2hor(lon[altaz_mask], lat[altaz_mask], frame, time[altaz_mask], location=location)

        # Update values.
        az = copy.deepcopy(lon)
        el = copy.deepcopy(lat)
        az[altaz_mask] = altaz.az.deg
        el[altaz_mask] = altaz.alt.deg
        elfac = np.cos(np.deg2rad(el))

        # Only apply if the QD data has not been flagged.
        qd_good = self["QD_BAD"] == 0

        if qd_good.sum() == 0:
            logger.info("All quadrant detector data has been flagged. Will not apply corrections.")
            return

        # Apply corrections.
        az[qd_good] -= (self["QD_XEL"] / elfac)[qd_good]
        el[qd_good] -= self["QD_EL"][qd_good]

        # Convert back to sky coordinates.
        lonlat = hor2eq(az[altaz_mask], el[altaz_mask], frame, time[altaz_mask], location=location)
        lon[altaz_mask] = lonlat.data.lon.to("deg").value
        lat[altaz_mask] = lonlat.data.lat.to("deg").value

        # Update the sky coordinates.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self["CRVAL2"] = lon
            self["CRVAL3"] = lat

        self._qd_corrected = True

    def _check_qd_jump(self):
        """
        Check the quadrant detector data for jumps.
        A jump is defined as a change by more than `threshold`
        between consecutive samples. For now the `threshold`
        is hardcoded to 10".
        """

        threshold = 10  # arcsec.
        # TODO: set this as a configurable parameter.
        return np.any(np.diff(self["QD_XEL"]) * 3600 > threshold) or np.any(np.diff(self["QD_EL"]) * 3600 > threshold)

    def _qd_mask(self, threshold: float = 0.5) -> bool:
        """
        Generate a mask of rows where the pointing errors registered by the quadrant detector (QD)
        are larger than `threshold` times the half power beam width.

        Parameters
        ----------
        threshold : float
            Fraction of the HPBW to mask as bad.
        """

        diam = 100 * u.m
        lmbd = (self["CRVAL1"].to_numpy() * u.Hz).to("m", equivalencies=u.spectral())
        hpbw = (1.2 * lmbd / diam * u.rad).to("deg").value
        qdtot = np.sqrt(self["QD_XEL"] ** 2 + self["QD_EL"] ** 2).to_numpy()
        qdbad = self["QD_BAD"].to_numpy() == 1

        mask = abs(qdtot) > threshold * hpbw

        mask[qdbad] = False

        return mask

    def qd_check(self, threshold: float = 0.5) -> None:
        """
        Check that the pointing errors registered by the quadrant detector (QD) are less than `threshold` times the beam width.

        Parameters
        ----------
        threshold : float
            Fraction of the beam width used to check.
        """

        mask = self._qd_mask(threshold)
        bad_frac = mask.sum() / len(mask) * 100
        logger.info(f"{bad_frac}% of the data has a pointing error of more than {threshold} times HPBW.")

    def qd_flag(self, threshold: float = 0.5) -> None:
        """
        Flag rows where the pointing errors registered by the quadrant detector (QD)
        are more than `threshold` times the half power beam width.

        Parameters
        ----------
        threshold : float
            Fraction of the beam width used to check.
        """

        mask = self._qd_mask(threshold)
        if mask.sum() == 0.0:
            # Nothing to flag.
            return
        flag_rows = np.where(mask == True)[0].tolist()  # noqa: E712
        self.flag(row=flag_rows)

    def _get_beam(self, scan, mask, bi=1):
        if mask.sum() == 0:
            raise ValueError(f"Scan {scan} does not have column 'PROCSCAN' with 'BEAM{bi}' values.")
        else:
            feed = set(self["FDNUM"][mask])
            if len(feed) > 1:
                raise ValueError(f"Scan {scan} contains more than one FDNUM for 'PROCSCAN'='BEAM{bi}'.")
            feed = next(iter(feed))
        return feed

    def get_nod_beams(self, scan):
        """
        Find the FDNUM values for two nodding beams.

        Parameters
        ----------
        scan : int
            Scan for which to find the nodding beams.

        Returns
        -------
        beams : list of two ints
            Feed numbers representing the nodding beams.
            The first item is the first nodding beam.

        Raises
        ------
        TypeError
            If `scan` is not an integer.
        ValueError
            If there is no 'SCAN'=`scan` in the index, or if it is not possible to determine the nodding beams.
        """
        if not isinstance(scan, numbers.Integral):
            raise TypeError("scan must be an integer.")
        scan = [scan]

        # Find the Nod scan pairs.
        scanlist = self._common_scan_list_selection(
            scan,
            self._index,
            prockey="PROCSEQN",
            procvals={"ON": 1, "OFF": 2},
            check=False,
        )
        scans = list(itertools.chain.from_iterable(list(scanlist.values())))

        # Make the index smaller, checkling at every step that the rules contain valid data.
        mask = self["SCAN"].isin(scans)
        if mask.sum() == 0:
            raise ValueError(f"Scan {scan} not found in index.")
        mask = mask & np.isclose(self["FEEDXOFF"], 0.0) & np.isclose(self["FEEDEOFF"], 0.0)
        if mask.sum() == 0:
            raise ValueError(f"Scan {scan} does not have a beam centered on the target.")
        if "BEAM1" in set(self["PROCSCAN"][mask]) and "BEAM2" in set(self["PROCSCAN"][mask]):
            mask1 = mask & (self["PROCSCAN"] == "BEAM1")
            mask2 = mask & (self["PROCSCAN"] == "BEAM2")
            feed1 = self._get_beam(scan, mask1, bi=1)
            feed2 = self._get_beam(scan, mask2, bi=2)
        elif len(set(self["FDNUM"][mask])) == 2:
            feeds = iter(set(self["FDNUM"][mask]))
            feed1 = next(feeds)
            feed2 = next(feeds)
        else:
            raise ValueError("Cannot determine nodding beams. Please set fdnum manually.")

        return [feed1, feed2]

    def calseq(
        self,
        scan,
        fdnum=0,
        ifnum=0,
        plnum=0,
        freq: u.Quantity = None,
        tcold: float | None = None,
        twarm: float | None = None,
        apply_flags: bool = True,
    ):
        """
        This routine returns the system temperature and gain for the selected W-band channel.

        The W-band receiver uses a CALSEQ where during a scan three different
        observations are made: sky, cold1 and cold2, from which the
        system temperature is derived.

        Parameters
        ----------
        scan : int or list of int
            Scan number(s) where CALSEQ is expected. See sdf.summary() to find the scan number(s).
            If multiple scans are used, an average Tsys is computed.
        fdnum : int, optional
            Feed to be used, 0 being the first.
            The default is 0.
        ifnum : int, optional
            IF to be used, 0 being the first.
            The default is 0.
        plnum : int, optional
            Polarization to be used, 0 being the first.
            The default is 0.
        freq : `~astropy.units.Quantity`, optional
            Set the frequency. By default the topocentric frequency of `ifnum` at `scan` will be used.
            It must have units of frequency.
        tcold : float, optional
            Set the cold temperature. By default it is computed as ``54 K - 0.6 K/GHz * (freq - 77 GHz)``.
        twarm : float, optional
            Set the warm temperature. By default it will use the value in the TWARM column of the SDFITS.
        apply_flags : bool, optional
            If True, apply flags before computing the system temperature.

        Returns
        -------
        tsys : float
            The system temperature, in K
        g : float
            The gain in K/counts

        Raises
        ------
        ValueError
            If `fdnum` is not 0 or 1.
        """

        tp_args = {
            "scan": scan,
            "ifnum": ifnum,
            "plnum": plnum,
            "fdnum": fdnum,
            "calibrate": True,
            "cal": False,
            "apply_flags": apply_flags,
        }
        vsky = self.gettp(CALPOSITION="Observing", **tp_args).timeaverage()
        vcold1 = self.gettp(CALPOSITION="Cold1", **tp_args).timeaverage()
        vcold2 = self.gettp(CALPOSITION="Cold2", **tp_args).timeaverage()

        # @todo ? there was a period when TWARM was recorded wrongly as 99 C, where TAMBIENT (in K) would be better.
        if twarm is None:
            twarm = vsky.meta["TWARM"]  # TWARM recorded in Kelvin when using Rcvr68_92 (W-Band).

        if freq is None:
            freq = vsky.spectral_axis.quantity.mean().to("GHz")

        if tcold is None:
            tcold = (54 - 0.6 * u.GHz**-1 * (freq - 77 * u.GHz)).value

        if fdnum == 0:
            g = (twarm - tcold) / mean_data(vcold2.data - vcold1.data)
        elif fdnum == 1:
            g = (twarm - tcold) / mean_data(vcold1.data - vcold2.data)
        else:
            raise ValueError(f"Illegal fdnum={fdnum} for a CALSEQ")

        tsys = mean_data(g * vsky.data)

        logger.debug(f"Twarm={twarm} Tcold={tcold}")
        logger.debug(f"IFNUM {ifnum} PLNUM {plnum} FDNUM {fdnum}")
        logger.debug(f"Tsys = {tsys}")
        logger.debug(f"Gain [K/counts] = {g}")

        return tsys, g

    def getvane(
        self,
        scan: int,
        fdnum: int,
        ifnum: int,
        plnum: int,
        t_cal: float | None = None,
        zenith_opacity: float | None = None,
        t_atm: float | None = None,
        t_warm: float | None = None,
        t_bkg: float = 2.725,
        apply_flags=True,
        **kwargs,
    ):
        """
        Return a `~dysh.spectra.vane.VaneSpectrum` used for calibrating observations with a vane.
        Uses the Equations provided in [1]_. For the most accurate results `zenith_opacity` and `tatm` should be provided.
        Otherwise, it will try to fetch these values from the GBT weather forecast scripts (only available at GBO).

        Parameters
        ----------
        scan : int
            Scan number for either the VANE object.
        fdnum : int
            The feed number.
        ifnum : int
            The intermediate frequency (IF) number.
        plnum : int
            The polarization number.
        t_cal : float, optional
            Calibration temperature. If no value is provided, but `zenith_opacity` and `tatm` are provided, then
            it will use Eq. (22) of [1]_. If `zenith_opacity` and `tatm` are not provided, it will first try to
            retrieve them using the weather forecasts (only available at GBO), if that fails it will use the
            ambient temperature, Eq. (23) of [1]_.
        zenith_opacity : float, optional
            Zenith opacity. If not provided it will try to fetch "Opacity" from the weather forecasts (only available at GBO).
        t_atm : float, optional
            Atmospheric temperature in K. If not provided it will try to fetch "Tatm" from the weather forecasts (only available at GBO).
        t_warm : float, optional
            Temperature of the VANE in K. If not provided it will use the value found in the "TWARM" column of the SDFITS for `scan`.
        t_bkg : float, optional
            Background temperature in K.
        apply_flags : bool, optional
            If True, apply flags before deriving the system temperature.

        Returns
        -------
        `~dysh.spectra.vane.VaneSpectrum`
            A `~dysh.spectra.vane.VaneSpectrum` object which can be used to calibrate observations with a vane.

        .. [1] `D. Frayer et al., "Calibration of Argus and the 4mm Receiver on the GBT" <https://ui.adsabs.harvard.edu/abs/2019nrao.reptE...1F/abstract>`_
        """

        vane = self.gettp(
            scan=scan,
            fdnum=fdnum,
            ifnum=ifnum,
            plnum=plnum,
            calibrate=True,
            cal=False,
            apply_flags=apply_flags,
            vane=True,
        ).timeaverage()

        return VaneSpectrum.from_spectrum(
            vane,
            scan,
            fdnum,
            ifnum,
            plnum,
            tcal=t_cal,
            zenith_opacity=zenith_opacity,
            tatm=t_atm,
            twarm=t_warm,
            tbkg=t_bkg,
        )

    def vanecal(
        self,
        scan,
        ifnum=0,
        plnum=0,
        fdnum=0,
        mode=2,
        tcal=None,
        zenith_opacity=None,
        tatm=None,
        twarm=None,
        tbkg=2.725,
        apply_flags=True,
        **kwargs,
    ):
        """
        Compute the system temperature from a VANE/SKY calibration cycle.
        Uses the Equations provided in [1]_. For the most accurate results `zenith_opacity` and `tatm` should be provided.

        Parameters
        ----------
        scan : int
            Scan number for either the SKY or VANE object.
            The pair will be found by the object name.
        ifnum : int
            The intermediate frequency (IF) number.
        plnum : int
            The polarization number.
        fdnum : int
            The feed number.
        mode : int, optional
            Mode of computing. See also `mean_tsys()`
            mode=0  Do the mean before the division
            mode=1  Do the mean after the division
            mode=2  Take a median of the inverse division
            The default is 2.
        tcal : float, optional
            Calibration temperature. If no value is provided, but `zenith_opacity` and `tatm` are provided, then
            it will use Eq. (22) of [1]_. If `zenith_opacity` and `tatm` are not provided, it will first try to
            retrieve them using the weather forecasts (only available at GBO), if that fails it will use the
            ambient temperature, Eq. (23) of [1]_.
        zenith_opacity : float, optional
            Zenith opacity. If not provided it will try to fetch "Opacity" from the weather forecasts (only available at GBO).
        tatm : float, optional
            Atmospheric temperature in K. If not provided it will try to fetch "Tatm" from the weather forecasts (only available at GBO).
        twarm : float, optional
            Temperature of the VANE in K. If not provided it will use the value found in the "TWARM" column of the SDFITS for `scan`.
        tbkg : float, optional
            Background temperature in K.
        apply_flags : bool, optional
            If True, apply flags before deriving the system temperature.

        Returns
        -------
        tsys : float
            System temperature in K.

        .. [1] `D. Frayer et al., "Calibration of Argus and the 4mm Receiver on the GBT" <https://ui.adsabs.harvard.edu/abs/2019nrao.reptE...1F/abstract>`_
        """
        t = set(self._index["OBJECT"][self._index["SCAN"] == scan])
        if len(t) > 1:
            raise TypeError(f"More than one OBJECT for scan {scan}")
        if t == {"VANE"}:
            vane_scan = scan
            sky_scan = scan + 1
        elif t == {"SKY"}:
            sky_scan = scan
            vane_scan = scan - 1

        vane = self.gettp(
            scan=vane_scan,
            fdnum=fdnum,
            ifnum=ifnum,
            plnum=plnum,
            calibrate=True,
            cal=False,
            apply_flags=apply_flags,
        ).timeaverage(use_wcs=False)
        sky = self.gettp(
            scan=sky_scan,
            fdnum=fdnum,
            ifnum=ifnum,
            plnum=plnum,
            calibrate=True,
            cal=False,
            apply_flags=apply_flags,
        ).timeaverage(use_wcs=False)

        if twarm is None:
            twarm = sky.meta["TWARM"] + 273.15  # TWARM is recorded in Celsius when using RcvrArray75_115 (Argus).

        if zenith_opacity is None:
            try:
                gbwf = GBTWeatherForecast()
                result = gbwf.fetch(
                    vartype="Opacity",
                    specval=sky.spectral_axis.quantity.mean(),
                    mjd=sky.obstime.mjd,
                )
                zenith_opacity = result[:, -1]
            except ValueError as e:
                logger.debug("Could not get forecasted zenith opacity ", e)

        if tatm is None:
            try:
                gbwf = GBTWeatherForecast()
                result = gbwf.fetch(vartype="Tatm", specval=sky.spectral_axis.quantity.mean(), mjd=sky.obstime.mjd)
                tatm = result[:, -1]
            except ValueError as e:
                logger.debug("Could not get forecasted atmospheric temperature ", e)

        if tcal is None:
            logger.info("No calibration temperature provided.")
            if zenith_opacity is None:
                tcal = sky.meta["TAMBIENT"]
                logger.info(
                    f"No zenith opacity provided. Will approximate the calibration temperature to the ambient temperature {tcal} K."
                )
            elif tatm is not None:
                airmass = GBTGainCorrection().airmass(sky.meta["ELEVATIO"] * u.deg, zd=False)
                tcal = (tatm - tbkg) + (twarm - tatm) * np.exp(zenith_opacity * airmass)

        match mode:
            case 0:
                mean_off = mean_data(sky.data)
                mean_dif = mean_data(vane.data - sky.data)
                tsys = tcal * mean_off / mean_dif
            case 1:
                tsys = tcal / mean_data((vane.data - sky.data) / sky.data)
            case 2:
                tsys = tcal / np.nanmedian((vane.data - sky.data) / sky.data)

        logger.debug(f"TCAL={tcal} K")
        logger.debug(f"mode={mode}")
        logger.debug(f"TSYS={tsys} K")

        return tsys

    def _get_bintable(self, df: pd.DataFrame) -> int:
        """
        Extracts the binary table from `df`.

        Parameters
        ----------
        df : `~pandas.DataFrame`
            The data frame to be used.

        Returns
        -------
        bintable : int
            The binary table index.

        Raises
        ------
        TypeError
            If there is more than one unique value in the "BINTABLE" column of `df`.
        """

        bintable = df["BINTABLE"].unique()
        # I do not know if this is possible, but just in case.
        if len(bintable) > 1:
            raise TypeError(
                "Selection crosses binary tables. Please provide more details during data selection (e.g., bintable=x)."
            )
        return bintable[0]

    def _get_refspec_tsys(self, refspec):
        """
        Find the system temperature in a `~dysh.spectra.spectrum.Spectrum`.
        It checks the meta attribute keys in the following order:
        "TSYS", "MEANTSYS", "WTTSYS"
        and returns the first not None value.
        """
        tsyskw = ["TSYS", "MEANTSYS", "WTTSYS"]
        for kw in tsyskw:
            tsys = refspec.meta.get(kw, None)
            if tsys is not None:
                break
        if tsys is None:
            raise ValueError(
                "Reference spectrum has no system temperature in its metadata.  Solve with refspec.meta['TSYS']=value or add parameter `t_sys` to getps/getsigref."
            )
        return tsys

    def _get_tcal(self, tcal):
        """
        Retrieve the value of TCAL.

        Raises
        ------
        ValueError
            If there's more than one value for TCAL.
        """
        tcal_set = tcal.unique()
        if len(tcal_set) > 1:
            raise ValueError(f"More than one value for TCAL: {tcal_set}")
        return tcal_set[0]

    def _vane_setup(self, vane, fdnum, ifnum, plnum, units, zenith_opacity, t_warm, t_atm, t_bkg, t_cal, apply_flags):
        """
        Set up a `~dysh.spectra.vane.VaneSpectrum` for use in the calibration routines.
        It also handles the hacks needed to get the units correctly when using a vane.
        """

        requested_units = copy.copy(units)  # Keep track of what the user wants.
        if units.lower() not in ["ta*", "flux"]:
            logger.info("Vane calibrated data will be calibrated to Ta* units by default.")
            units = "Ta"  # Set to Ta to disable scaling during calibration. Vane calibrates to Ta* by default.
            requested_units = (
                "Ta*"  # Force to Ta* if the input was Ta. This will be used at the end to scale the ScanBase.
            )
        if isinstance(vane, VaneSpectrum):
            if (
                zenith_opacity is not None
                and vane._zenith_opacity is not None
                and zenith_opacity != vane._zenith_opacity
            ):
                vane._zenith_opacity = zenith_opacity
                logger.info(
                    f"Zenith opacity provided and present in vane, but they do not match. Will use the value provided ({zenith_opacity} nepers)"
                )
            elif zenith_opacity is None and vane._zenith_opacity is not None:
                zenith_opacity = vane._zenith_opacity
                logger.info(f"Will use a zenith opacity of {zenith_opacity} nepers. Taken from vane.")
            if t_warm is not None:
                logger.info(
                    "t_warm provided, but not used. To change this value, please create a new VaneSpectrum or call with vane as a scan number."
                )
            if t_atm is not None:
                logger.info(
                    "t_atm provided, but not used. To change this value, please create a new VaneSpectrum or call with vane as a scan number."
                )
            if t_cal is not None:
                logger.info(
                    "t_cal provided, but not used. To change this value, please create a new VaneSpectrum or call with vane as a scan number."
                )
        elif isinstance(vane, int):
            vane = self.getvane(
                scan=vane,
                fdnum=fdnum,
                ifnum=ifnum,
                plnum=plnum,
                zenith_opacity=zenith_opacity,
                t_warm=t_warm,
                t_atm=t_atm,
                t_bkg=t_bkg,
                t_cal=t_cal,
                apply_flags=apply_flags,
            )
        else:
            raise TypeError(f"vane must be an int or VaneSpectrum. Got a {type(vane)} instead.")

        return vane, units, requested_units, zenith_opacity

    def _check_vane_and_t_sys_args(self, vane, t_sys):
        """Check if both arguments were provided."""
        if t_sys is not None and vane is not None:
            logger.warning("Both t_sys and vane provided. Ignoring t_sys.")
            t_sys = None

    def _set_scale_vane(self, scan, units, zenith_opacity):
        """
        Force scale to be Ta* and then scale as needed.
        This is used for calibration with a vane, because
        it calibrates to Ta* and we do not have (?) a
        way of handling this without this kludge.

        Parameters
        ----------
        scan : `~dysh.spectra.scan.ScanBase`
            Scan to have its scale updated.
        units : str
            Units of the updated `scan`.
        zenith_opacity : float
            Zenith opacity in nepers.
        """
        scan._tscale_fac[:] = 1.0
        scan._tscale = "ta*"
        scan._update_scale_meta()
        scan.scale(units, zenith_opacity=zenith_opacity)


class GBTOffline(GBTFITSLoad):
    """
    GBTOffline('foo')   connects to a GBT project 'foo' using GBTFITSLoad

    Note project directories are assumed to exist in /home/sdfits
    or whereever dysh_data thinks your /home/sdfits lives.

    Also note, as in GBTIDL, one can use SDFITS_DATA instead of DYSH_DATA

    Use dysh_data('?') to display all filenames in the "sdfits" area.

    """

    @log_call_to_history
    def __init__(self, fileobj, *args, **kwargs):
        self._offline = fileobj
        self._filename = dysh_data(fileobj)
        GBTFITSLoad.__init__(self, self._filename, *args, **kwargs)


# NOTE: if GBTFITSLoad has new functions added, they may need to be added in GBTOnline() as well
#       these two variables with _check_functions() will warn in runtime, but fail in pytest
#       If you add more to _skip_functions, deduct the number in _need_functions
_skip_functions = ["velocity_convention", "velocity_frame"]
_need_functions = 55


def _check_functions(verbose=False):
    """
    check if number of functions in GBTFITSLoad() didn't change from
    the last time we (manually) recorded this.
    """
    fns = inspect.getmembers(GBTFITSLoad, predicate=inspect.isfunction)
    need = _need_functions
    n = 0
    for i in range(len(fns)):
        fn = fns[i][0]
        if fn[0] == "_":
            continue
        if fn in _skip_functions:
            continue
        n = n + 1
        if verbose:
            print(n, fn)
    if n != need:
        # this means GBTOnline may need to have the new
        logger.warning(f"GBTOnline: parent GBTFITSLoad() was expected have {need} functions, but found {n}.")
    # return values for tests
    return (need, n)


class GBTOnline(GBTFITSLoad):
    """
    GBTOnline('foo')   monitors project 'foo' as if it could be online
    GBTOnline()        monitors for new projects and connects, and refreshes when updated

    Note project directories are assumed to exist in /home/sdfits
    or whereever dysh_data thinks your /home/sdfits lives.

    Also note, as in GBTIDL, one can use SDFITS_DATA instead of DYSH_DATA

    Use dysh_data('?') to display all filenames in the "sdfits" area.

    """

    @log_call_to_history
    def __init__(self, fileobj=None, *args, **kwargs):
        self._online = fileobj
        self._args = args
        self._kwargs = kwargs
        self._platform = platform.system()  # cannot update in "Windows", see #447
        _check_functions()  # check if # functions if GBTFITSLoad didn't change
        if fileobj is not None:
            self._online_mode = 1  # monitor this file
            if os.path.isdir(fileobj):
                GBTFITSLoad.__init__(self, fileobj, *args, **kwargs)
            else:
                self._online = dysh_data(fileobj)
                GBTFITSLoad.__init__(self, self._online, *args, **kwargs)
            logger.info(f"Connecting to explicit file: {self._online} - will be monitoring this")

        else:
            self._online_mode = 2  #  monitor all files?
            logger.debug("Testing online mode, finding most recent file")
            if "SDFITS_DATA" in os.environ:
                logger.debug("warning: using SDITS_DATA")
                sdfits_root = os.environ["SDFITS_DATA"]
            elif "DYSH_DATA" in os.environ:
                sdfits_root = os.environ["DYSH_DATA"] + "/sdfits"
                logger.debug("warning: using DYSH_DATA")
            else:
                sdfits_root = "/home/sdfits"
            logger.debug(f"Using SDFITS_DATA {sdfits_root}")

            if not os.path.isdir(sdfits_root):  # @todo shouldn't this be an exception?
                logger.info(f"Cannot find {sdfits_root}")
                return None

            # 1. check the status_file ?
            status_file = "sdfitsStatus.txt"
            if os.path.exists(sdfits_root + "/" + status_file):
                logger.debug(f"Found {status_file} but not using it yet")

            # 2. visit each directory where the final leaf contains fits files, and find the most recent one
            n = 0
            mtime_max = 0
            for dirname, subdirs, files in os.walk(sdfits_root):
                if len(subdirs) == 0:
                    n = n + 1
                    for fname in files:
                        if fname.split(".")[-1] == "fits":
                            mtime = os.path.getmtime(dirname + "/" + fname)
                            if mtime > mtime_max:
                                mtime_max = mtime
                                project = dirname
                            break
            if n == 0:
                return None

            self._online = project
            GBTFITSLoad.__init__(self, self._online, *self._args, **self._kwargs)
            self._mtime = os.path.getmtime(self.filenames()[0])

        # we only test the first filename in the list, assuming they're all being written

        self._mtime = os.path.getmtime(self.filenames()[0])
        for f in self.filenames():
            self._mtime = max(self._mtime, os.path.getmtime(f))
        delta = (time.time() - self._mtime) / 60.0

        logger.info(f"Connected to: {self._online}")
        logger.info(f"Data has not been updated in {delta:.2f} minutes.")
        # end of __init__

    def _reload(self, force=False):
        """force a reload of the latest"""
        if self._platform == "Windows":
            logger.warning("Cannot reload on Windows, see issue #447")
            return
        if not force:
            for f in self.filenames():
                mtime = max(self._mtime, os.path.getmtime(f))
            if mtime > self._mtime:
                self._mtime = mtime
                logger.debug("NEW MTIME:", self._mtime)
                force = True
        if force:
            logger.info(f"Reload {self._online}")
            GBTFITSLoad.__init__(self, self._online, *self._args, **self._kwargs)
        return force

    # examples of catchers for reloading

    def summary(self, *args, **kwargs):
        self._reload()
        return super().summary(*args, **kwargs)

    def get_summary(self, *args, **kwargs):
        self._reload()
        return super().get_summary(*args, **kwargs)

    def write(self, *args, **kwargs):
        self._reload()
        return super().write(*args, **kwargs)

    def gettp(self, *args, **kwargs):
        self._reload()
        return super().gettp(*args, **kwargs)

    def getsigref(self, *args, **kwargs):
        self._reload()
        return super().getsigref(*args, **kwargs)

    def getps(self, *args, **kwargs):
        self._reload()
        return super().getps(*args, **kwargs)

    def getnod(self, *args, **kwargs):
        self._reload()
        return super().getnod(*args, **kwargs)

    def getfs(self, *args, **kwargs):
        self._reload()
        return super().getfs(*args, **kwargs)

    def subbeamnod(self, *args, **kwargs):
        self._reload()
        return super().subbeamnod(*args, **kwargs)

    def vanecal(self, *args, **kwargs):
        self._reload()
        return super().vanecal(*args, **kwargs)

    def calseq(self, *args, **kwargs):
        self._reload()
        return super().calseq(*args, **kwargs)

    def gettcal(self, *args, **kwargs):
        self._reload()
        return super().gettcal(*args, **kwargs)


def _parse_tsys(tsys: float | np.ndarray | list | dict, scans: list) -> dict:
    """
    Parse the system temperatures for a list of scans.

    Parameters
    ----------
    tsys : float or `~numpy.ndarray` or list or dict
        The system temperature(s)
    scans : list
        list of scan numbers associated with the system temperature(s)

    Raises
    ------
    TypeError
        If there is a mismatch between number of system temperatures and number of scans

    Returns
    -------
    dict
        Dictionary of system temperatures with scan number as keys
    """
    if isinstance(tsys, numbers.Real):
        tsys = _tsys_1Darray_to_dict(tsys, scans)
    if isinstance(tsys, list):
        tsys = np.array(tsys)
    if isinstance(tsys, np.ndarray):
        if tsys.ndim <= 1:
            tsys = _tsys_1Darray_to_dict(tsys, scans)
        elif tsys.ndim == 2:
            tsys = _tsys_2Darray_to_dict(tsys, scans)
    if isinstance(tsys, dict):
        # Check that there is one entry for every scan.
        if list(tsys.keys()) != list(scans):
            missing = set(scans) - set(tsys.keys())
            raise TypeError(f"Missing system temperature for scan(s): {','.join(map(str, missing))}")
        tsys = _tsys_dict_to_dict(tsys, scans)

    return tsys


def _tsys_1Darray_to_dict(tsys, scans):
    """Convert 1D array of system temperatures to a dictionary with scan number as keys"""
    tsys_dict = {}
    for scan in scans:
        tsys_dict[scan] = np.vstack((tsys, tsys))
    return tsys_dict


def _tsys_2Darray_to_dict(tsys, scans):
    tsys_dict = {}
    for scan in scans:
        tsys_dict[scan] = np.vstack((tsys[0], tsys[1]))
    return tsys_dict


def _tsys_dict_to_dict(tsys, scans):
    tsys_dict = {}
    for scan in scans:
        try:
            len(tsys[scan])
        except TypeError:
            tsys[scan] = [tsys[scan]]
        if len(tsys[scan]) < 2:
            tsys_dict[scan] = np.vstack((tsys[scan], tsys[scan]))
        else:
            tsys_dict[scan] = tsys[scan]
    return tsys_dict
