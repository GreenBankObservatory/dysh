"""Load SDFITS files produced by the Green Bank Telescope"""

import copy
import os
import platform
import time
import warnings
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits

from astropy import units as u
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
    SubBeamNodScan,
    TPScan,
)
from ..util import (
    Flag,
    Selection,
    consecutive,
    convert_array_to_mask,
    eliminate_flagged_rows,
    keycase,
    select_from,
    uniq,
)
from ..util.files import dysh_data
from ..util.selection import Flag, Selection
from . import conf
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

    """

    @log_call_to_history
    def __init__(self, fileobj, source=None, hdu=None, skipflags=False, **kwargs):
        kwargs_opts = {"index": True, "verbose": False, "fix_ka": True}  # only set to False for performance testing.
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
        elif path.is_dir():
            logger.debug(f"Treating given path {path} as a directory")
            # Find all the FITS files in the directory and sort alphabetically
            # because e.g., VEGAS does A,B,C,D,E
            for f in sorted(path.glob("*.fits")):
                logger.debug(f"Selecting {f} to load")
                if kwargs.get("verbose", None):
                    print(f"Loading {f}")
                self._sdf.append(SDFITSLoad(f, source, hdu, **kwargs_opts))
            if len(self._sdf) == 0:  # fixes issue 381
                raise Exception(f"No FITS files found in {fileobj}.")
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
            self._create_index_if_needed(skipflags)
            self._update_radesys()
        if kwargs_opts["fix_ka"]:
            self._fix_ka_rx_if_needed()
        # We cannot use this to get mmHg as it will disable all default astropy units!
        # https://docs.astropy.org/en/stable/api/astropy.units.cds.enable.html#astropy.units.cds.enable
        # cds.enable()  # to get mmHg

        # ushow/udata depend on the index being present, so check that index is created.
        if kwargs.get("verbose", None) and kwargs_opts["index"]:
            print("==GBTLoad %s" % fileobj)
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

    def index(self, hdu=None, bintable=None, fitsindex=None):
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

    def summary(self, scans=None, verbose=False, show_index=True):  # selected=False
        # From GBTIDL:
        # Intended to work with un-calibrated GBT data and is
        # likely to give confusing results for other data.  For other data,
        # list is usually more useful.   @todo what's the dysh eqv. of list ?
        #
        # @todo perhaps return as a astropy.Table then we can have units
        """
        Create a summary of the input dataset.
        If `verbose=False` (default), some numeric data
        (e.g., RESTFREQ, AZIMUTH, ELEVATIO) are
        averaged over the records with the same scan number.

        Parameters
        ----------
        scans : int or 2-tuple
            The scan(s) to use. A 2-tuple represents (beginning, ending) scans. Default: show all scans
        verbose: bool
            If True, list every record, otherwise return a compact summary.
        show_index: bool
            If True, show the DataFrame index column. Default: False

        Returns
        -------
        summary : `~pandas.DataFrame`
            Summary of the data as a DataFrame.

        """
        # @todo allow user to change show list
        # @todo set individual format options on output by
        # changing these to dicts(?)
        #
        pd.set_option("display.max_rows", conf.summary_max_rows)
        # 'show' is fragile because anything we might need to query in 'uf' below in
        # order to do a calculation,  whether we want to show it, or not must be in 'show.'
        # (e.g. PROCSIZE is needed to calculate n_integrations).
        show = [
            "SCAN",
            "OBJECT",
            "VELOCITY",
            "PROC",
            "PROCSEQN",
            "PROCSIZE",
            "RESTFREQ",
            "DOPFREQ",
            "IFNUM",
            "FEED",
            "AZIMUTH",
            "ELEVATIO",
            "FDNUM",
            "INTNUM",
            "PLNUM",
            "SIG",
            "CAL",
            "DATE-OBS",
        ]
        comp_colnames = [
            "SCAN",
            "OBJECT",
            "VELOCITY",
            "PROC",
            "PROCSEQN",
            "RESTFREQ",
            "DOPFREQ",
            "# IF",
            "# POL",
            "# INT",
            "# FEED",
            "AZIMUTH",
            "ELEVATIO",
        ]
        # In the process, some columns get cast to floats or others. Make sure we cast them
        # back to an appropriate data type before return.
        col_dtypes = {"SCAN": int, "PROCSEQN": int}
        uncompressed_df = None
        self._create_index_if_needed()
        # make a copy here because we can't guarantee if this is a
        # view or a copy without it. See https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#returning-a-view-versus-a-copy
        _df = self[show].copy()
        _df.loc[:, "VELOCITY"] /= 1e3  # convert to km/s
        _df["RESTFREQ"] = _df["RESTFREQ"] / 1.0e9  # convert to GHz
        _df["DOPFREQ"] = _df["DOPFREQ"] / 1.0e9  # convert to GHz
        if scans is not None:
            if type(scans) == int:
                scans = [scans]
            if len(scans) == 1:
                scans = [scans[0], scans[0]]  # or should this be [scans[0],lastscan]?
            _df = self._select_scans(scans, _df).filter(show)
            if uncompressed_df is None:
                uncompressed_df = _df
            else:  # no longer used
                uncompressed_df = pd.concat([uncompressed_df, _df])
        else:
            if uncompressed_df is None:
                uncompressed_df = _df.filter(show)
            else:  # no longer used
                uncompressed_df = pd.concat([uncompressed_df, _df.filter(show)])
        if verbose:
            uncompressed_df = uncompressed_df.astype(col_dtypes)
            return uncompressed_df
        # do the work to compress the info
        # in the dataframe on a scan basis
        compressed_df = pd.DataFrame(columns=comp_colnames)
        scanset = set(uncompressed_df["SCAN"])
        avg_cols = ["SCAN", "VELOCITY", "PROCSEQN", "RESTFREQ", "DOPFREQ", "AZIMUTH", "ELEVATIO"]
        for s in scanset:
            uf = select_from("SCAN", s, uncompressed_df)
            # for some columns we will display
            # the mean value
            ser = uf.filter(avg_cols).mean(numeric_only=True)
            ser.rename("filtered ser")
            # for others we will count how many there are
            nIF = uf["IFNUM"].nunique()
            nPol = uf["PLNUM"].nunique()
            nfeed = uf["FEED"].nunique()
            # For counting integrations, take care of out-of-sync samplers by just
            # looking at the first instance of FEED, PLNUM, and IFNUM.
            uf_int = select_from("FEED", uf["FEED"].iloc[0], uf)
            uf_int = select_from("PLNUM", uf_int["PLNUM"].iloc[0], uf_int)
            uf_int = select_from("IFNUM", uf_int["IFNUM"].iloc[0], uf_int)
            nint = len(set(uf_int["DATE-OBS"]))  # see gbtidl io/line_index__define.pro
            obj = list(set(uf["OBJECT"]))[0]  # We assume they are all the same!
            proc = list(set(uf["PROC"]))[0]  # We assume they are all the same!
            s2 = pd.Series(
                [obj, proc, nIF, nPol, nint, nfeed],
                name="uniqued data",
                index=["OBJECT", "PROC", "# IF", "# POL", "# INT", "# FEED"],
            )
            ser = pd.concat([ser, s2]).reindex(comp_colnames)
            ser.rename("appended ser")
            compressed_df = pd.concat([compressed_df, ser.to_frame().T], ignore_index=True)
        compressed_df = compressed_df.astype(col_dtypes)
        if not show_index:
            print(compressed_df.to_string(index=False))
            # return compressed_df.style.hide(axis="index")
        else:
            return compressed_df

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
        (convention, frame) = decode_veldef(veldef)
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
        (convention, frame) = decode_veldef(veldef)
        return frame

    def _select_scans(self, scans, df):
        return df[(df["SCAN"] >= scans[0]) & (df["SCAN"] <= scans[1])]

    # @todo maybe move all selection/flag methods to sdfitsload after adding Selection/Flag
    # to sdfitsload
    # @todo maybe write a Delegator class to autopass to Selection.
    # See, e.g., https://michaelcho.me/article/method-delegation-in-python/
    @log_call_to_history
    def select(self, tag=None, **kwargs):
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
            key : str
                The key  (SDFITS column name or other supported key)
            value : any
                The value to select

        """
        self._selection.select(tag=tag, **kwargs)

    @log_call_to_history
    def select_range(self, tag=None, **kwargs):
        """
        Select a range of inclusive values for a given key(s).
        e.g., `key1 = (v1,v2), key2 = (v3,v4), ...`
        will select data  `v1 <= data1 <= v2, v3 <= data2 <= v4, ... `
        Upper and lower limits may be given by setting one of the tuple values
        to None. e.g., `key1 = (None,v1)` for an upper limit `data1 <= v1` and
        `key1 = (v1,None)` for a lower limit `data >=v1`.  Lower
        limits may also be specified by a one-element tuple `key1 = (v1,)`.
        See `~dysh.util.selection.Selection`.

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
        self._selection.select_range(tag=tag, **kwargs)

    @log_call_to_history
    def select_within(self, tag=None, **kwargs):
        """
        Select a value within a plus or minus for a given key(s).
        e.g. `key1 = [value1,epsilon1], key2 = [value2,epsilon2], ...`
        Will select data
        `value1-epsilon1 <= data1 <= value1+epsilon1,`
        `value2-epsilon2 <= data2 <= value2+epsilon2,...`

        See `~dysh.util.selection.Selection`.

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
        self._selection.select_within(tag=tag, **kwargs)

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
            # tuples also work, though can be harder for a human to read
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
    def flag(self, tag=None, **kwargs):
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
        key : str
            The key  (SDFITS column name or other supported key)
        value : any
            The value to select

        """
        self._flag.flag(tag=tag, **kwargs)

    @log_call_to_history
    def flag_range(self, tag=None, **kwargs):
        """
        Flag a range of inclusive values for a given key(s).
        e.g., `key1 = (v1,v2), key2 = (v3,v4), ...`
        will select data  `v1 <= data1 <= v2, v3 <= data2 <= v4, ...`

        Upper and lower limits may be given by setting one of the tuple values
        to None. e.g., `key1 = (None,v1)` for an upper limit `data1 <= v1` and
        `key1 = (v1,None)` for a lower limit `data >=v1`.  Lower
        limits may also be specified by a one-element tuple `key1 = (v1,)`.
        See `~dysh.util.selection.Flag`.

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
        self._flag.flag_range(tag=tag, **kwargs)

    @log_call_to_history
    def flag_within(self, tag=None, **kwargs):
        """
        Flag a value within a plus or minus for a given key(s).
        e.g. `key1 = [value1,epsilon1], key2 = [value2,epsilon2], ...`
        Will select data
        `value1-epsilon1 <= data1 <= value1+epsilon1,`
        `value2-epsilon2 <= data2 <= value2+epsilon2,...`

        See `~dysh.util.selection.Flag`.

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
        self._flag.flag_within(tag=tag, **kwargs)

    @log_call_to_history
    def flag_channel(self, channel, tag=None):
        """
        Select channels and/or channel ranges. These are NOT used in :meth:`final`
        but rather will be used to create a mask for
        flagging. Single arrays/tuples will be treated as channel lists;
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
            # tuples also work, though can be harder for a human to read
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
            for i, ((fi, bi), g) in enumerate(dfs):
                chan_mask = convert_array_to_mask(chan, self._sdf[fi].nchan(bi))
                rows = g["ROW"].to_numpy()
                logger.debug(f"Applying {chan} to {rows=}")
                logger.debug(f"{np.where(chan_mask)}")
                # print(f"Applying {chan} to {rows=}")
                # print(f"{np.where(chan_mask)}")
                self._sdf[fi]._flagmask[bi][rows] |= chan_mask

    @log_call_to_history
    def clear_flags(self):
        """Clear all flags for these data"""
        for sdf in self._sdf:
            sdf._init_flags()
        self._flag.clear()

    def _create_index_if_needed(self, skipflags=False):
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
                self.flags.read(flagfile, fitsindex=fi)
                found_flags = True
        if found_flags:
            print("Flags were created from existing flag files. Use GBTFITSLoad.flags.show() to see them.")

    def _construct_procedure(self):
        """
        Construct the procedure string (PROC) from OBSMODE and add it to the index (i.e., a new SDFITS column).
        OBSTYPE and SUBOBSMODE are also created here.  OBSMODE has the form like 'PROC:OBSTYPE:SUBOBSMODE', e.g.
        OnOff:PSWITCHON:TPWCAL.

        """
        if self._selection is None:
            warnings.warn("Couldn't construct procedure string: index is not yet created.")
            return
        if "OBSMODE" not in self._index:
            warnings.warn("Couldn't construct procedure string: OBSMODE is not in index.")
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
            warnings.warn("Couldn't construct integration number: index is not yet created.")
            return

        # check it hasn't been constructed before.
        if "INTNUM" in self._index:
            return
        # check that GBTIDL didn't write it out at some point.
        if "INT" in self._index:
            self._index.rename(columns={"INT": "INTNUM"}, inplace=True)
            for s in self._sdf:
                s._rename_binary_table_column("int", "intnum")
            return

        intnumarray = np.empty(len(self._index), dtype=int)
        # Leverage pandas to group things by scan and observing time.
        dfs = self._index.groupby(["SCAN"])
        for name, group in dfs:
            dfst = group.groupby("DATE-OBS")
            intnums = np.arange(0, len(dfst.groups))
            for i, (n, g) in enumerate(dfst):
                idx = g.index
                intnumarray[idx] = intnums[i]
        self._index["INTNUM"] = intnumarray
        self._flag["INTNUM"] = intnumarray

    def info(self):
        """Return information on the HDUs contained in this object. See :meth:`~astropy.HDUList/info()`"""
        for s in self._sdf:
            s.info()

    def _common_selection(self, ifnum, plnum, fdnum, **kwargs):
        """Do selection and flag application common to all calibration methods.
        Flags are not applied unless selection results in non-zero length data selection.

        Parameters
        ----------
        fdnum: int
            The feed number
        ifnum : int
            The IF number
        plnum : int
            The polarization number

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
        elif type(scans) == int:
            scans = list([scans])
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

    @log_call_to_result
    def gettp(
        self,
        fdnum,
        ifnum,
        plnum,
        sig=None,
        cal=None,
        calibrate=True,
        bintable=None,
        smoothref=1,
        apply_flags=True,
        **kwargs,
    ):
        """
        Get a total power scan, optionally calibrating it.

        Parameters
        ----------
        fdnum: int
            The feed number
        ifnum : int
            The IF number
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
        bintable : int, optional
            Limit to the input binary table index. The default is None which means use all binary tables.
        smooth_ref: int, optional
            the number of channels in the reference to boxcar smooth prior to calibration
        apply_flags : boolean, optional.  If True, apply flags before calibration.
            See :meth:`apply_flags`. Default: True
        **kwargs : dict
            Optional additional selection  keyword arguments, typically
            given as key=value, though a dictionary works too.
            e.g., `source="NGC132", intnum=range(20)` etc.

        Returns
        -------
        data : `~spectra.scan.ScanBlock`
            A ScanBlock containing one or more `~spectra.scan.TPScan`

        """
        (scans, _sf) = self._common_selection(fdnum=fdnum, ifnum=ifnum, plnum=plnum, apply_flags=apply_flags, **kwargs)
        TF = {True: "T", False: "F"}
        scanblock = ScanBlock()
        calrows = {}
        for i in range(len(self._sdf)):
            _df = select_from("FITSINDEX", i, _sf)
            for scan in scans:
                _sifdf = select_from("SCAN", scan, _df)
                dfcalT = select_from("CAL", "T", _sifdf)
                dfcalF = select_from("CAL", "F", _sifdf)
                calrows["ON"] = list(dfcalT["ROW"])
                calrows["OFF"] = list(dfcalF["ROW"])
                # print("PJT CALROWS: ",calrows["ON"] ,calrows["OFF"])
                if len(calrows["ON"]) != len(calrows["OFF"]):
                    if len(calrows["ON"]) > 0:
                        raise Exception(f'unbalanced calrows {len(calrows["ON"])} != {len(calrows["OFF"])}')
                    # else: print("Warning: hacking gettp with no calrows")
                # sig and cal are treated specially since
                # they are not in kwargs and in SDFITS header
                # they are not booleans but chars
                if sig is not None:
                    _sifdf = select_from("SIG", TF[sig], _sifdf)
                # if cal is not None:
                #    df = select_from("CAL", TF[cal], df)
                # the rows with the selected sig state and all cal states
                tprows = list(_sifdf["ROW"])
                logger.debug(f"TPROWS len={len(tprows)}")
                logger.debug(f"CALROWS on len={len(calrows['ON'])}")
                logger.debug(f"fitsindex={i}")
                # print("PJT TPROWS", tprows)
                if len(tprows) == 0:
                    continue
                g = TPScan(
                    self._sdf[i],
                    scan,
                    sig,
                    cal,
                    tprows,
                    calrows,
                    fdnum=fdnum,
                    ifnum=ifnum,
                    plnum=plnum,
                    bintable=bintable,
                    calibrate=calibrate,
                    smoothref=smoothref,
                    apply_flags=apply_flags,
                )
                g.merge_commentary(self)
                scanblock.append(g)
        if len(scanblock) == 0:
            raise Exception("Didn't find any scans matching the input selection criteria.")
        scanblock.merge_commentary(self)
        return scanblock
        # end of gettp()

    @log_call_to_result
    def getps(
        self,
        fdnum,
        ifnum,
        plnum,
        calibrate=True,
        bintable=None,
        smoothref: int = 1,
        apply_flags: str = True,
        bunit: str = "ta",
        zenith_opacity: float = None,
        **kwargs,
    ):
        """
        Retrieve and calibrate position-switched data.

        Parameters
        ----------
        fdnum: int
            The feed number
        ifnum : int
            The IF number
        plnum : int
            The polarization number
        calibrate : boolean, optional
            Calibrate the scans. The default is True.
        bintable : int, optional
            Limit to the input binary table index. The default is None which means use all binary tables.
            (This keyword should eventually go away)
        smooth_ref: int, optional
            the number of channels in the reference to boxcar smooth prior to calibration
        apply_flags : boolean, optional.  If True, apply flags before calibration.
            See :meth:`apply_flags`. Default: True
        bunit : str, optional
            The brightness scale unit for the output scan, must be one of (case-insensitive)
                    - 'ta'  : Antenna Temperature
                    - 'ta*' : Antenna temperature corrected to above the atmosphere
                    - 'jy'  : flux density in Jansky
            If 'ta*' or 'jy' the zenith opacity must also be given. Default:'ta'
        zenith_opacity: float, optional
            The zenith opacity to use in calculating the scale factors for the integrations.  Default:None

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
        scanblock : `~spectra.scan.ScanBlock`
            ScanBlock containing one or more `~spectra.scan.PSScan`.

        """
        ScanBase._check_bunit(bunit)
        if bunit.lower() != "ta" and zenith_opacity is None:
            raise ValueError("Can't scale the data without a valid zenith opacity")

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
        # @todo pjt  two additions in this merge ?
        if True:
            som = uniq(_sf["SUBOBSMODE"])
            if len(som) > 1:
                raise Exception(f"Multiple SUBOBSMODE present, cannot deal with this yet {som}")
            if som[0] == "TPNOCAL":
                self._tpnocal = True
                raise Exception("Cannot deal with TPNOCAL yet")
        scanblock = ScanBlock()
        for i in range(len(self._sdf)):
            _df = select_from("FITSINDEX", i, _sf)
            scanlist = self._common_scan_list_selection(scans, _df, prockey=prockey, procvals=procvals, check=False)
            if len(scanlist["ON"]) == 0 or len(scanlist["OFF"]) == 0:
                logger.debug(f"scans {scans} not found, continuing")
                continue
            rows = {}
            # loop over scan pairs
            c = 0
            for on, off in zip(scanlist["ON"], scanlist["OFF"]):
                _ondf = select_from("SCAN", on, _df)
                _offdf = select_from("SCAN", off, _df)
                # rows["ON"] = list(_ondf.index)
                # rows["OFF"] = list(_offdf.index)
                rows["ON"] = list(_ondf["ROW"])
                rows["OFF"] = list(_offdf["ROW"])
                for key in rows:
                    if len(rows[key]) == 0:
                        raise Exception(f"{key} scans not found in scan list {scans}")
                # do not pass scan list here. We need all the cal rows. They will
                # be intersected with scan rows in PSScan
                calrows = {}
                dfcalT = select_from("CAL", "T", _df)
                dfcalF = select_from("CAL", "F", _df)
                # calrows["ON"] = list(dfcalT.index)
                # calrows["OFF"] = list(dfcalF.index)
                calrows["ON"] = list(dfcalT["ROW"])
                calrows["OFF"] = list(dfcalF["ROW"])
                d = {"ON": on, "OFF": off}
                g = PSScan(
                    self._sdf[i],
                    scan=d,
                    scanrows=rows,
                    calrows=calrows,
                    fdnum=fdnum,
                    ifnum=ifnum,
                    plnum=plnum,
                    bintable=bintable,
                    calibrate=calibrate,
                    smoothref=smoothref,
                    apply_flags=apply_flags,
                    bunit=bunit,
                    zenith_opacity=zenith_opacity,
                )
                g.merge_commentary(self)
                scanblock.append(g)
                c = c + 1
        if len(scanblock) == 0:
            raise Exception("Didn't find any scans matching the input selection criteria.")
        scanblock.merge_commentary(self)
        return scanblock
        # end of getps()

    @log_call_to_result
    def getnod(
        self,
        ifnum,
        plnum,
        fdnum=None,
        calibrate=True,
        bintable=None,
        smoothref=1,
        apply_flags=True,
        t_sys=None,
        nocal=False,
        bunit="ta",
        zenith_opacity=None,
        **kwargs,
    ):
        """
        Retrieve and calibrate nodding data.

        Parameters
        ----------
        ifnum : int
            The IF number
        plnum : int
            The polarization number
        fdnum:  2-tuple, optional
            The feed numbers. A pair of feed numbers may be given to choose different nodding beams than were used to obtain the observations.  Default: None which means use the beams found in the data.
        calibrate : boolean, optional
            Calibrate the scans.
            The default is True.
        bintable : int, optional
            Limit to the input binary table index. The default is None which means use all binary tables.
            (This keyword should eventually go away)
        smooth_ref: int, optional
            the number of channels in the reference to boxcar smooth prior to calibration
        apply_flags : boolean, optional.  If True, apply flags before calibration.
            See :meth:`apply_flags`. Default: True
        t_sys : float, optional
            System temperature. If provided, it overrides the value computed using the noise diode.
            If no noise diode is fired, and `t_sys=None`, then the column "TSYS" will be used instead.
        nocal : bool, optional
            Is the noise diode being fired? False means the noise diode was firing.
            By default it will figure this out by looking at the "CAL" column.
            It can be set to True to override this. Default: False
        **kwargs : dict
            Optional additional selection keyword arguments, typically
            given as key=value, though a dictionary works too.
            e.g., `ifnum=1, plnum=2` etc.
            For multi-beam with more than 2 beams, fdnum=[BEAM1,BEAM2] must be selected,
            unless the data have been properly taggeed using PROCSCAN which BEAM1 and BEAM2 are.

        Raises
        ------
        Exception
            If scans matching the selection criteria are not found.

        Returns
        -------
        scanblock : `~spectra.scan.ScanBlock`
            ScanBlock containing one or more `~spectra.scan.NodScan`.

        """

        def get_nod_beams(sdf):
            """find the two nodding beams if user did not specify them"""
            kb = ["DATE-OBS", "SCAN", "IFNUM", "PLNUM", "FDNUM", "PROCSCAN", "FEED", "SRFEED", "FEEDXOFF", "FEEDEOFF"]
            a = sdf._index[kb]
            b = a.loc[a["FEEDXOFF"] == 0.0]
            c = b.loc[b["FEEDEOFF"] == 0.0]
            d1 = c.loc[c["PROCSCAN"] == "BEAM1"]
            d2 = c.loc[c["PROCSCAN"] == "BEAM2"]
            if len(d1["FDNUM"].unique()) == 1 and len(d2["FDNUM"].unique()) == 1:
                beam1 = d1["FDNUM"].unique()[0]
                beam2 = d2["FDNUM"].unique()[0]
                return [beam1, beam2]
            else:
                # one more attempt (this can happen if PROCSCAN contains "Unknown")
                # ugh, is it possible that BEAM1 and BEAM2 are switched here, given how we unique() ?
                if len(c["FEED"].unique()) == 2:
                    logger.debug("get_nod_beams rescued")
                    b = c["FEED"].unique() - 1
                    return list(b)
                return []

        ScanBase._check_bunit(bunit)
        if bunit.lower() != "ta" and zenith_opacity is None:
            raise ValueError("Can't scale the data without a valid zenith opacity")

        nod_beams = get_nod_beams(self)
        feeds = fdnum
        if fdnum is None:
            logger.info(f"Found nodding beams {nod_beams}")
            feeds = nod_beams
        else:
            if nod_beams != fdnum:
                logger.info(f"Using beams {fdnum} instead of nodding beams {nod_beams} found in the data")
        if type(feeds) is int or len(feeds) != 2:
            raise Exception(f"fdnum={feeds} not valid, need a list with two feeds")
        logger.debug(f"getnod: using fdnum={feeds}")
        prockey = "PROCSEQN"
        procvals = {"ON": 1, "OFF": 2}
        (scans, _sf) = self._common_selection(
            fdnum=feeds,
            ifnum=ifnum,
            plnum=plnum,
            apply_flags=apply_flags,
            prockey=prockey,
            procvals=procvals,
            **kwargs,
        )

        beam1_selected = True
        scanblock = ScanBlock()

        for i in range(len(self._sdf)):
            df0 = select_from("FITSINDEX", i, _sf)
            for f in feeds:
                _df = select_from("FDNUM", f, df0)
                if len(_df) == 0:  # skip IF's and beams not part of the nodding pair.
                    continue
                # scanlist = self._nod_scan_list_selection(scans, _df, feeds, check=False)
                scanlist = self._common_scan_list_selection(scans, _df, prockey=prockey, procvals=procvals, check=False)
                if len(scanlist["ON"]) == 0 or len(scanlist["OFF"]) == 0:
                    logger.debug(f"Some of scans {scans} not found, continuing")
                    continue

                beam1_selected = f == feeds[0]
                logger.debug(f"SCANLIST {scanlist}")
                logger.debug(f"FEED {f} {beam1_selected} {feeds[0]}")
                logger.debug(f"PROCSEQN {set(_df['PROCSEQN'])}")
                logger.debug(f"Sending dataframe with scans {set(_df['SCAN'])}")
                logger.debug(f"and PROC {set(_df['PROC'])}")
                rows = {}
                # Loop over scan pairs.
                c = 0
                for on, off in zip(scanlist["ON"], scanlist["OFF"]):
                    _ondf = select_from("SCAN", on, _df)
                    _offdf = select_from("SCAN", off, _df)
                    rows["ON"] = list(_ondf["ROW"])
                    rows["OFF"] = list(_offdf["ROW"])
                    for key in rows:
                        if len(rows[key]) == 0:
                            raise Exception(f"{key} scans not found in scan list {scans}")
                    # Do not pass scan list here. We need all the cal rows. They will
                    # be intersected with scan rows in NodScan.
                    calrows = {}
                    dfcalT = select_from("CAL", "T", _df)
                    dfcalF = select_from("CAL", "F", _df)
                    calrows["ON"] = list(dfcalT["ROW"])
                    calrows["OFF"] = list(dfcalF["ROW"])
                    d = {"ON": on, "OFF": off}
                    # Check if there is a noise diode.
                    if len(calrows["ON"]) == 0 or nocal:
                        nocal = True
                        if t_sys is None:
                            dfoncalF = select_from("CAL", "F", _ondf)
                            t_sys = dfoncalF["TSYS"].to_numpy()
                            logger.info("Using TSYS column")
                    logger.debug(f"{i, f, c} SCANROWS {rows}")
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
                        bintable=bintable,
                        calibrate=calibrate,
                        smoothref=smoothref,
                        apply_flags=apply_flags,
                        nocal=nocal,
                        tsys=t_sys,
                        bunit=bunit,
                        zenith_opacity=zenith_opacity,
                    )
                    g.merge_commentary(self)
                    scanblock.append(g)
                    c = c + 1
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
        fdnum,
        ifnum,
        plnum,
        calibrate=True,
        fold=True,
        shift_method="fft",
        use_sig=True,
        bintable=None,
        smoothref=1,
        apply_flags=True,
        bunit="ta",
        zenith_opacity=None,
        observer_location=Observatory["GBT"],
        **kwargs,
    ):
        """
        Retrieve and calibrate frequency-switched data.

        Parameters
        ----------
        fdnum: int
            The feed number
        ifnum : int
            The IF number
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
        bintable : int, optional
            Limit to the input binary table index. The default is None which means use all binary tables.
        smooth_ref: int, optional
            the number of channels in the reference to boxcar smooth prior to calibration
        apply_flags : boolean, optional.  If True, apply flags before calibration.
            See :meth:`apply_flags`. Default: True
        bunit : str, optional
            The brightness scale unit for the output scan, must be one of (case-insensitive)
                    - 'ta'  : Antenna Temperature
                    - 'ta*' : Antenna temperature corrected to above the atmosphere
                    - 'jy'  : flux density in Jansky
            If 'ta*' or 'jy' the zenith opacity must also be given. Default:'ta'
        zenith_opacity: float, optional
                The zenith opacity to use in calculating the scale factors for the integrations.  Default:None
        observer_location : `~astropy.coordinates.EarthLocation`
            Location of the observatory. See `~dysh.coordinates.Observatory`.
            This will be transformed to `~astropy.coordinates.ITRS` using the time of
            observation DATE-OBS or MJD-OBS in
            the SDFITS header.  The default is the location of the GBT.
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
        scanblock : `~spectra.scan.ScanBlock`
            ScanBlock containing one or more`~spectra.scan.FSScan`.

        """
        debug = kwargs.pop("debug", False)
        logger.debug(kwargs)

        ScanBase._check_bunit(bunit)
        if bunit.lower() != "ta" and zenith_opacity is None:
            raise ValueError("Can't scale the data without a valid zenith opacity")

        (scans, _sf) = self._common_selection(ifnum=ifnum, plnum=plnum, fdnum=fdnum, apply_flags=apply_flags, **kwargs)

        scanblock = ScanBlock()

        for i in range(len(self._sdf)):
            logger.debug(f"Processing file {i}: {self._sdf[i].filename}")

            df = select_from("FITSINDEX", i, _sf)
            # loop over scans:
            for scan in scans:
                logger.debug(f"doing scan {scan}")
                calrows = {}
                _df = select_from("SCAN", scan, df)
                dfcalT = select_from("CAL", "T", _df)
                dfcalF = select_from("CAL", "F", _df)
                sigrows = {}
                dfsigT = select_from("SIG", "T", _df)
                dfsigF = select_from("SIG", "F", _df)
                #
                calrows["ON"] = list(dfcalT["ROW"])
                calrows["OFF"] = list(dfcalF["ROW"])
                sigrows["ON"] = list(dfsigT["ROW"])
                sigrows["OFF"] = list(dfsigF["ROW"])
                print(f"{scan=} {sigrows=} {calrows=}")
                g = FSScan(
                    self._sdf[i],
                    scan=scan,
                    sigrows=sigrows,
                    calrows=calrows,
                    fdnum=fdnum,
                    ifnum=ifnum,
                    plnum=plnum,
                    bintable=bintable,
                    calibrate=calibrate,
                    fold=fold,
                    shift_method=shift_method,
                    use_sig=use_sig,
                    observer_location=observer_location,
                    smoothref=1,
                    apply_flags=apply_flags,
                    bunit=bunit,
                    zenith_opacity=zenith_opacity,
                    debug=debug,
                )
                g.merge_commentary(self)
                scanblock.append(g)
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
        rx = set(self["FRONTEND"])
        if "Rcvr26_40" not in rx:
            return

        self._fix_column("FDNUM", 1, {"FRONTEND": "Rcvr26_40", "PLNUM": 1})
        logger.info(f"Fixing FDNUM mislabel for Rcvr26_40. FDNUM 0 changed to 1")
        self._fix_column("FDNUM", 0, {"FRONTEND": "Rcvr26_40", "PLNUM": 0})
        logger.info(f"Fixing FDNUM mislabel for Rcvr26_40. FDNUM 1 changed to 0")

    # @todo sig/cal no longer needed?
    @log_call_to_result
    def subbeamnod(
        self,
        fdnum,
        ifnum,
        plnum,
        method="cycle",
        sig=None,
        cal=None,
        calibrate=True,
        timeaverage=True,
        polaverage=False,
        weights="tsys",
        bintable=None,
        smoothref=1,
        apply_flags=True,
        bunit="ta",
        zenith_opacity=None,
        observer_location=Observatory["GBT"],
        **kwargs,
    ):
        """Get a subbeam nod power scan, optionally calibrating it.

        Parameters
        ----------
        fdnum: int
            The feed number
        ifnum : int
            The IF number
        plnum : int
            The polarization number
        method: str
            Method to use when processing. One of 'cycle' or 'scan'.  'cycle' is more accurate and averages data in each SUBREF_STATE cycle. 'scan' reproduces GBTIDL's snodka function which has been shown to be less accurate.  Default:'cycle'
        sig : bool
            True to indicate if this is the signal scan, False if reference
        cal: bool
            True if calibration (diode) is on, False if off.
        calibrate: bool
            whether or not to calibrate the data.  If `True`, the data will be (calon - caloff)*0.5, otherwise it will be SDFITS row data. Default:True
        timeaverage : boolean, optional
            Average the scans in time. The default is True.
        polaverage : boolean, optional
            Average the scans in polarization. The default is False.
        weights: str or None
            None to indicate equal weighting or 'tsys' to indicate tsys weighting to use in time averaging. Default: 'tsys'
        bintable : int, optional
            Limit to the input binary table index. The default is None which means use all binary tables.
        smooth_ref: int, optional
            the number of channels in the reference to boxcar smooth prior to calibration
        apply_flags : boolean, optional.  If True, apply flags before calibration.
            See :meth:`apply_flags`. Default: True
        bunit : str, optional
            The brightness scale unit for the output scan, must be one of (case-insensitive)
                    - 'ta'  : Antenna Temperature
                    - 'ta*' : Antenna temperature corrected to above the atmosphere
                    - 'jy'  : flux density in Jansky
            If 'ta*' or 'jy' the zenith opacity must also be given. Default:'ta'
        zenith_opacity: float, optional
                The zenith opacity to use in calculating the scale factors for the integrations.  Default:None
        observer_location : `~astropy.coordinates.EarthLocation`
            Location of the observatory. See `~dysh.coordinates.Observatory`.
            This will be transformed to `~astropy.coordinates.ITRS` using the time of
            observation DATE-OBS or MJD-OBS in
            the SDFITS header.  The default is the location of the GBT.
        **kwargs : dict
            Optional additional selection keyword arguments, typically
            given as key=value, though a dictionary works too.
            e.g., `ifnum=1, plnum=[2,3]` etc.

        Returns
        -------
        data : `~spectra.scan.ScanBlock`
            A ScanBlock containing one or more `~spectra.scan.SubBeamNodScan`
        """

        ScanBase._check_bunit(bunit)
        if bunit.lower() != "ta" and zenith_opacity is None:
            raise ValueError("Can't scale the data without a valid zenith opacity")

        (scans, _sf) = self._common_selection(ifnum=ifnum, plnum=plnum, fdnum=fdnum, apply_flags=apply_flags, **kwargs)
        scanblock = ScanBlock()

        if method == "cycle":
            # Calibrate each cycle individually and then
            # average the calibrated data.
            for sdfi in range(len(self._sdf)):
                _df = select_from("FITSINDEX", sdfi, _sf)
                for scan in scans:
                    reftp = []
                    sigtp = []
                    fulltp = []
                    logger.debug(f"doing scan {scan}")
                    df = select_from("SCAN", scan, _df)
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
                    # Loop over cycles, calibrating each independently.
                    groups_zip = zip(ref_on_groups, sig_on_groups, ref_off_groups, sig_off_groups)

                    for i, (rgon, sgon, rgoff, sgoff) in enumerate(groups_zip):
                        # Do it the dysh way.
                        calrows = {"ON": rgon, "OFF": rgoff}
                        tprows = np.sort(np.hstack((rgon, rgoff)))
                        reftp.append(
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
                                bintable=bintable,
                                calibrate=calibrate,
                                smoothref=smoothref,
                                apply_flags=apply_flags,
                            )
                        )
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
                                bintable=bintable,
                                calibrate=calibrate,
                                smoothref=smoothref,
                                apply_flags=apply_flags,
                            )
                        )
                    sb = SubBeamNodScan(
                        sigtp,
                        reftp,
                        fdnum=fdnum,
                        ifnum=ifnum,
                        plnum=plnum,
                        calibrate=calibrate,
                        weights=weights,
                        smoothref=smoothref,
                        apply_flags=apply_flags,
                        bunit=bunit,
                        zenith_opacity=zenith_opacity,
                    )
                    scanblock.append(sb)
        elif method == "scan":
            for sdfi in range(len(self._sdf)):
                # Process the whole scan as a single block.
                # This is less accurate, but might be needed if
                # the scan was aborted and there are not enough
                # sig/ref cycles to do a per cycle calibration.
                for scan in scans:
                    reftp = []
                    sigtp = []
                    fulltp = []
                    tpon = self.gettp(
                        fdnum=fdnum,
                        ifnum=ifnum,
                        plnum=plnum,
                        scan=scan,
                        sig=None,
                        cal=None,
                        bintable=bintable,
                        subref=-1,
                        calibrate=calibrate,
                        smoothref=smoothref,
                        apply_flags=apply_flags,
                    )
                    sigtp.append(tpon[0])
                    tpoff = self.gettp(
                        fdnum=fdnum,
                        ifnum=ifnum,
                        plnum=plnum,
                        scan=scan,
                        sig=None,
                        cal=None,
                        bintable=bintable,
                        subref=1,
                        calibrate=calibrate,
                        smoothref=smoothref,
                        apply_flags=apply_flags,
                    )
                    reftp.append(tpoff[0])
                    # in order to reproduce gbtidl tsys, we need to do a normal
                    # total power scan
                    ftp = self.gettp(
                        fdnum=fdnum,
                        ifnum=ifnum,
                        plnum=plnum,
                        scan=scan,
                        sig=None,
                        cal=None,
                        bintable=bintable,
                        calibrate=calibrate,
                        smoothref=smoothref,
                        apply_flags=apply_flags,
                    )
                    fulltp.append(ftp[0])
                sb = SubBeamNodScan(
                    sigtp,
                    reftp,
                    fdnum=fdnum,
                    ifnum=ifnum,
                    plnum=plnum,
                    calibrate=calibrate,
                    weights=weights,
                    smoothref=smoothref,
                    apply_flags=apply_flags,
                    bunit=bunit,
                    zenith_opacity=zenith_opacity,
                )
                sb.merge_commentary(self)
                scanblock.append(sb)
        if len(scanblock) == 0:
            raise Exception("Didn't find any unflagged scans matching the input selection criteria.")
        scanblock.merge_commentary(self)
        return scanblock

    def _common_scan_list_selection(self, scans, selection, prockey, procvals, check=False):
        s = {"ON": [], "OFF": []}
        df2 = selection[selection["SCAN"].isin(scans)]
        procset = set(df2["PROC"])
        lenprocset = len(procset)
        if lenprocset == 0:
            # This is ok since not all files in a set have all the polarizations, feeds, or IFs
            return s
        if lenprocset > 1:
            raise Exception(f"Found more than one PROCTYPE in the requested scans: {procset}")
        proc = list(procset)[0]
        dfon = select_from(prockey, procvals["ON"], selection)
        dfoff = select_from(prockey, procvals["OFF"], selection)
        onscans = uniq(list(dfon["SCAN"]))  # wouldn't set() do this too?
        offscans = uniq(list(dfoff["SCAN"]))
        # pol1 = set(dfon["PLNUM"])
        # pol2 = set(dfoff["PLNUM"])
        # scans = list(selection["SCAN"])
        # The companion scan will always be +/- 1 depending if procseqn is 1(ON) or 2(OFF).
        # First check the requested scan number(s) are in the ONs or OFFs of this bintable.
        seton = set(onscans)
        setoff = set(offscans)
        onrequested = seton.intersection(scans)
        offrequested = setoff.intersection(scans)
        if len(onrequested) == 0 and len(offrequested) == 0:
            raise ValueError(f"Scans {scans} not found in ONs or OFFs")
        # Then check that for each requested ON/OFF there is a matching OFF/ON
        # and build the final matched list of ONs and OFfs.
        sons = list(onrequested.copy())
        soffs = list(offrequested.copy())
        missingoff = []
        missingon = []
        # Special case position switch calibration
        if procvals["ON"] == "PSWITCHON":
            # Figure out the companion scan
            if proc == "OnOff":
                offdelta = 1
                ondelta = -1
            elif proc == "OffOn":
                offdelta = -1
                ondelta = 1
            else:
                raise Exception(f"I don't know how to handle PROCTYPE {proc} for the requested scan operation")
        else:
            # Nod data
            offdelta = 1
            ondelta = -1
        for i in onrequested:
            expectedoff = i + offdelta
            if len(setoff.intersection([expectedoff])) == 0:
                missingoff.append(expectedoff)
            else:
                soffs.append(expectedoff)
        for i in offrequested:
            expectedon = i + ondelta
            if len(seton.intersection([expectedon])) == 0:
                missingon.append(expectedon)
            else:
                sons.append(expectedon)
        if check:
            s["OFF"] = sorted(set(soffs).union(missingoff))
            s["ON"] = sorted(set(sons).union(missingon))
        else:
            if len(missingoff) > 0:
                raise ValueError(
                    f"For the requested ON scans {onrequested}, the OFF scans {missingoff} were not present"
                )
            if len(missingon) > 0:
                raise ValueError(
                    f"For the requested OFF scans {offrequested}, the ON scans {missingon} were not present"
                )
            s["ON"] = sorted(set(sons))
            s["OFF"] = sorted(set(soffs))
            if len(s["ON"]) != len(s["OFF"]):
                raise Exception(f'ON and OFF scan list lengths differ {len(s["ON"])} != {len(s["OFF"])}')
        return s

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
        fi = list(set(_final["FITSINDEX"]))
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
                bintables = list(set(df.BINTABLE))
                for b in bintables:  # loop over the bintables in this fitsfile
                    rows = list(set(df.ROW))
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
                    print(f"Writing {this_rows_written} rows to {outfile}.")
                outhdu.writeto(outfile, output_verify=output_verify, overwrite=overwrite, checksum=checksum)
            if verbose:
                print(f"Total of {total_rows_written} rows written to files.")
        else:
            hdu = self._sdf[fi[0]]._hdu[0].copy()
            outhdu = fits.HDUList(hdu)
            for k in fi:
                df = select_from("FITSINDEX", k, _final)
                bintables = list(set(df.BINTABLE))
                for b in bintables:
                    rows = list(set(df.ROW))
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
            else:
                if verbose:
                    print(f"Writing {total_rows_written} to {fileobj}")
            # outhdu.update_extend()  # possibly unneeded
            outhdu.writeto(fileobj, output_verify=output_verify, overwrite=overwrite, checksum=checksum)
            outhdu.close()

    def _update_radesys(self):
        """
        Updates the 'RADESYS' column of the index for cases when it is empty.
        """

        radesys = {"AzEl": "AltAz", "HADec": "hadec", "Galactic": "galactic"}

        warning_msg = (
            lambda scans, a, coord, limit: f"""Scan(s) {scans} have {a} {coord} below {limit}. The GBT does not go that low. Any operations that rely on the sky coordinates are likely to be inaccurate (e.g., switching velocity frames)."""
        )

        # Elevation below the GBT elevation limit (5 degrees) warning.
        low_el_mask = self["ELEVATIO"] < 5
        if low_el_mask.sum() > 0:
            low_el_scans = map(str, set(self._index.loc[low_el_mask, "SCAN"]))
            warnings.warn(warning_msg(",".join(low_el_scans), "an", "elevation", "5 degrees"))

        # Azimuth and elevation case.
        self._fix_column("RADESYS", radesys["AzEl"], {"CTYPE2": "AZ", "CTYPE3": "EL"})

        # Hour angle and declination case.
        self._fix_column("RADESYS", radesys["HADec"], {"CTYPE2": "HA"})

        # Galactic coordinates.
        self._fix_column("RADESYS", radesys["Galactic"], {"CTYPE2": "GLON"})

    def _fix_column(self, column, new_val, mask_dict):
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
            warnings.warn(f"Changing an existing SDFITS column {items}")
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
            warnings.warn(
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
        flag_rows = np.where(mask == True)[0].tolist()
        self.flag(row=flag_rows)

    def getbeam(self, debug=False):
        """
        Find the two nodding beams based on on a given FDNUM, FEED
        needs PROCSCAN='BEAM1' or 'BEAM2'

        Parameters
        ----------
        sdf : `GBTFITSLoad`
            data handle, containing one or more SDFITS files specific to GBT
        debug : boolean, optional
            Add more debugging output. @todo use logger
            The default is False.

        Returns
        -------
        beams : list of two ints representing the nodding beams (0 = first beam)

        """
        # list of columns needed to differentiate and find the nodding beams
        kb = ["FEEDXOFF", "FEEDEOFF", "PROCSCAN", "FDNUM", "FEED"]
        a0 = self._index[kb]
        b1 = a0.loc[a0["FEEDXOFF"] == 0.0]
        b2 = b1.loc[b1["FEEDEOFF"] == 0.0]
        d1 = b2.loc[b2["PROCSCAN"] == "BEAM1"]
        d2 = b2.loc[b2["PROCSCAN"] == "BEAM2"]
        #
        if len(d1["FDNUM"].unique()) == 1 and len(d2["FDNUM"].unique()) == 1:
            beam1 = d1["FDNUM"].unique()[0]
            beam2 = d2["FDNUM"].unique()[0]
            fdnum1 = d1["FEED"].unique()[0]
            fdnum2 = d2["FEED"].unique()[0]
            if debug:
                print("beams: ", beam1, beam2, fdnum1, fdnum2)
            return [beam1, beam2]
        else:
            # try one other thing
            if len(b2["FEED"].unique()) == 2:
                print("getbeam rescued")
                b = b2["FEED"].unique() - 1
                return list(b)
            print("too many in beam1:", d1["FDNUM"].unique())
            print("too many in beam2:", d2["FDNUM"].unique())
            return []

    def calseq(self, scan, tcold=54, fdnum=0, ifnum=0, plnum=0, freq=None, verbose=False):
        """
        This routine returns the Tsys and gain for the selected W-band channel.

        W-band receivers use a CALSEQ where during a scan three different
        observations are made: sky, cold1 and cold2, from which the
        system temperature is derived.


        Parameters
        ----------
        sdf : `GBTFITSLoad`
            data handle, containing one or more SDFITS files specific to GBT
        scan : int or list of int
            Scan number(s) where CALSEQ is expected. See sdf.summary() to find the scan number(s).
            If multiple scans are used, an average Tsys is computed.
        tcold : float, optional
            Set the cold temperature. See also freq= for an alternative computation.
            The default is 54.
        fdnum : int, optional
            Feed to be used, 0 being the first.
            The default is 0.
        ifnum : int, optional
            IF to be used, 0 being the first.
            The default is 0.
        plnum : int, optional
            Polarization to be used, 0 being the first.
            The default is 0.
        freq : float, optional
            Observing frequency if Tcold to be set different from the default:
            Tcold = 54 - 0.6*(FREQ-77)      FREQ in GHz
            The default is None.
        verbose : boolean, optional
            Add more information mimicking the GBTIDL outout of VANECAL.
            The default is False

        Returns
        -------
        tsys : float
            The system temperature, in K
        g : float
            The gain in K/counts

        """
        if freq is not None:
            # see eq.(13) in GBT memo 302
            tcold = 54 - 0.6 * (freq - 77)
            print(f"Warning: calseq using freq={freq} GHz and setting tcold={tcold} K")

        twarm = self._index["TWARM"].mean()
        # @todo ? there was a period when TWARM was recorded wrongly as 99C, wwhere TAMBIENT (in K) would be better

        tp_args = {"scan": scan, "ifnum": ifnum, "plnum": plnum, "fdnum": fdnum, "calibrate": True, "cal": False}
        vsky = self.gettp(CALPOSITION="Observing", **tp_args).timeaverage()
        vcold1 = self.gettp(CALPOSITION="Cold1", **tp_args).timeaverage()
        vcold2 = self.gettp(CALPOSITION="Cold2", **tp_args).timeaverage()

        if fdnum == 0:
            g = (twarm - tcold) / mean_data(vcold2.data - vcold1.data)
        elif fdnum == 1:
            g = (twarm - tcold) / mean_data(vcold1.data - vcold2.data)
        else:
            print(f"Illegal fdnum={fdnum} for a CALSEQ")
            return None
        tsys = mean_data(g * vsky.data)

        if verbose:
            print(f"Twarm={twarm} Tcold={tcold}")
            print(f"IFNUM {ifnum} PLNUM {plnum} FDNUM {fdnum}")
            print(f"Tsys = {tsys}")
            print(f"Gain [K/counts] = {g}")

        return tsys, g

    # @todo PJT feeds->fdnum and add other standard args
    def vanecal(self, vane_sky, ifnum, plnum, feeds=range(16), mode=2, tcal=None, verbose=False, **kwargs):
        """
        Return Tsys calibration values for all or selected beams of the Argus
        VANE/SKY calibration cycle.


        Parameters
        ----------
        sdf : `GBTFITSLoad`
            data handle, containing one or more SDFITS files specific to GBT
        vane_sky : list of two ints
            The first designates the VANE scan, the second the SKY scan.
            Normally the SKY scan is directly followed by the VANE scan.
            @todo if one scan given, assume sky is vane+1
        ifnum : int
            The IF number
        plnum : int
            The polarization number
        feeds : list of ints, optional
            The default is range(16), i.e. using all Argus beams.
        mode : int, optional
            Mode of computing. See also `mean_tsys()`
            mode=0  Do the mean before the division
            mode=1  Do the mean after the division
            mode=2  Take a median of the inverse division
            The default is 2.
        tcal : float, optional
            Tcal value for normalization. Normally obtained from the
            environment, but offsite cannot be done.
            @todo fix this, but right now it is advised to manually enter tcal.
            The default is None.
        verbose : boolean, optional
            Add more information mimicking the GBTIDL outout of VANECAL.
            The default is False

        Returns
        -------
        tsys : list of floats
            Values of Tsys for each of the `feeds` given.

        """
        vane = vane_sky[0]
        sky = vane_sky[1]
        if len(feeds) == 0:
            print("Warning, no feeds= given")
            return None
        tsys = np.zeros(len(feeds), dtype=float)

        #  for VANE/CAL data usually tcal=1
        if tcal is None:
            tcal = self._index["TCAL"].mean()
            if tcal == 1.0:
                # until we figure this out via getatmos  @todo warn here
                tcal = 100.0  # set to 100K for now

        i = 0
        for f in feeds:
            v = self.gettp(scan=vane, fdnum=f, ifnum=ifnum, plnum=plnum, calibrate=True, cal=False).timeaverage()
            s = self.gettp(scan=sky, fdnum=f, ifnum=ifnum, plnum=plnum, calibrate=True, cal=False).timeaverage()
            if mode == 0:
                mean_off = mean_data(s.data)
                mean_dif = mean_data(v.data - s.data)
                tsys[i] = tcal * mean_off / mean_dif
            elif mode == 1:
                tsys[i] = tcal / mean_data((v.data - s.data) / s.data)
            elif mode == 2:
                tsys[i] = tcal / np.nanmedian((v.data - s.data) / s.data)
            #  vanecal.pro seems to do    tcal / median( (v-s)/s)
            #  as well as not take off the edges
            i = i + 1
        if verbose:
            for i in range(len(feeds)):
                print(f"fdnum,Tsys   {feeds[i]:2d}  {tsys[i]:10.5f}")
            print(f"<Tsys>  {np.nanmean(tsys):.5f} +/- {np.nanstd(tsys):.5f}")
            print(f"mode={mode}")
            print("TCAL=", tcal)

        return tsys

    def _getnod(self, scans, beams, ifnum=0, plnum=0, tsys=None):
        """
        fake getnod() based on alternating gettp() with averaging done internally
        use the real sdf.getnod() for final analysis.
        @todo   this should be replaced by an improved proper getnod()
        sdf:   the sdfits handle
        scans: list of two scans for the nodding
        beams: list of two beams for the nodding
        ifnum: the ifnum to use
        plnum: the plnum to use
        Returns the two nodding spectra, caller is responsible for averaging them, e.g.
             sp1.average(sp2)

        Parameters
        ----------
        sdf : GBTFITSLoad`
            data handle, containing one or more SDFITS files specific to GBT
        scans : list of 2 ints
            list of two scans for the nodding
        beams : list of 2 ints
            list of two beams for the nodding
        ifnum : int, optional
            IF number. The default is 0.
        plnum : int, optional
            Polarization number. The default is 0.
        tsys : float or list of two floats, optional
            Sytem temperature in K. The default is None.

        Returns
        -------
        (sp1, sp2) : tuple of `~spectra.spectrum.Spectrum`
            the two nodding spectra, caller is responsible for averaging them, e.g. `sp1.average(sp2)`
        """

        if tsys is None:
            tsys = np.array([1.0, 1.0])
        if np.isscalar(tsys):
            tsys = np.array([tsys, tsys])
        if len(tsys) == 1:
            tsys = np.array([tsys, tsys])  # because np.isscalar(np.array([1])) is False !

        ps1_on = self.gettp(
            scan=scans[0], fdnum=beams[0], ifnum=ifnum, plnum=plnum, calibrate=True, cal=False
        ).timeaverage()
        ps1_off = self.gettp(
            scan=scans[1], fdnum=beams[0], ifnum=ifnum, plnum=plnum, calibrate=True, cal=False
        ).timeaverage()
        sp1 = (ps1_on - ps1_off) / ps1_off * tsys[0]

        ps2_on = self.gettp(
            scan=scans[1], fdnum=beams[1], ifnum=ifnum, plnum=plnum, calibrate=True, cal=False
        ).timeaverage()
        ps2_off = self.gettp(
            scan=scans[0], fdnum=beams[1], ifnum=ifnum, plnum=plnum, calibrate=True, cal=False
        ).timeaverage()
        sp2 = (ps2_on - ps2_off) / ps2_off * tsys[1]

        sp1.meta["TSYS"] = tsys[0]
        sp2.meta["TSYS"] = tsys[1]

        return (sp1, sp2)


class GBTOffline(GBTFITSLoad):
    """
    GBTOffline('foo')   connects to a GBT project 'foo' using GBTFITSLoad

    Note project directories are assumed to exist in /home/sdfits
    or whereever dysh_data thinks your /home/sdfits lives.

    Also note in GBTIDL one can use SDFITS_DATA instead of DYSH_DATA

    """

    @log_call_to_history
    def __init__(self, fileobj, *args, **kwargs):
        self._offline = fileobj
        self._filename = dysh_data(fileobj)
        GBTFITSLoad.__init__(self, self._filename, *args, **kwargs)


class GBTOnline(GBTFITSLoad):
    """
    GBTOnline('foo')   monitors project 'foo' as if it could be online
    GBTOnline()        monitors for new projects and connects, and refreshes when updated

    Note project directories are assumed to exist in /home/sdfits
    or whereever dysh_data thinks your /home/sdfits lives.

    Also note in GBTIDL one can use SDFITS_DATA instead of DYSH_DATA

    Use dysh_data('?') as a method to get all filenames in SDFITS_DATA

    GBTIDL says:  Connecting to file: .....
                  File has not been updated in xxx.xx minutes.
    """

    @log_call_to_history
    def __init__(self, fileobj=None, *args, **kwargs):
        self._online = fileobj
        self._platform = platform.system()  # cannot update in "Windows":
        # print("GBTOnline not supported on Windows yet, see issue #447")
        if fileobj is not None:
            self._online_mode = 1  # monitor this file
            if os.path.isdir(fileobj):
                GBTFITSLoad.__init__(self, fileobj, *args, **kwargs)
            else:
                self._online = dysh_data(fileobj)
                GBTFITSLoad.__init__(self, self._online, *args, **kwargs)
            print(f"Connecting to explicit file: {self._online} - will be monitoring this")

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

            if not os.path.isdir(sdfits_root):
                print("Cannot find ", sdfits_root)
                return None

            # 1. check the status_file ?
            status_file = "sdfitsStatus.txt"
            if os.path.exists(sdfits_root + "/" + status_file):
                print(f"Warning, found {status_file} but not using it yet")

            # 2. visit each directory where the final leaf contains fits files, and find the most recent one
            n = 0
            mtime_max = 0
            for dirname, subdirs, files in os.walk(sdfits_root):
                # print("dirname",dirname,"subdirs",subdirs)
                if len(subdirs) == 0:
                    n = n + 1
                    # print("===dirname",dirname)
                    for fname in files:
                        if fname.split(".")[-1] == "fits":
                            mtime = os.path.getmtime(dirname + "/" + fname)
                            # print(mtime,fname)
                            if mtime > mtime_max:
                                mtime_max = mtime
                                project = dirname
                            break
            # print(f"Found {n} under {sdfits_root}")
            if n == 0:
                return None

            self._online = project
            GBTFITSLoad.__init__(self, self._online, *args, **kwargs)
            self._mtime = os.path.getmtime(self.filenames()[0])

        # we only test the first filename in the list, assuming they're all being written

        self._mtime = os.path.getmtime(self.filenames()[0])
        # print("MTIME:",self._mtime)
        delta = (time.time() - self._mtime) / 60.0

        print(f"Connected to file: {self._online}")
        print(f"File has not been updated in {delta:.2f} minutes.")
        # end of __init__

    def _reload(self, force=False):
        """force a reload of the latest"""
        if self._platform == "Windows":
            print("warning, cannot reload on Windows, see issue #447")
            return
        if not force:
            mtime = os.path.getmtime(self.filenames()[0])
            if mtime > self._mtime:
                self._mtime = mtime
                print("NEW MTIME:", self._mtime)
                force = True
        if force:
            print(f"Reload {self._online}")
            GBTFITSLoad.__init__(self, self._online)
        return force

    # examples of catchers for reloading

    def summary(self, **kwargs):
        """reload, if need be"""
        self._reload()
        return super().summary(**kwargs)

    def gettp(self, **kwargs):
        self._reload()
        return super().gettp(**kwargs)

    def getps(self, **kwargs):
        self._reload()
        return super().getps(**kwargs)

    def getnod(self, **kwargs):
        self._reload()
        return super().getnod(**kwargs)

    def getfs(self, **kwargs):
        self._reload()
        return super().getfs(**kwargs)

    def subbeamnod(self, **kwargs):
        self._reload()
        return super().subbeamnod(**kwargs)

    def vanecal(self, **kwargs):
        self._reload()
        return super().vanecal(**kwargs)
