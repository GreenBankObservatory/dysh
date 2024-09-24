"""Load SDFITS files produced by the Green Bank Telescope"""

import copy
import warnings
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits

from dysh.log import logger

from ..coordinates import Observatory, decode_veldef
from ..log import HistoricalBase, log_call_to_history, log_call_to_result
from ..spectra.scan import FSScan, PSScan, ScanBlock, SubBeamNodScan, TPScan
from ..util import consecutive, indices_where_value_changes, keycase, select_from, uniq
from ..util.selection import Flag, Selection
from .sdfitsload import SDFITSLoad

calibration_kwargs = {
    "calibrate": True,
    "timeaverage": False,
    "polaverage": False,
    "tsys": None,
    "weights": "tsys",
}

# from GBT IDL users guide Table 6.7
# @todo what about the Track/OnOffOn in e.g. AGBT15B_287_33.raw.vegas  (EDGE HI data)
_PROCEDURES = ["Track", "OnOff", "OffOn", "OffOnSameHA", "Nod", "SubBeamNod"]


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

    """

    @log_call_to_history
    def __init__(self, fileobj, source=None, hdu=None, **kwargs):
        kwargs_opts = {
            "fix": False,  # fix non-standard header elements
            "index": True,  # only set to False for performance testing.
            "verbose": False,
        }
        HistoricalBase.__init__(self)
        kwargs_opts.update(kwargs)
        path = Path(fileobj)
        self._sdf = []
        self._selection = None
        self._flag = None
        self.GBT = Observatory["GBT"]
        if path.is_file():
            logger.debug(f"Treating given path {path} as a file")
            self._sdf.append(SDFITSLoad(fileobj, source, hdu, **kwargs_opts))
        elif path.is_dir():
            logger.debug(f"Treating given path {path} as a directory")
            # Find all the FITS files in the directory and sort alphabetically
            # because e.g., VEGAS does A,B,C,D,E
            for f in sorted(path.glob("*.fits")):
                logger.debug(f"Selecting {f} to load")
                if kwargs.get("verbose", None):
                    print(f"doing {f}")
                self._sdf.append(SDFITSLoad(f, source, hdu, **kwargs_opts))
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
            self._create_index_if_needed()
            self._update_radesys()
        # We cannot use this to get mmHg as it will disable all default astropy units!
        # https://docs.astropy.org/en/stable/api/astropy.units.cds.enable.html#astropy.units.cds.enable
        # cds.enable()  # to get mmHg

        if kwargs.get("verbose", None):
            print("==GBTLoad %s" % fileobj)
            self.ushow("OBJECT", 0)
            self.ushow("SCAN", 0)
            self.ushow("SAMPLER", 0)
            self.ushow("PLNUM")
            self.ushow("IFNUM")
            self.ushow("SIG", 0)
            self.ushow("CAL", 0)
            self.ushow("PROCSEQN", 0)
            self.ushow("PROCSIZE", 0)
            self.ushow("OBSMODE", 0)
            self.ushow("SIDEBAND", 0)

        lsdf = len(self._sdf)
        if lsdf > 1:
            print(f"Loaded {lsdf} FITS files")
        self.add_history(f"Project ID: {self.projectID}", add_time=True)

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
        ~pandas.Index
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
        ~dysh.util.Selection
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
        ~pandas.DataFrame
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
        ~dysh.util.Flag
            The Flag object

        """
        return self._flag

    @property
    def final_flags(self):
        """
        The merged flag rules in the Flag object.
        See :meth:`~dysh.util.SelectionBase.final`

        Returns
        -------
        ~pandas.DataFrame
            The final merged flags

        """
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
        index : ~pandas.DataFrame
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
    def rawspectra(self, bintable, fitsindex):
        """
        Get the raw (unprocessed) spectra from the input bintable.

        Parameters
        ----------
        bintable :  int
            The index of the `bintable` attribute
        fitsindex: int
            the index of the FITS file contained in this GBTFITSLoad.  Default:0

        Returns
        -------
        rawspectra : ~numpy.ndarray
            The DATA column of the input bintable

        """
        return self._sdf[fitsindex].rawspectra(bintable)

    def rawspectrum(self, i, bintable=0, fitsindex=0):
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

        Returns
        -------
        rawspectrum : ~numpy.ndarray
            The i-th row of DATA column of the input bintable

        """
        return self._sdf[fitsindex].rawspectrum(i, bintable)

    def getspec(self, i, bintable=0, observer_location=Observatory["GBT"], fitsindex=0):
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
        Returns
        -------
        s : `~dysh.spectra.spectrum.Spectrum`
            The Spectrum object representing the data row.

        """
        return self._sdf[fitsindex].getspec(i, bintable, observer_location)

    def summary(self, scans=None, verbose=False, show_index=True):  # selected=False
        # From GBTIDL:
        # Intended to work with un-calibrated GBT data and is
        # likely to give confusing results for other data.  For other data,
        # list is usually more useful.   @todo what's the dysh eqv. of list ?
        #
        # @todo perhaps return as a astropy.Table then we can have units
        """
        Create a summary list of the input dataset.
        If `verbose=False` (default), some numeric data
        (e.g., RESTFREQ, AZIMUTH, ELEVATIO) are
        averaged over the records with the same scan number.

        Parameters
        ----------
        scans : int or 2-tuple
            The scan(s) to use. A 2-tuple represents (beginning, ending) scans. Default: show all scans

        verbose: bool
            If True, list every record, otherwise return a compact summary
        show_index: bool
            If True, show the DataFrame index column.  Default: False

        Returns
        -------
        summary - `~pandas.DataFrame`
            Summary of the data as a DataFrame.

        """
        # @todo allow user to change show list
        # @todo set individual format options on output by
        # changing these to dicts(?)
        #
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

    # def _select_onoff(self, df):
    #    return df[(df["PROC"] == "OnOff") | (df["PROC"] == "OffOn")]

    # def select_track(self, df):
    #    return df[(df["PROC"] == "Track")]

    # @todo move all selection methods to sdfitsload after adding Selection
    # to sdfitsload
    # @todo write a Delegator class to autopass to Selection. See, e.g., https://michaelcho.me/article/method-delegation-in-python/
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

    @log_call_to_history
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

        See `~dysh.util.selection.Selection`.

        Parameters
        ----------
        chan : number, or array-like
            The channels to select

        Returns
        -------
        None.
        """
        self._selection.select_channel(tag=tag, chan=chan)

    @log_call_to_history
    def flag(self, tag=None, **kwargs):
        """Add one or more exact flag rules, e.g., `key1 = value1, key2 = value2, ...`
        If `value` is array-like then a match to any of the array members will be selected.
        For instance `flag(object=['3C273', 'NGC1234'])` will flag data for either of those
        objects and `flag(ifnum=[0,2])` will flag IF number 0 or IF number 2.
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
        self._selection.select(tag=tag, **kwargs)

    @log_call_to_history
    def flag_range(self, tag=None, **kwargs):
        """
        Flag a range of inclusive values for a given key(s).
        e.g., `key1 = (v1,v2), key2 = (v3,v4), ...`
        will select data  `v1 <= data1 <= v2, v3 <= data2 <= v4, ... `
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
        self._selection.select_range(tag=tag, **kwargs)

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

    @log_call_to_history
    def flag_channel(self, chan, tag=None):
        """
        Select channels and/or channel ranges. These are NOT used in :meth:`final`
        but rather will be used to create a mask for
        flagging. Single arrays/tuples will be treated as channel lists;
        nested arrays will be treated as ranges, for instance

        ``
        # flags channels 1 and 10
        flag_channel([1,10])
        # flagss channels 1 thru 10 inclusive
        flag_channel([[1,10]])
        # flags channel ranges 1 thru 10 and 47 thru 56 inclusive, and channel 75
        flag_channel([[1,10], [47,56], 75)])
        # tuples also work, though can be harder for a human to read
        flag_channel(((1,10), [47,56], 75))
        ``

        See `~dysh.util.selection.Flag`.

        Parameters
        ----------
        chan : number, or array-like
            The channels to flag

        Returns
        -------
        None.
        """
        self._selection.select_channel(tag=tag, chan=chan)

    def _create_index_if_needed(self):
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
        self._index["PROC"] = df[0]
        # Assign these to something that might be useful later,
        # since we have them
        self._index["OBSTYPE"] = df[1]
        self._index["SUBOBSMODE"] = df[2]
        for sdf in self._sdf:
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

        scan_changes = indices_where_value_changes("SCAN", self._index)
        time_changes = indices_where_value_changes("DATE-OBS", self._index)
        # there is probably some super clever pythonic way to do this in one line
        # but I am not that clever, so brute force it.
        intnumarray = []
        for i in self._index.index:
            if i in scan_changes:
                intnum = 0
                # scindex += 1
            else:
                if i in time_changes:
                    intnum += 1
            intnumarray.append(intnum)
        self._index["INTNUM"] = intnumarray

        # Here need to add it as a new column in the BinTableHDU,
        # but we have to sort out FITSINDEX.
        # s.add_col("INTNUM",intnumarray)
        fits_index_changes = indices_where_value_changes("FITSINDEX", self._index)
        lf = len(fits_index_changes)
        for i in range(lf):
            fic = fits_index_changes[i]
            if i + 1 < lf:
                fici = fits_index_changes[i + 1]
            else:
                fici = -1
            fi = self["FITSINDEX"][fic]
            # @todo fix this MWP
            # self._sdf[fi].add_col("INTNUM", intnumarray[fic:fici])  # bintable index???

    def info(self):
        """Return information on the HDUs contained in this object. See :meth:`~astropy.HDUList/info()`"""
        for s in self._sdf:
            s.info()

    @log_call_to_result
    def getfs(
        self,
        calibrate=True,
        fold=True,
        shift_method="fft",
        use_sig=True,
        timeaverage=True,
        polaverage=False,
        weights="tsys",
        bintable=None,
        smoothref=1,
        observer_location=Observatory["GBT"],
        **kwargs,
    ):
        """
        Retrieve and calibrate frequency-switched data.

        Parameters
        ----------
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
        timeaverage : boolean, optional
            Average the scans in time.
            The default is True.
        polaverage : boolean, optional
            Average the scans in polarization.
            The default is False.
        weights : str or None, optional
            How to weight the spectral data when averaging.  'tsys' means use system
            temperature weighting (see e.g., :meth:`~spectra.scan.FSScan.timeaverage`);
            None means uniform weighting.
            The default is 'tsys'.
        bintable : int, optional
            Limit to the input binary table index. The default is None which means use all binary tables.
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
            ScanBlock containing the individual `~spectra.scan.FSScan`s

        """
        debug = kwargs.pop("debug", False)
        logger.debug(kwargs)
        # either the user gave scans on the command line (scans !=None) or pre-selected them
        # with self.selection.selectXX()
        if len(self._selection._selection_rules) > 0:
            _final = self._selection.final
        else:
            _final = self._index
        scans = kwargs.get("scan", None)
        kwargs = keycase(kwargs)
        if type(scans) is int:
            scans = [scans]
        preselected = {}
        for kw in ["SCAN", "IFNUM", "PLNUM", "FDNUM"]:
            preselected[kw] = uniq(_final[kw])
        if scans is None:
            scans = preselected["SCAN"]
        for k, v in preselected.items():
            if k not in kwargs:
                kwargs[k] = v
        logger.debug("scans/w sel:", scans, self._selection)
        fs_selection = copy.deepcopy(self._selection)
        # now downselect with any additional kwargs
        logger.debug(f"SELECTION FROM MIXED KWARGS {kwargs}")
        fs_selection._select_from_mixed_kwargs(**kwargs)
        logger.debug(fs_selection.show())
        _sf = fs_selection.final
        if len(_sf) == 0:
            raise Exception("Didn't find any scans matching the input selection criteria.")
        # _sf = fs_selection.merge(how='inner')   ## ??? PJT
        ifnum = set(_sf["IFNUM"])
        plnum = set(_sf["PLNUM"])
        scans = set(_sf["SCAN"])
        logger.debug(f"using SCANS {scans} IF {ifnum} PL {plnum}")
        scanblock = ScanBlock()

        for i in range(len(self._sdf)):
            df = select_from("FITSINDEX", i, _sf)
            for k in ifnum:
                _ifdf = select_from("IFNUM", k, df)  # one FSScan per ifnum
                logger.debug(f"POLS {set(df['PLNUM'])}")
                logger.debug(f"Sending dataframe with scans {set(_ifdf['SCAN'])}")
                logger.debug(f"and PROC {set(_ifdf['PROC'])}")
                # loop over scans:
                for scan in scans:
                    logger.debug(f"doing scan {scan}")
                    calrows = {}
                    _df = select_from("SCAN", scan, _ifdf)
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
                    g = FSScan(
                        self._sdf[i],
                        scan=scan,
                        sigrows=sigrows,
                        calrows=calrows,
                        bintable=bintable,
                        calibrate=calibrate,
                        fold=fold,
                        shift_method=shift_method,
                        use_sig=use_sig,
                        observer_location=observer_location,
                        smoothref=1,
                        debug=debug,
                    )
                    g.merge_commentary(self)
                    scanblock.append(g)
        if len(scanblock) == 0:
            raise Exception("Didn't find any scans matching the input selection criteria.")
        scanblock.merge_commentary(self)
        return scanblock
        # end of getfs()

    @log_call_to_result
    def getps(
        self, calibrate=True, timeaverage=True, polaverage=False, weights="tsys", bintable=None, smoothref=1, **kwargs
    ):
        """
        Retrieve and calibrate position-switched data.

        Parameters
        ----------
        calibrate : boolean, optional
            Calibrate the scans. The default is True.
        timeaverage : boolean, optional
            Average the scans in time. The default is True.
        polaverage : boolean, optional
            Average the scans in polarization. The default is False.
        weights : str or None, optional
            How to weight the spectral data when averaging.  'tsys' means use system
            temperature weighting (see e.g., :meth:`~spectra.scan.PSScan.timeaverage`);
            None means uniform weighting. The default is 'tsys'.
        bintable : int, optional
            Limit to the input binary table index. The default is None which means use all binary tables.
            (This keyword should eventually go away)
        **kwargs : dict
            Optional additional selection keyword arguments, typically
            given as key=value, though a dictionary works too.
            e.g., `ifnum=1, plnum=[2,3]` etc.

        Raises
        ------
        Exception
            If scans matching the selection criteria are not found.

        Returns
        -------
        scanblock : `~spectra.scan.ScanBlock`
            ScanBlock containing the individual `~spectra.scan.PSScan`s

        """
        # either the user gave scans on the command line (scans !=None) or pre-selected them
        # with select_fromion.selectXX(). In either case make sure the matching ON or OFF
        # is in the starting selection.
        if len(self._selection._selection_rules) > 0:
            _final = self._selection.final
        else:
            _final = self._index
        scans = kwargs.pop("scan", None)
        # debug = kwargs.pop("debug", False)
        kwargs = keycase(kwargs)
        if type(scans) is int:
            scans = [scans]
        preselected = {}
        for kw in ["SCAN", "IFNUM", "PLNUM"]:
            preselected[kw] = uniq(_final[kw])
        if scans is None:
            scans = preselected["SCAN"]
        missing = self._onoff_scan_list_selection(scans, _final, check=True)
        scans_to_add = set(missing["ON"]).union(missing["OFF"])
        logger.debug(f"after check scans_to_add={scans_to_add}")
        # now remove any scans that have been pre-selected by the user.
        # scans_to_add -= scans_preselected
        logger.debug(f"after removing preselected {preselected['SCAN']}, scans_to_add={scans_to_add}")
        ps_selection = copy.deepcopy(self._selection)
        logger.debug(f"SCAN {scans}")
        logger.debug(f"TYPE {type(ps_selection)}")
        if len(scans_to_add) != 0:
            # add a rule selecting the missing scans :-)
            logger.debug(f"adding rule scan={scans_to_add}")
            kwargs["SCAN"] = list(scans_to_add)
        for k, v in preselected.items():
            if k not in kwargs:
                kwargs[k] = v
        # now downselect with any additional kwargs
        ps_selection._select_from_mixed_kwargs(**kwargs)
        _sf = ps_selection.final
        logger.debug(f"{_sf = }")
        if len(_sf) == 0:
            raise Exception("Didn't find any scans matching the input selection criteria.")
        ifnum = uniq(_sf["IFNUM"])
        plnum = uniq(_sf["PLNUM"])
        scans = uniq(_sf["SCAN"])
        logger.debug(f"FINAL i {ifnum} p {plnum} s {scans}")
        scanblock = ScanBlock()
        for i in range(len(self._sdf)):
            df = select_from("FITSINDEX", i, _sf)
            for k in ifnum:
                _df = select_from("IFNUM", k, df)
                # @todo Calling this method every loop may be expensive. If so, think of
                # a way to tighten it up.
                scanlist = self._onoff_scan_list_selection(scans, _df, check=False)

                if len(scanlist["ON"]) == 0 or len(scanlist["OFF"]) == 0:
                    logger.debug("scans not found, continuing")
                    continue
                logger.debug(f"SCANLIST {scanlist}")
                logger.debug(f"POLS {set(df['PLNUM'])}")
                logger.debug(f"Sending dataframe with scans {set(_df['SCAN'])}")
                logger.debug(f"and PROC {set(_df['PROC'])}")
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
                    logger.debug(f"{i, k, c} SCANROWS {rows}")
                    logger.debug(f"POL ON {set(_ondf['PLNUM'])} POL OFF {set(_offdf['PLNUM'])}")
                    g = PSScan(
                        self._sdf[i],
                        scans=d,
                        scanrows=rows,
                        calrows=calrows,
                        bintable=bintable,
                        calibrate=calibrate,
                        smoothref=smoothref,
                    )
                    g.merge_commentary(self)
                    scanblock.append(g)
                    c = c + 1
        if len(scanblock) == 0:
            raise Exception("Didn't find any scans matching the input selection criteria.")
        scanblock.merge_commentary(self)
        return scanblock

    @log_call_to_result
    def gettp(
        self,
        sig=None,
        cal=None,
        calibrate=True,
        timeaverage=True,
        polaverage=False,
        weights="tsys",
        bintable=None,
        smoothref=1,
        **kwargs,
    ):
        """
        Get a total power scan, optionally calibrating it.

        Parameters
        ----------
        sig : bool or None
            True to use only integrations where signal state is True, False to use reference state (signal state is False). None to use all integrations.
        cal: bool or None
            True to use only integrations where calibration (diode) is on, False if off. None to use all integrations regardless calibration state.
            The system temperature will be calculated from both states regardless of the value of this variable.
        calibrate: bool
            whether or not to calibrate the data.  If `True`, the data will be (calon - caloff)*0.5, otherwise it will be SDFITS row data. Default:True
        timeaverage : boolean, optional
            Average the scans in time. The default is True.
        polaverage : boolean, optional
            Average the scans in polarization. The default is False.
        weights: str or None
            None or 'tsys' to indicate equal weighting or tsys weighting to use in time averaging. Default: 'tsys'
        bintable : int, optional
            Limit to the input binary table index. The default is None which means use all binary tables.
        **kwargs : dict
            Optional additional selection  keyword arguments, typically
            given as key=value, though a dictionary works too.
            e.g., `ifnum=1, plnum=[2,3]` etc.

        Returns
        -------
        data : `~spectra.scan.ScanBlock`
            A ScanBlock containing one or more `~spectra.scan.TPScan`

        """
        TF = {True: "T", False: "F"}

        if len(self._selection._selection_rules) > 0:
            _final = self._selection.final
        else:
            _final = self._index
        scans = kwargs.get("scan", None)
        debug = kwargs.pop("debug", False)
        kwargs = keycase(kwargs)
        if type(scans) is int:
            scans = [scans]
        preselected = {}
        for kw in ["SCAN", "IFNUM", "PLNUM"]:
            preselected[kw] = uniq(_final[kw])
        if scans is None:
            scans = preselected["SCAN"]
        ps_selection = copy.deepcopy(self._selection)
        for k, v in preselected.items():
            if k not in kwargs:
                kwargs[k] = v
        # now downselect with any additional kwargs
        ps_selection._select_from_mixed_kwargs(**kwargs)
        _sf = ps_selection.final
        logger.debug("SF=", _sf)
        ifnum = uniq(_sf["IFNUM"])
        plnum = uniq(_sf["PLNUM"])
        scans = uniq(_sf["SCAN"])
        feeds = uniq(_sf["FDNUM"])
        logger.debug(f"FINAL i {ifnum} p {plnum} s {scans} f {feeds}")
        scanblock = ScanBlock()
        calrows = {}
        # @todo loop over feeds too?
        for i in range(len(self._sdf)):
            df = select_from("FITSINDEX", i, _sf)
            for k in ifnum:
                _ifdf = select_from("IFNUM", k, df)
                for scan in scans:
                    df = select_from("SCAN", scan, _ifdf)
                    dfcalT = select_from("CAL", "T", df)
                    dfcalF = select_from("CAL", "F", df)
                    calrows["ON"] = list(dfcalT["ROW"])
                    calrows["OFF"] = list(dfcalF["ROW"])
                    if len(calrows["ON"]) != len(calrows["OFF"]):
                        raise Exception(f'unbalanced calrows {len(calrows["ON"])} != {len(calrows["OFF"])}')
                    # sig and cal are treated specially since
                    # they are not in kwargs and in SDFITS header
                    # they are not booleans but chars
                    if sig is not None:
                        df = select_from("SIG", TF[sig], df)
                    # if cal is not None:
                    #    df = select_from("CAL", TF[cal], df)
                    # the rows with the selected sig state and all cal states
                    tprows = list(df["ROW"])
                    logger.debug("TPROWS len=", len(tprows))
                    logger.debug("CALROWS on len=", len(calrows["ON"]))
                    logger.debug("fitsindex=", i)
                    if len(tprows) == 0:
                        continue
                    g = TPScan(
                        self._sdf[i],
                        scan,
                        sig,
                        cal,
                        tprows,
                        calrows,
                        bintable,
                        calibrate,
                        smoothref=smoothref,
                    )
                    g.merge_commentary(self)
                    scanblock.append(g)
        if len(scanblock) == 0:
            raise Exception("Didn't find any scans matching the input selection criteria.")
        scanblock.merge_commentary(self)
        return scanblock

    # @todo sig/cal no longer needed?
    @log_call_to_result
    def subbeamnod(
        self,
        method="cycle",
        sig=None,
        cal=None,
        calibrate=True,
        timeaverage=True,
        polaverage=False,
        weights="tsys",
        bintable=None,
        smoothref=1,
        **kwargs,
    ):
        """Get a subbeam nod power scan, optionally calibrating it.

        Parameters
        ----------
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
        **kwargs : dict
            Optional additional selection keyword arguments, typically
            given as key=value, though a dictionary works too.
            e.g., `ifnum=1, plnum=[2,3]` etc.

        Returns
        -------
        data : `~spectra.scan.ScanBlock`
            A ScanBlock containing one or more `~spectra.scan.SubBeamNodScan`
        """
        if len(self._selection._selection_rules) > 0:
            _final = self._selection.final
        else:
            _final = self._index
        scans = kwargs.get("scan", None)
        debug = kwargs.pop("debug", False)
        kwargs = keycase(kwargs)
        logger.debug(kwargs)

        if type(scans) is int:
            scans = [scans]
        preselected = {}
        for kw in ["SCAN", "IFNUM", "PLNUM", "FDNUM"]:
            preselected[kw] = uniq(_final[kw])
        if scans is None:
            scans = preselected["SCAN"]
        for k, v in preselected.items():
            if k not in kwargs:
                kwargs[k] = v
        # Check if we are dealing with Ka data before the beam switch.
        rx = np.unique(_final["FRONTEND"])
        if len(rx) > 1:
            raise TypeError("More than one receiver for the selected scan.")
        elif rx[0] == "Rcvr26_40":  # and df["DATE-OBS"][-1] < xxxx
            # Switch the polarizations to match the beams,
            # for this receiver only because it has had its feeds
            # mislabelled since $DATE.
            # For the rest of the receivers the method should use
            # the same polarization for the selected feeds.
            # See also issue #160
            # NOTE THIS "FIX" FAILS if kwargs["FDNUM"] has multiple values
            # e.g.  kwargs["FDNUM"]=[0,1]
            if kwargs["FDNUM"] == 0:
                kwargs["PLNUM"] = 1
            elif kwargs["FDNUM"] == 1:
                kwargs["PLNUM"] = 0
        # now downselect with any additional kwargs
        ps_selection = copy.deepcopy(self._selection)
        ps_selection._select_from_mixed_kwargs(**kwargs)
        _sf = ps_selection.final
        ifnum = uniq(_sf["IFNUM"])
        plnum = uniq(_sf["PLNUM"])
        scans = uniq(_sf["SCAN"])
        fdnum = uniq(_sf["FDNUM"])
        logger.debug(f"FINAL i {ifnum} p {plnum} s {scans} f {fdnum}")
        scanblock = ScanBlock()

        if method == "cycle":
            # Calibrate each cycle individually and then
            # average the calibrated data.
            for sdfi in range(len(self._sdf)):
                _df = select_from("FITSINDEX", sdfi, _sf)
                for k in ifnum:
                    # Row selection.
                    _ifdf = select_from("IFNUM", k, _df)
                    for scan in scans:
                        reftp = []
                        sigtp = []
                        fulltp = []
                        logger.debug(f"doing scan {scan}")
                        df = select_from("SCAN", scan, _ifdf)
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
                                    bintable,
                                    calibrate=calibrate,
                                    smoothref=smoothref,
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
                                    bintable,
                                    calibrate=calibrate,
                                    smoothref=smoothref,
                                )
                            )
                        sb = SubBeamNodScan(sigtp, reftp, calibrate=calibrate, weights=weights, smoothref=smoothref)
                        scanblock.append(sb)
        elif method == "scan":
            for sdfi in range(len(self._sdf)):
                # Process the whole scan as a single block.
                # This is less accurate, but might be needed if
                # the scan was aborted and there are not enough
                # sig/ref cycles to do a per cycle calibration.
                for k in ifnum:
                    for fn in fdnum:
                        for scan in scans:
                            reftp = []
                            sigtp = []
                            fulltp = []
                            tpon = self.gettp(
                                scan=scan,
                                sig=None,
                                cal=None,
                                bintable=bintable,
                                fdnum=fn,
                                plnum=plnum,
                                ifnum=k,
                                subref=-1,
                                weights=weights,
                                calibrate=calibrate,
                                smoothref=smoothref,
                            )
                            sigtp.append(tpon[0])
                            tpoff = self.gettp(
                                scan=scan,
                                sig=None,
                                cal=None,
                                bintable=bintable,
                                fdnum=fn,
                                plnum=plnum,
                                ifnum=k,
                                subref=1,
                                weights=weights,
                                calibrate=calibrate,
                                smoothref=smoothref,
                            )
                            reftp.append(tpoff[0])
                            # in order to reproduce gbtidl tsys, we need to do a normal
                            # total power scan
                            ftp = self.gettp(
                                scan=scan,
                                sig=None,
                                cal=None,
                                bintable=bintable,
                                fdnum=fn,
                                plnum=plnum,
                                ifnum=k,
                                weights=weights,
                                calibrate=calibrate,
                                smoothref=smoothref,
                            )  # .timeaverage(weights=w)
                            fulltp.append(ftp[0])
                        sb = SubBeamNodScan(
                            sigtp,
                            reftp,
                            calibrate=calibrate,
                            weights=weights,
                            smoothref=smoothref,
                        )
                        sb.merge_commentary(self)
                        scanblock.append(sb)
        if len(scanblock) == 0:
            raise Exception("Didn't find any scans matching the input selection criteria.")
        scanblock.merge_commentary(self)
        return scanblock

    def _onoff_scan_list_selection(self, scans, selection, check=False):
        """
        Get the scans for position-switch data sorted
        by ON and OFF state using the current selection

        Parameters
        ----------
        scans : array-like
            list of one or more scans

        selection : `~pandas.DataFrame`
            selection object

        check : boolean
            If True, when scans are mising, return the missing scans in the ON, OFF dict.
            If False, return the normal scanlist and except if scans are missing

        Returns
        -------
        rows : dict
            A dictionary with keys 'ON' and 'OFF' giving the scan numbers of ON and OFF data for the input scan(s)
        """
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
        dfon = select_from("OBSTYPE", "PSWITCHON", selection)
        dfoff = select_from("OBSTYPE", "PSWITCHOFF", selection)
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
        # Figure out the companion scan
        if proc == "OnOff":
            offdelta = 1
            ondelta = -1
        elif proc == "OffOn":
            offdelta = -1
            ondelta = 1
        else:
            raise Exception(
                f"I don't know how to handle PROCTYPE {self._selection.final['PROC']} for the requested scan operation"
            )
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
                raise Exception('ON and OFF scan list lengths differ {len(s["ON"])} != {len(s["OFF"]}')
        return s

    def onoff_scan_list(self, scans=None, ifnum=0, plnum=0, bintable=None, fitsindex=0):
        # need to change to selection kwargs, allow fdnum etc and allow these values to be None
        """Get the scans for position-switch data sorted
           by ON and OFF state

        Parameters
        ----------
        scans : int or list-like
            The scan numbers to find the rows of
        ifnum : int
            the IF index
        plnum : int
            the polarization index


        Returns
        -------
        rows : dict
            A dictionary with keys 'ON' and 'OFF' giving the scan numbers of ON and OFF data for the input scan(s)
        """
        self._create_index_if_needed()
        s = {"ON": [], "OFF": []}
        if type(scans) == int:
            scans = [scans]
        df = self.index(bintable=bintable, fitsindex=fitsindex)
        if plnum is not None:
            df = df[df["PLNUM"] == plnum]
        if ifnum is not None:
            df = df[df["IFNUM"] == ifnum]
        # don't want to limit scans yet since only on or off scan scan numbers may have been
        # passed in, but do need to ensure that a single PROCTYPE is in the given scans
        # Alterative is to this check at the end (which is probably better)
        df2 = df[df["SCAN"].isin(scans)]
        procset = set(uniq(df2["PROC"]))
        lenprocset = len(procset)
        if lenprocset == 0:
            # This is ok since not all files in a set have all the polarizations, feeds, or IFs
            return s
        if lenprocset > 1:
            raise Exception(f"Found more than one PROCTYPE in the requested scans: {procset}")
        proc = list(procset)[0]
        dfon = select_from("OBSTYPE", "PSWITCHON", df)
        dfoff = select_from("OBSTYPE", "PSWITCHOFF", df)
        onscans = uniq(list(dfon["SCAN"]))  # wouldn't set() do this too?
        offscans = uniq(list(dfoff["SCAN"]))
        if scans is not None:
            # The companion scan will always be +/- 1 depending if procseqn is 1(ON) or 2(OFF).
            # First check the requested scan number(s) are in the ONs or OFFs of this bintable.
            seton = set(onscans)
            setoff = set(offscans)
            onrequested = seton.intersection(scans)
            offrequested = setoff.intersection(scans)
            if len(onrequested) == 0 and len(offrequested) == 0:
                raise ValueError(f"Scans {scans} not found in ONs or OFFs of bintable {bintable}")
            # Then check that for each requested ON/OFF there is a matching OFF/ON
            # and build the final matched list of ONs and OFfs.
            sons = list(onrequested.copy())
            soffs = list(offrequested.copy())
            missingoff = []
            missingon = []
            # Figure out the companion scan
            if proc == "OnOff":
                offdelta = 1
                ondelta = -1
            elif proc == "OffOn":
                offdelta = -1
                ondelta = 1
            else:
                raise Exception(f"I don't know how to handle PROCTYPE {df['PROC']} for the requested scan operation")
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
            if len(missingoff) > 0:
                raise ValueError(
                    f"For the requested ON scans {onrequested}, the OFF scans {missingoff} were not present in bintable"
                    f" {bintable}"
                )
            if len(missingon) > 0:
                raise ValueError(
                    f"For the requested OFF scans {offrequested}, the ON scans {missingon} were not present in bintable"
                    f" {bintable}"
                )
            s["ON"] = sorted(set(sons))
            s["OFF"] = sorted(set(soffs))
        else:
            s["ON"] = uniq(list(dfon["SCAN"]))
            s["OFF"] = uniq(list(dfoff["SCAN"]))

        return s

    def calonoff_rows(self, scans=None, bintable=None, fitsindex=0, **kwargs):
        """
        Get individual scan row numbers  sorted by whether the calibration (diode) was on or off, and selected by ifnum,plnum, fdnum,subref,bintable.

        Parameters
        ----------
        scans : int or list-like
            The scan numbers to find the rows of
        ifnum : int
            the IF index
        plnum : int
            the polarization index
        fdnum : int
            the feed index
        subref : int
            the subreflector state (-1,0,1)
        bintable : int
            the index for BINTABLE containing the scans
        fitsindex: int
            the index of the FITS file contained in this GBTFITSLoad.  Default:0

        Returns
        -------
        rows : dict
            A dictionary with keys 'ON' and 'OFF' giving the row indices of CALON and CALOFF data for the input scan(s)

        """
        self._create_index_if_needed()
        s = {"ON": [], "OFF": []}
        ifnum = kwargs.get("ifnum", None)
        plnum = kwargs.get("plnum", None)
        fdnum = kwargs.get("fdnum", None)
        subref = kwargs.get("subref", None)
        if type(scans) == int:
            scans = [scans]
        df = self.index(bintable=bintable, fitsindex=fitsindex)
        if scans is not None:
            df = df[df["SCAN"].isin(scans)]
        dfon = select_from("CAL", "T", df)
        dfoff = select_from("CAL", "F", df)
        if ifnum is not None:
            dfon = select_from("IFNUM", ifnum, dfon)
            dfoff = select_from("IFNUM", ifnum, dfoff)
        if plnum is not None:
            dfon = select_from("PLNUM", plnum, dfon)
            dfoff = select_from("PLNUM", plnum, dfoff)
        if fdnum is not None:
            dfon = select_from("FDNUM", fdnum, dfon)
            dfoff = select_from("FDNUM", fdnum, dfoff)
        if subref is not None:
            dfon = select_from("SUBREF_STATE", subref, dfon)
            dfoff = select_from("SUBREF_STATE", subref, dfoff)
        s["ON"] = list(dfon.index)
        s["OFF"] = list(dfoff.index)
        return s

    # def _onoff_rows_selection(self, scanlist):
    #    """
    #    Get individual ON/OFF (position switch) scan row numbers selected by ifnum,plnum, bintable.
    #
    #   Parameters
    #    scanlist : dict
    #        dictionary of ON and OFF scans
    #    bintable : int
    #        the index for BINTABLE in `sdfits` containing the scans. Default:None
    #    fitsindex: int
    #         the index of the FITS file contained in this GBTFITSLoad.  Default:0
    #
    #     Returns
    #     -------
    #     rows : dict
    #         A dictionary with keys 'ON' and 'OFF' giving the row indices of the ON and OFF data for the input scan(s)

    #    """
    #      rows = {"ON": [], "OFF": []}
    #     # scans is now a dict of "ON" "OFF
    #     for key in scanlist:
    #         rows[key] = self.scan_rows(scanlist[key], ifnum, plnum, bintable, fitsindex=fitsindex)
    #     return rows

    def onoff_rows(self, scans=None, ifnum=0, plnum=0, bintable=None, fitsindex=0):
        """
        Get individual ON/OFF (position switch) scan row numbers selected by ifnum,plnum, bintable.

        Parameters
        ----------
        scans : int or list-like
            The scan numbers to find the rows of
        ifnum : int
            the IF index
        plnum : int
            the polarization index
        bintable : int
            the index for BINTABLE in `sdfits` containing the scans
        fitsindex: int
            the index of the FITS file contained in this GBTFITSLoad.  Default:0

        Returns
        -------
        rows : dict
            A dictionary with keys 'ON' and 'OFF' giving the row indices of the ON and OFF data for the input scan(s)

        """
        # @TODO deal with mulitple bintables
        # @TODO rename this sigref_rows?
        # keep the bintable keyword and allow iteration over bintables if requested (bintable=None)
        rows = {"ON": [], "OFF": []}
        if type(scans) is int:
            scans = [scans]
        _scans = self.onoff_scan_list(scans, ifnum, plnum, bintable, fitsindex=fitsindex)
        # scans is now a dict of "ON" "OFF
        for key in _scans:
            rows[key] = self.scan_rows(_scans[key], ifnum, plnum, bintable, fitsindex=fitsindex)
        return rows

    def scan_rows(self, scans, ifnum=0, plnum=0, bintable=None, fitsindex=0):
        """
        Get scan rows selected by ifnum,plnum, bintable.

        Parameters
        ----------
        scans : int or list-like
            The scan numbers to find the rows of
        ifnum : int
            the IF index
        plnum : int
            the polarization index
        bintable : int
            the index for BINTABLE in `sdfits` containing the scans. Default:None
        fitsindex: int
            the index of the FITS file contained in this GBTFITSLoad.  Default:0

        Returns
        -------
        rows : list
            Lists of the rows in each bintable that contain the scans. Index of `rows` is the bintable index number

        """
        # scans is a list
        self._create_index_if_needed()
        if scans is None:
            raise ValueError("Parameter 'scans' cannot be None. It must be int or list of int")
        df = self.index(bintable=bintable, fitsindex=fitsindex)
        df = df[df["SCAN"].isin(scans)]
        if plnum is not None:
            df = df[df["PLNUM"] == plnum]
        if ifnum is not None:
            df = df[df["IFNUM"] == ifnum]
        rows = list(df["ROW"])
        if len(rows) == 0:
            raise Exception(f"Scans {scans} not found in bintable {bintable}")
        return rows

    def _scan_rows_all(self, scans):
        """
        Get scan rows regardless of ifnum,plnum, bintable.

        Parameters
        ----------
        scans : int or list-like
            The scan numbers to find the rows of

        Returns
        -------
        rows : list
            Lists of the rows in each bintable that contain the scans. Index of `rows` is the bintable index number

        """
        if scans is None:
            raise ValueError("Parameter 'scans' cannot be None. It must be int or list of int")
        rows = []
        scanidx = self[self["SCAN"].isin(scans)]
        bt = self.udata("BINTABLE")
        for j in bt:
            df = scanidx[scanidx["BINTABLE"] == j]
            rows.append(list(df.index))
        return rows

    def write(
        self,
        fileobj,
        multifile=True,
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
        debug = kwargs.pop("debug", False)
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
                for b in bintables:
                    rows = list(set(df.ROW))
                    lr = len(rows)
                    if lr > 0:
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
            hdu = self._sdf[fi[0]]._hdu[fi[0]].copy()
            outhdu = fits.HDUList(hdu)
            for k in fi:
                df = select_from("FITSINDEX", k, _final)
                bintables = list(set(df.BINTABLE))
                for b in bintables:
                    rows = list(set(df.ROW))
                    lr = len(rows)
                    if lr > 0:
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

        radesys = {"AzEl": "AltAz", "HADec": "hadec"}

        warning_msg = (
            lambda scans, a, coord, limit: f"""Scan(s) {scans} have {a} {coord} below {limit}. The GBT does not go that low. Any operations that rely on the sky coordinates are likely to be inaccurate (e.g., switching velocity frames)."""
        )

        # Elevation below the GBT elevation limit (5 degrees) warning.
        low_el_mask = self["ELEVATIO"] < 5
        if low_el_mask.sum() > 0:
            low_el_scans = map(str, set(self._index.loc[low_el_mask, "SCAN"]))
            warnings.warn(warning_msg(",".join(low_el_scans), "an", "elevation", "5 degrees"))

        # Azimuth and elevation case.
        azel_mask = (self["CTYPE2"] == "AZ") & (self["CTYPE3"] == "EL")
        # Update self._index.
        self._index.loc[azel_mask, "RADESYS"] = radesys["AzEl"]
        # Update SDFITSLoad.index.
        sdf_idx = set(self["FITSINDEX"][azel_mask])
        for i in sdf_idx:
            sdfi = self._sdf[i].index()
            azel_mask = (sdfi["CTYPE2"] == "AZ") & (sdfi["CTYPE3"] == "EL")
            sdfi.loc[azel_mask, "RADESYS"] = radesys["AzEl"]

        # Hour angle and declination case.
        hadec_mask = self["CTYPE2"] == "HA"
        self._index.loc[hadec_mask, "RADESYS"] = radesys["HADec"]
        sdf_idx = set(self["FITSINDEX"][hadec_mask])
        for i in sdf_idx:
            sdfi = self._sdf[i].index()
            hadec_mask = sdfi["CTYPE2"] == "HA"
            sdfi.loc[hadec_mask, "RADESYS"] = radesys["HADec"]

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
            warnings.warn("Changing an existing SDFITS column")
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
