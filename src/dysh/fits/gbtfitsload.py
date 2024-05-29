"""Load SDFITS files produced by the Green Bank Telescope"""

import copy
import warnings
from pathlib import Path

# import astropy.units as u
import numpy as np
import pandas as pd
from astropy.io import fits

from ..coordinates import Observatory, decode_veldef
from ..spectra.scan import FSScan, PSScan, ScanBlock, SubBeamNodScan, TPScan
from ..util import consecutive, keycase, select_from, uniq
from ..util.selection import Selection
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


class GBTFITSLoad(SDFITSLoad):
    """
    GBT-specific container to reprensent one or more SDFITS files

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

    def __init__(self, fileobj, source=None, hdu=None, **kwargs):
        kwargs_opts = {
            "fix": False,  # fix non-standard header elements
            "index": True,  # only set to False for performance testing.
            "verbose": False,
        }
        kwargs_opts.update(kwargs)
        path = Path(fileobj)
        self._sdf = []
        self._index = None
        self.GBT = Observatory["GBT"]
        if path.is_file():
            self._sdf.append(SDFITSLoad(fileobj, source, hdu, **kwargs_opts))
        elif path.is_dir():
            # Find all the FITS files in the directory and sort alphabetically
            # because e.g., VEGAS does A,B,C,D,E
            for f in sorted(path.glob("*.fits")):
                if kwargs.get("verbose", None):
                    print(f"doing {f}")
                self._sdf.append(SDFITSLoad(f, source, hdu, **kwargs_opts))
        else:
            raise Exception(f"{fileobj} is not a file or directory path")
        if kwargs_opts["index"]:
            self._create_index_if_needed()

        # We cannot use this to get mmHg as it will disable all default astropy units!
        # https://docs.astropy.org/en/stable/api/astropy.units.cds.enable.html#astropy.units.cds.enable
        # cds.enable()  # to get mmHg
        self._compute_proc()
        if kwargs.get("verbose", None):
            print("==GBTLoad %s" % fileobj)
            self.ushow(0, "OBJECT")
            self.ushow(0, "SCAN")
            self.ushow(0, "SAMPLER")
            self.ushow("PLNUM")
            self.ushow("IFNUM")
            self.ushow(0, "SIG")
            self.ushow(0, "CAL")
            self.ushow(0, "PROCSEQN")
            self.ushow(0, "PROCSIZE")
            self.ushow(0, "OBSMODE")
            self.ushow(0, "SIDEBAND")
        self._selection = Selection(self)
        lsdf = len(self._sdf)
        if lsdf > 1:
            warnings.warn(f"Found {lsdf} FITS files")  # or maybe just print()

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

    def _compute_proc(self):
        """
        Compute the procedure string from obsmode and add to index

        """
        df = self._index["OBSMODE"].str.split(":", expand=True)
        self._index["PROC"] = df[0]
        # Assign these to something that might be useful later,
        # since we have them
        self._index["_OBSTYPE"] = df[1]
        self._index["_SUBOBSMODE"] = df[2]
        for sdf in self._sdf:
            df = sdf._index["OBSMODE"].str.split(":", expand=True)
            sdf._index["PROC"] = df[0]
            sdf._index["_OBSTYPE"] = df[1]
            sdf._index["_SUBOBSMODE"] = df[2]

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
            df = self._index
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
        _df = self._index[show].copy()
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
            # print(f"Uniq data for scan {s}: {nint} {nIF} {nPol} {nfeed} {obj} {proc}")
            s2 = pd.Series(
                [obj, proc, nIF, nPol, nint, nfeed],
                name="uniqued data",
                index=["OBJECT", "PROC", "# IF", "# POL", "# INT", "# FEED"],
            )
            ser = pd.concat([ser, s2]).reindex(comp_colnames)
            ser.rename("appended ser")
            # print("append series data",ser)
            # print("append series index ",ser.index)
            # print("df cols",compressed_df.columns)
            # print("SAME? ",all(ser.index == compressed_df.columns))
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

    def select_track(self, df):
        return df[(df["PROC"] == "Track")]

    # @todo move all selection methods to sdfitsload after adding Selection
    # to sdfitsload
    # @todo write a Delegator class to autopass to Selection. See, e.g., https://michaelcho.me/article/method-delegation-in-python/
    # the problem is I would rather use __getattr__ to allow us to do stuff like sdf["COLUMNAME"] to return a column via _index.
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

    def _create_index_if_needed(self):
        i = 0
        if self._index is None:
            for s in self._sdf:
                if s._index is None:
                    s._create_index()
                # add a FITSINDEX column
                s._index["FITSINDEX"] = i * np.ones(len(s._index))
                if self._index is None:
                    self._index = s._index
                else:
                    self._index = pd.concat([self._index, s._index], axis=0, ignore_index=True)
                i = i + 1

    def info(self):
        """Return information on the HDUs contained in this object. See :meth:`~astropy.HDUList/info()`"""
        for s in self._sdf:
            s.info()

    def getfs(
        self,
        calibrate=True,
        fold=True,
        use_sig=True,
        timeaverage=True,
        polaverage=False,
        weights="tsys",
        bintable=None,
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
        if debug:
            print(kwargs)
        # either the user gave scans on the command line (scans !=None) or pre-selected them
        # with self.selection.selectXX()
        if len(self._selection._selection_rules) > 0:
            _final = self._selection.final
        else:
            _final = self._index
        scans = kwargs.pop("scan", None)
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
        if debug:
            print("scans/w sel:", scans, self._selection)
        fs_selection = copy.deepcopy(self._selection)
        # now downselect with any additional kwargs
        if debug:
            print(f"SELECTION FROM MIXED KWARGS {kwargs}")
        fs_selection._select_from_mixed_kwargs(**kwargs)
        if debug:
            print(fs_selection.show())
        _sf = fs_selection.final
        if len(_sf) == 0:
            raise Exception("Didn't find any scans matching the input selection criteria.")
        # _sf = fs_selection.merge(how='inner')   ## ??? PJT
        ifnum = set(_sf["IFNUM"])
        plnum = set(_sf["PLNUM"])
        scans = set(_sf["SCAN"])
        if debug:
            print(f"using SCANS {scans} IF {ifnum} PL {plnum}")
        scanblock = ScanBlock()

        for i in range(len(self._sdf)):
            df = select_from("FITSINDEX", i, _sf)
            for k in ifnum:
                _ifdf = select_from("IFNUM", k, df)  # one FSScan per ifnum
                if debug:
                    # print(f"SCANLIST {scanlist}")
                    print(f"POLS {set(df['PLNUM'])}")
                    print(f"Sending dataframe with scans {set(_ifdf['SCAN'])}")
                    print(f"and PROC {set(_ifdf['PROC'])}")
                # loop over scans:
                for scan in scans:
                    if debug:
                        print(f"doing scan {scan}")
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
                        use_sig=use_sig,
                        observer_location=observer_location,
                        debug=debug,
                    )
                    scanblock.append(g)
        if len(scanblock) == 0:
            raise Exception("Didn't find any scans matching the input selection criteria.")
        return scanblock
        # end of getfs()

    def getps(self, calibrate=True, timeaverage=True, polaverage=False, weights="tsys", bintable=None, **kwargs):
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
        # print(kwargs)
        # either the user gave scans on the command line (scans !=None) or pre-selected them
        # with select_fromion.selectXX(). In either case make sure the matching ON or OFF
        # is in the starting selection.
        if len(self._selection._selection_rules) > 0:
            _final = self._selection.final
        else:
            _final = self._index
        # print(kwargs)
        scans = kwargs.pop("scan", None)
        debug = kwargs.pop("debug", False)
        kwargs = keycase(kwargs)
        # print(f"case kwargs {kwargs}")
        if type(scans) is int:
            scans = [scans]
        preselected = {}
        for kw in ["SCAN", "IFNUM", "PLNUM"]:
            preselected[kw] = uniq(_final[kw])
        if scans is None:
            scans = preselected["SCAN"]
        missing = self._onoff_scan_list_selection(scans, _final, check=True)
        scans_to_add = set(missing["ON"]).union(missing["OFF"])
        if debug:
            print(f"after check scans_to_add={scans_to_add}")
        # now remove any scans that have been pre-selected by the user.
        # scans_to_add -= scans_preselected
        if debug:
            print(f"after removing preselected {preselected['SCAN']}, scans_to_add={scans_to_add}")
        ps_selection = copy.deepcopy(self._selection)
        if debug:
            print("SCAN ", scans)
            print("TYPE: ", type(ps_selection))
        if len(scans_to_add) != 0:
            # add a rule selecting the missing scans :-)
            if debug:
                print(f"adding rule scan={scans_to_add}")
            kwargs["SCAN"] = list(scans_to_add)
        for k, v in preselected.items():
            if k not in kwargs:
                kwargs[k] = v
        # now downselect with any additional kwargs
        if debug:
            print(f"SELECTION FROM MIXED KWARGS {kwargs}")
            print(ps_selection.show())
        ps_selection._select_from_mixed_kwargs(**kwargs)
        if debug:
            print("AFTER")
            print(ps_selection.show())
        _sf = ps_selection.final
        if len(_sf) == 0:
            raise Exception("Didn't find any scans matching the input selection criteria.")
        ifnum = uniq(_sf["IFNUM"])
        plnum = uniq(_sf["PLNUM"])
        scans = uniq(_sf["SCAN"])
        if debug:
            print(f"FINAL i {ifnum} p {plnum} s {scans}")
        scanblock = ScanBlock()
        for i in range(len(self._sdf)):
            df = select_from("FITSINDEX", i, _sf)
            for k in ifnum:
                _df = select_from("IFNUM", k, df)
                # @todo Calling this method every loop may be expensive. If so, think of
                # a way to tighten it up.
                scanlist = self._onoff_scan_list_selection(scans, _df, check=False)

                if len(scanlist["ON"]) == 0 or len(scanlist["OFF"]) == 0:
                    # print("scans not found, continuing")
                    continue
                if debug:
                    print(f"SCANLIST {scanlist}")
                    print(f"POLS {set(df['PLNUM'])}")
                    print(f"Sending dataframe with scans {set(_df['SCAN'])}")
                    print(f"and PROC {set(_df['PROC'])}")
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
                    # print(f"Sending PSScan({d},ROWS:{rows},CALROWS:{calrows},BT: {bintable}")
                    if debug:
                        print(f"{i, k, c} SCANROWS {rows}")
                        print(f"POL ON {set(_ondf['PLNUM'])} POL OFF {set(_offdf['PLNUM'])}")
                    g = PSScan(
                        self._sdf[i],
                        scans=d,
                        scanrows=rows,
                        calrows=calrows,
                        bintable=bintable,
                        calibrate=calibrate,
                    )
                    scanblock.append(g)
                    c = c + 1
        if len(scanblock) == 0:
            raise Exception("Didn't find any scans matching the input selection criteria.")
        return scanblock

    def gettp(
        self,
        sig=None,
        cal=None,
        calibrate=True,
        timeaverage=True,
        polaverage=False,
        weights="tsys",
        bintable=None,
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
        sigstate = {True: "SIG", False: "REF", None: "BOTH"}
        calstate = {True: "ON", False: "OFF", None: "BOTH"}
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
        ifnum = uniq(_sf["IFNUM"])
        plnum = uniq(_sf["PLNUM"])
        scans = uniq(_sf["SCAN"])
        feeds = uniq(_sf["FDNUM"])
        if debug:
            print(f"FINAL i {ifnum} p {plnum} s {scans} f {feeds}")
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
                    if cal is not None:
                        df = select_from("CAL", TF[cal], df)
                    tprows = list(df["ROW"])
                    if debug:
                        print("TPROWS len=", len(tprows))
                        print("CALROWS on len=", len(calrows["ON"]))
                        print("fitsindex=", i)
                    if len(tprows) == 0:
                        continue
                    g = TPScan(self._sdf[i], scan, sigstate[sig], calstate[cal], tprows, calrows, bintable, calibrate)
                    scanblock.append(g)
        if len(scanblock) == 0:
            raise Exception("Didn't find any scans matching the input selection criteria.")
        return scanblock

    # @todo sig/cal no longer needed?
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
        if debug:
            print(kwargs)

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
        if debug:
            print(f"FINAL i {ifnum} p {plnum} s {scans} f {fdnum}")
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
                        if debug:
                            print(f"doing scan {scan}")
                        df = select_from("SCAN", scan, _ifdf)
                        df_on = df[df["CAL"] == "T"]
                        df_off = df[df["CAL"] == "F"]
                        df_on_sig = df_on[df_on["SUBREF_STATE"] == -1]
                        df_on_ref = df_on[df_on["SUBREF_STATE"] == 1]
                        df_off_sig = df_off[df_off["SUBREF_STATE"] == -1]
                        df_off_ref = df_off[df_off["SUBREF_STATE"] == 1]
                        if debug:
                            print(f"SCANs in df_on_sig {set(df_on_sig['SCAN'])}")
                            print(f"SCANs in df_on_ref {set(df_on_ref['SCAN'])}")
                            print(f"SCANs in df_off_sig {set(df_off_sig['SCAN'])}")
                            print(f"SCANs in df_off_ref {set(df_off_ref['SCAN'])}")
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
                        # print("GROUPS ", ref_on_groups, sig_on_groups, ref_off_groups, sig_off_groups)
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
                                    "BOTH",
                                    "BOTH",
                                    tprows,
                                    calrows,
                                    bintable,
                                    calibrate=calibrate,
                                )
                            )
                            calrows = {"ON": sgon, "OFF": sgoff}
                            tprows = np.sort(np.hstack((sgon, sgoff)))
                            sigtp.append(
                                TPScan(
                                    self._sdf[sdfi],
                                    scan,
                                    "BOTH",
                                    "BOTH",
                                    tprows,
                                    calrows,
                                    bintable,
                                    calibrate=calibrate,
                                )
                            )
                        sb = SubBeamNodScan(sigtp, reftp, method=method, calibrate=calibrate, weights=weights)
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
                            )  # .timeaverage(weights=w)
                            fulltp.append(ftp[0])
                        sb = SubBeamNodScan(sigtp, reftp, fulltp, method=method, calibrate=calibrate, weights=weights)
                        scanblock.append(sb)
        if len(scanblock) == 0:
            raise Exception("Didn't find any scans matching the input selection criteria.")
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
        dfon = select_from("_OBSTYPE", "PSWITCHON", selection)
        dfoff = select_from("_OBSTYPE", "PSWITCHOFF", selection)
        onscans = uniq(list(dfon["SCAN"]))  # wouldn't set() do this too?
        offscans = uniq(list(dfoff["SCAN"]))
        # pol1 = set(dfon["PLNUM"])
        # pol2 = set(dfoff["PLNUM"])
        # print(f"polON {pol1} polOFF {pol2}")
        # scans = list(selection["SCAN"])
        # The companion scan will always be +/- 1 depending if procseqn is 1(ON) or 2(OFF).
        # First check the requested scan number(s) are in the ONs or OFFs of this bintable.
        seton = set(onscans)
        setoff = set(offscans)
        # print(f"SETON {seton} SETOFF {setoff}")
        onrequested = seton.intersection(scans)
        offrequested = setoff.intersection(scans)
        if len(onrequested) == 0 and len(offrequested) == 0:
            raise ValueError(f"Scans {scans} not found in ONs or OFFs")
        # Then check that for each requested ON/OFF there is a matching OFF/ON
        # and build the final matched list of ONs and OFfs.
        sons = list(onrequested.copy())
        soffs = list(offrequested.copy())
        # print(f"SONS {sons} SOFFS {soffs}")
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
            # print(f"DOING ONREQUESTED {i}, looking for off {expectedoff}")
            if len(setoff.intersection([expectedoff])) == 0:
                missingoff.append(expectedoff)
            else:
                soffs.append(expectedoff)
        for i in offrequested:
            expectedon = i + ondelta
            # print(f"DOING OFFEQUESTED {i}, looking for on {expectedon}")
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
        # print(f"onoff_scan_list(scans={scans},if={ifnum},pl={plnum},bintable={bintable},fitsindex={fitsindex})")
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
        dfon = select_from("_OBSTYPE", "PSWITCHON", df)
        dfoff = select_from("_OBSTYPE", "PSWITCHOFF", df)
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
                # print(f"DOING ONREQUESTED {i}, looking for off {expectedoff}")
                if len(setoff.intersection([expectedoff])) == 0:
                    missingoff.append(expectedoff)
                else:
                    soffs.append(expectedoff)
            for i in offrequested:
                expectedon = i + ondelta
                # print(f"DOING OFFEQUESTED {i}, looking for on {expectedon}")
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
        # print(f"onoff_rows(scans={scans},ifnum={ifnum},plnum={plnum},bintable={bintable},fitsindex={fitsindex}")
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
        # print(f"scan_rows(scans={scans},ifnum={ifnum},plnum={plnum},bintable={bintable},fitsindex={fitsindex}")
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
        df_out = []
        rows = []
        scanidx = self._index[self._index["SCAN"].isin(scans)]
        bt = self.udata("BINTABLE")
        for j in bt:
            df = scanidx[scanidx["BINTABLE"] == j]
            rows.append(list(df.index))
        return rows

    def __repr__(self):
        return str(self.files)

    def write_scans(self, fileobj, scans, output_verify="exception", overwrite=False, checksum=False):
        """
        Write specific scans of the `GBTFITSLoad` to a new file.
        TBD: How does this work for multiple files??

        Parameters
        ----------
        fileobj : str, file-like or `pathlib.Path`
            File to write to.  If a file object, must be opened in a
            writeable mode.

        scans: int or list-like
            Range of scans to write out. e.g. 0, [14,25,32].

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
        """
        # get the rows that contain the scans in all bintables
        rows = self._scan_rows_all(scans)
        # @TODO deal with multipl sdfs
        # copy the PrimaryHDU, but not the BinTableHDU
        hdu0 = self._sdf[0]._hdu[0].copy()
        outhdu = fits.HDUList(hdu0)
        # get the bintables rows as new bintables.
        for i in range(len(rows)):
            ob = self._sdf[0]._bintable_from_rows(rows[i], i)
            # print(f"bintable {i} #rows {len(rows[i])} data length {len(ob.data)}")
            if len(ob.data) > 0:
                outhdu.append(ob)
        # print(outhdu.info())
        # write it out!
        outhdu.update_extend()  # possibly unneeded
        outhdu.writeto(fileobj, output_verify=output_verify, overwrite=overwrite, checksum=checksum)
