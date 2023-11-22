"""Load SDFITS files produced by the Green Bank Telescope"""

import copy
import os
import sys
import warnings
from pathlib import Path

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.units import cds
from astropy.wcs import WCS

from ..spectra.core import tsys_weight
from ..spectra.scan import PSScan, ScanBlock, SubBeamNodScan, TPScan
from ..spectra.spectrum import Spectrum
from ..util import consecutive, uniq
from .sdfitsload import SDFITSLoad

# from GBT IDL users guide Table 6.7
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
            "index": True,
            "verbose": False,
        }
        kwargs_opts.update(kwargs)
        path = Path(fileobj)
        self._sdf = []
        self._index = None
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

        cds.enable()  # to get mmHg
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

    def files(self):
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

    def summary(self, scans=None, verbose=False, show_index=True):  # selected=False
        # From GBTIDL:
        # Intended to work with un-calibrated GBT data and is
        # likely to give confusing results for other data.  For other data,
        # list is usually more useful.
        #
        # @TODO perhaps return as a astropy.Table then we can have units
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
            _df = self.select_scans(scans, _df).filter(show)
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
            uf = self.select("SCAN", s, uncompressed_df)
            # for some columns we will display
            # the mean value
            ser = uf.filter(avg_cols).mean(numeric_only=True)
            ser.rename("filtered ser")
            # for others we will count how many there are
            nIF = uf["IFNUM"].nunique()
            nPol = uf["PLNUM"].nunique()
            nfeed = uf["FEED"].nunique()
            nint = len(set(uf["DATE-OBS"]))  # see gbtidl io/line_index__define.pro
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

    def velocity_convention(self, veldef, velframe):
        # GBT uses VELDEF and VELFRAME incorrectly.
        return "doppler_radio"

    def select_scans(self, scans, df):
        return df[(df["SCAN"] >= scans[0]) & (df["SCAN"] <= scans[1])]

    def select_onoff(self, df):
        return df[(df["PROC"] == "OnOff") | (df["PROC"] == "OffOn")]

    def select(self, key, value, df):
        """
        Select data where key=value.

        Parameters
        ----------
        key : str
            The key value (SDFITS column name)
        value : any
            The value to match
        df : `~pandas.DataFrame`
            The DataFrame to search

        Returns
        -------
        df : `~pandas.DataFrame`
            The subselected DataFrame

        """
        return df[(df[key] == value)]

    def _create_index_if_needed(self):
        if self._index is None:
            for s in self._sdf:
                if s._index is None:
                    s._create_index()
                if self._index is None:
                    self._index = s._index
                else:
                    self._index = pd.concat([self._index, s._index], axis=0, ignore_index=True)

    def info(self):
        """Return the `~astropy.HDUList` info()"""
        for s in self._sdf:
            s.info()

    #        TODO: figure how to allow [startscan, endscan]
    #            [sampler], ap_eff [if requested units are Jy]
    def getps(self, scans=None, bintable=None, **kwargs):
        """
        Get the rows that contain position-switched data.  These include ONs and OFFs.

        kwargs: plnum, feed, ifnum, integration, calibrate=T/F, average=T/F, tsys, weights

        Parameters
        ----------
        scans : int or 2-tuple
            Single scan number or list of scan numbers to use. Default: all scans.
            Scan numbers can be Ons or Offs
        weights: str
            'equal' or 'tsys' to indicate equal weighting or tsys weighting to use in time averaging. Default: 'tsys'

        Returns
        -------
        data : `~spectra.scan.ScanBlock`
            A ScanBlock containing one or more `~spectra.scan.PSScan`

        """
        # all ON/OFF scans
        kwargs_opts = {
            "ifnum": 0,
            "plnum": None,  # I prefer "pol"
            "fdnum": None,
            "calibrate": True,
            "timeaverage": False,
            "polaverage": False,
            "tsys": None,
            "weights": "tsys",
        }
        kwargs_opts.update(kwargs)

        ifnum = kwargs_opts["ifnum"]
        plnum = kwargs_opts["plnum"]
        # todo apply_selection(kwargs_opts)
        scanblock = ScanBlock()
        for i in range(len(self._sdf)):
            scanlist = self.onoff_scan_list(scans, ifnum=ifnum, plnum=plnum, fitsindex=i)
            if len(scanlist["ON"]) == 0 or len(scanlist["OFF"]) == 0:
                # print("scans not found, continuing")
                continue
            # add ifnum,plnum
            rows = self.onoff_rows(scans, ifnum=ifnum, plnum=plnum, bintable=bintable, fitsindex=i)
            # do not pass scan list here. We need all the cal rows. They will
            # be intersected with scan rows in PSScan
            # add ifnum,plnum
            calrows = self.calonoff_rows(scans=None, bintable=bintable, fitsindex=i)
            # print(f"Sending PSScan({scanlist},{rows},{calrows},{bintable}")
            g = PSScan(self._sdf[i], scanlist, rows, calrows, bintable)
            scanblock.append(g)
        if len(scanblock) == 0:
            raise Exception("Didn't find any scans matching the input selection criteria.")
        # warnings.warn("Didn't find any scans matching the input selection criteria.")
        return scanblock

    def gettp(self, scan, sig=None, cal=None, bintable=None, **kwargs):
        """
        Get a total power scan, optionally calibrating it.

        Parameters
        ----------
        scan: int
            scan number
        sig : bool or None
            True to use only integrations where signal state is True, False to use reference state (signal state is False). None to use all integrations.
        cal: bool or None
            True to use only integrations where calibration (diode) is on, False if off. None to use all integrations regardless calibration state. The system temperature will be calculated from both states regardless of the value of this variable.
        bintable : int
            the index for BINTABLE in `sdfits` containing the scans
        calibrate: bool
            whether or not to calibrate the data.  If `True`, the data will be (calon - caloff)*0.5, otherwise it will be SDFITS row data. Default:True
        weights: str
            'equal' or 'tsys' to indicate equal weighting or tsys weighting to use in time averaging. Default: 'tsys'

        scan args - ifnum, plnum, fdnum, subref

        Returns
        -------
        data : `~spectra.scan.ScanBlock`
            A ScanBlock containing one or more `~spectra.scan.TPScan`

        """
        kwargs_opts = {
            "ifnum": 0,
            "plnum": 0,
            "fdnum": None,
            "subref": None,  # subreflector position
            "timeaverage": True,
            "polaverage": True,
            "weights": "tsys",  # or 'tsys' or ndarray
            "calibrate": True,
            "debug": False,
        }
        kwargs_opts.update(kwargs)
        TF = {True: "T", False: "F"}
        sigstate = {True: "SIG", False: "REF", None: "BOTH"}
        calstate = {True: "ON", False: "OFF", None: "BOTH"}

        ifnum = kwargs_opts["ifnum"]
        plnum = kwargs_opts["plnum"]
        fdnum = kwargs_opts["fdnum"]
        subref = kwargs_opts["subref"]
        scanblock = ScanBlock()
        for i in range(len(self._sdf)):
            df = self._sdf[i]._index
            df = df[(df["SCAN"] == scan)]
            if sig is not None:
                sigch = TF[sig]
                df = df[(df["SIG"] == sigch)]
                if kwargs_opts["debug"]:
                    print("S ", len(df))
            if cal is not None:
                calch = TF[cal]
                df = df[df["CAL"] == calch]
                if kwargs_opts["debug"]:
                    print("C ", len(df))
            if ifnum is not None:
                df = df[df["IFNUM"] == ifnum]
                if kwargs_opts["debug"]:
                    print("I ", len(df))
            if plnum is not None:
                df = df[df["PLNUM"] == plnum]
                if kwargs_opts["debug"]:
                    print("P ", len(df))
            if fdnum is not None:
                df = df[df["FDNUM"] == fdnum]
                if kwargs_opts["debug"]:
                    print("F ", len(df))
            if subref is not None:
                df = df[df["SUBREF_STATE"] == subref]
                if kwargs_opts["debug"]:
                    print("SR ", len(df))
            # TBD: if ifnum is none then we will have to sort these by ifnum, plnum and store separate arrays or something.
            tprows = list(df.index)
            # data = self.rawspectra(bintable)[tprows]
            calrows = self.calonoff_rows(scans=scan, bintable=bintable, fitsindex=i, **kwargs_opts)
            if kwargs_opts["debug"]:
                print("TPROWS len=", len(tprows))
                print("CALROWS on len=", len(calrows["ON"]))
                print("fitsindex=", i)
            if len(tprows) == 0:
                continue
            g = TPScan(
                self._sdf[i], scan, sigstate[sig], calstate[cal], tprows, calrows, bintable, kwargs_opts["calibrate"]
            )
            scanblock.append(g)
        if len(scanblock) == 0:
            raise Exception("Didn't find any scans matching the input selection criteria.")
        return scanblock

    # Inspired by Dave Frayer's snodka: /users/dfrayer/gbtidlpro/snodka
    # TODO: Figure out when, if done, the fix to the Ka beam labelling took place.
    # TODO: sig and ref parameters no longer needed? Fix description of calibrate keyword, or possibly eliminate it.
    def _subbeamnod_old(self, scan, bintable=None, **kwargs):
        """Get a subbeam nod power scan, optionally calibrating it.

        Parameters
        ----------
        scan: int
            scan number
        method: str
            Method to use when processing. One of 'cycle' or 'scan'.  'cycle' is more accurate and averages data in each SUBREF_STATE cycle. 'scan' reproduces GBTIDL's snodka function which has been shown to be less accurate.  Default:'cycle'
        sig : bool
            True to indicate if this is the signal scan, False if reference
        cal: bool
            True if calibration (diode) is on, False if off.
        bintable : int
            the index for BINTABLE in `sdfits` containing the scans, None means use all bintables
        calibrate: bool
            whether or not to calibrate the data.  If `True`, the data will be (calon - caloff)*0.5, otherwise it will be SDFITS row data. Default:True
        weights: str
            'equal' or 'tsys' to indicate equal weighting or tsys weighting to use in time averaging. Default: 'tsys'

        scan args - ifnum, fdnum, subref  (plnum depends on fdnum)

        Returns
        -------
        data : `~spectra.spectrum.Spectrum`
            A Spectrum object containing the data

        """
        kwargs_opts = {
            "ifnum": 0,
            "fdnum": 0,
            "timeaverage": True,
            "weights": "tsys",  # or None or ndarray
            "calibrate": True,
            "method": "cycle",
            "debug": False,
        }
        kwargs_opts.update(kwargs)
        ifnum = kwargs_opts["ifnum"]
        fdnum = kwargs_opts["fdnum"]
        docal = kwargs_opts["calibrate"]
        w = kwargs_opts["weights"]
        method = kwargs_opts["method"]

        # Check if we are dealing with Ka data before the beam switch.
        df = self.index(bintable=bintable)
        df = df[df["SCAN"].isin([scan])]
        rx = np.unique(df["FRONTEND"])
        if len(rx) > 1:
            raise TypeError("More than one receiver for the selected scan.")
        elif rx[0] == "Rcvr26_40":  # and df["DATE-OBS"][-1] < xxxx
            # Switch the polarizations to match the beams.
            if fdnum == 0:
                plnum = 1
            elif fdnum == 1:
                plnum = 0

        scanblock = ScanBlock()
        if method == "cycle":
            # Calibrate each cycle individually and then
            # average the calibrated data.

            for sdfi in range(len(self._sdf)):
                # Row selection.
                df = self._sdf[sdfi].index(bintable=bintable)
                df = df[df["SCAN"].isin([scan])]
                df = df[df["IFNUM"].isin([ifnum])]
                df = df[df["FDNUM"].isin([fdnum])]
                df = df[df["PLNUM"].isin([plnum])]
                df_on = df[df["CAL"] == "T"]
                df_off = df[df["CAL"] == "F"]
                df_on_sig = df_on[df_on["SUBREF_STATE"] == -1]
                df_on_ref = df_on[df_on["SUBREF_STATE"] == 1]
                df_off_sig = df_off[df_off["SUBREF_STATE"] == -1]
                df_off_ref = df_off[df_off["SUBREF_STATE"] == 1]
                sig_on_rows = df_on_sig.index.to_numpy()
                ref_on_rows = df_on_ref.index.to_numpy()
                sig_off_rows = df_off_sig.index.to_numpy()
                ref_off_rows = df_off_ref.index.to_numpy()

                # Define how large of a gap between rows we will tolerate to consider
                # a row as part of a cycle.
                # Thinking about it, we should use the SUBREF_STATE=0 as delimiter rather
                # than this.
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

                # Define the calibrated data array, and
                # variables to store weights and exposure times.
                # @TODO: using TDIM7 is fragile if the headers ever change.
                #  nchan should be gotten from data length
                # e.g. len(self._hdu[1].data[:]["DATA"][row])
                # where row is a row number associated with this scan number
                df = self._sdf[sdfi].index(bintable=bintable)
                nchan = int(df["TDIM7"][0][1:-1].split(",")[0])
                ta = np.empty((len(sig_on_groups)), dtype=object)
                ta_avg = np.zeros(nchan, dtype="d")
                wt_avg = 0.0  # A single value for now, but it should be an array once we implement vector TSYS.
                tsys_wt = 0.0
                tsys_avg = 0.0
                exposure = 0.0

                # print("GROUPS ", ref_on_groups, sig_on_groups, ref_off_groups, sig_off_groups)
                # Loop over cycles, calibrating each independently.
                groups_zip = zip(ref_on_groups, sig_on_groups, ref_off_groups, sig_off_groups)

                for i, (rgon, sgon, rgoff, sgoff) in enumerate(groups_zip):
                    # Do it the dysh way.
                    calrows = {"ON": rgon, "OFF": rgoff}
                    tprows = np.sort(np.hstack((rgon, rgoff)))
                    reftp = TPScan(
                        self._sdf[sdfi], scan, "BOTH", "BOTH", tprows, calrows, bintable, kwargs_opts["calibrate"]
                    )
                    ref_avg = reftp.timeaverage(weights=w)
                    calrows = {"ON": sgon, "OFF": sgoff}
                    tprows = np.sort(np.hstack((sgon, sgoff)))
                    sigtp = TPScan(
                        self._sdf[sdfi], scan, "BOTH", "BOTH", tprows, calrows, bintable, kwargs_opts["calibrate"]
                    )
                    sig_avg = sigtp.timeaverage(weights=w)
                    # Combine sig and ref.
                    ta[i] = copy.deepcopy(sig_avg)
                    ta[i] = ta[i].new_flux_unit("K", suppress_conversion=True)
                    ta[i]._data = ((sig_avg - ref_avg) / ref_avg).flux.value * ref_avg.meta["WTTSYS"] * u.K
                    # TUNIT7 is fragile as locations of specific header data could change in the future.
                    ta[i].meta["TUNIT7"] = "Ta"
                    ta[i].meta["TSYS"] = ref_avg.meta["WTTSYS"]

                    # Add to the average, a-la accum.
                    # @TODO possibly replace with a call to util.sq_weighted_avg
                    wt_avg += ta[i].meta["TSYS"] ** -2.0
                    ta_avg[:] += ta[i].flux.value * ta[i].meta["TSYS"] ** -2.0
                    tsys_wt_ = tsys_weight(ta[i].meta["EXPOSURE"], ta[i].meta["CDELT1"], ta[i].meta["TSYS"])
                    tsys_wt += tsys_wt_
                    tsys_avg += ta[i].meta["TSYS"] * tsys_wt_
                    exposure += ta[i].meta["EXPOSURE"]

                # Divide by the sum of the weights.
                ta_avg /= wt_avg
                tsys_avg /= tsys_wt

                # Set up a Spectrum1D object to return.
                data = copy.deepcopy(ta[-1])
                data._data = ta_avg
                data.meta["TSYS"] = tsys_avg
                data.meta["EXPOSURE"] = exposure

                return data

        elif method == "scan":
            # Process the whole scan as a single block.
            # This is less accurate, but might be needed if
            # the scan was aborted and there are not enough
            # sig/ref cycles to do a per cycle calibration.

            tpon = self.gettp(
                scan,
                sig=None,
                cal=None,
                bintable=bintable,
                fdnum=fdnum,
                plnum=plnum,
                ifnum=ifnum,
                subref=-1,
                weight=w,
                calibrate=docal,
            )
            tpoff = self.gettp(
                scan,
                sig=None,
                cal=None,
                bintable=bintable,
                fdnum=fdnum,
                plnum=plnum,
                ifnum=ifnum,
                subref=1,
                weight=w,
                calibrate=docal,
            )

            on = tpon.timeaverage(weights=w)
            off = tpoff.timeaverage(weights=w)

            # in order to reproduce gbtidl tsys, we need to do a normal
            # total power scan
            fulltp = self.gettp(
                scan,
                sig=None,
                cal=None,
                bintable=bintable,
                fdnum=fdnum,
                plnum=plnum,
                ifnum=ifnum,
                weight=w,
                calibrate=docal,
            ).timeaverage(weights=w)
            tsys = fulltp.meta["TSYS"]
            data = tsys * (on - off) / off
            data.meta["MEANTSYS"] = 0.5 * np.mean((on.meta["TSYS"] + off.meta["TSYS"]))
            data.meta["WTTSYS"] = tsys
            data.meta["TSYS"] = data.meta["WTTSYS"]
            # if kwargs_opts['debug']:
            #    data._on = tpon
            #    data._off = tpoff
            return data

    def subbeamnod(self, scan, bintable=None, **kwargs):
        # TODO fix sig/cal -- no longer needed?
        """Get a subbeam nod power scan, optionally calibrating it.

        Parameters
        ----------
        scan: int
            scan number
        method: str
            Method to use when processing. One of 'cycle' or 'scan'.  'cycle' is more accurate and averages data in each SUBREF_STATE cycle. 'scan' reproduces GBTIDL's snodka function which has been shown to be less accurate.  Default:'cycle'
        sig : bool
            True to indicate if this is the signal scan, False if reference
        cal: bool
            True if calibration (diode) is on, False if off.
        bintable : int
            the index for BINTABLE in `sdfits` containing the scans, None means use all bintables
        calibrate: bool
            whether or not to calibrate the data.  If `True`, the data will be (calon - caloff)*0.5, otherwise it will be SDFITS row data. Default:True
        weights: str
            'equal' or 'tsys' to indicate equal weighting or tsys weighting to use in time averaging. Default: 'tsys'

            scan args - ifnum, fdnum, subref  (plnum depends on fdnum)

        Returns
        -------
        data : `~spectra.spectrum.Spectrum`
            A Spectrum object containing the data

        """
        kwargs_opts = {
            "ifnum": 0,
            "fdnum": 0,
            "timeaverage": True,
            "weights": "tsys",  # or None or ndarray
            "calibrate": True,
            "method": "cycle",
            "debug": False,
        }
        kwargs_opts.update(kwargs)
        ifnum = kwargs_opts["ifnum"]
        fdnum = kwargs_opts["fdnum"]
        docal = kwargs_opts["calibrate"]
        w = kwargs_opts["weights"]
        method = kwargs_opts["method"]

        # Check if we are dealing with Ka data before the beam switch.
        df = self.index(bintable=bintable)
        df = df[df["SCAN"].isin([scan])]
        rx = np.unique(df["FRONTEND"])
        if len(rx) > 1:
            raise TypeError("More than one receiver for the selected scan.")
        elif rx[0] == "Rcvr26_40":  # and df["DATE-OBS"][-1] < xxxx
            # Switch the polarizations to match the beams.
            if fdnum == 0:
                plnum = 1
            elif fdnum == 1:
                plnum = 0

        scanblock = ScanBlock()
        reftp = []
        sigtp = []
        fulltp = []
        if method == "cycle":
            # Calibrate each cycle individually and then
            # average the calibrated data.

            for sdfi in range(len(self._sdf)):
                # Row selection.
                df = self._sdf[sdfi].index(bintable=bintable)
                df = df[df["SCAN"].isin([scan])]
                df = df[df["IFNUM"].isin([ifnum])]
                df = df[df["FDNUM"].isin([fdnum])]
                df = df[df["PLNUM"].isin([plnum])]
                df_on = df[df["CAL"] == "T"]
                df_off = df[df["CAL"] == "F"]
                df_on_sig = df_on[df_on["SUBREF_STATE"] == -1]
                df_on_ref = df_on[df_on["SUBREF_STATE"] == 1]
                df_off_sig = df_off[df_off["SUBREF_STATE"] == -1]
                df_off_ref = df_off[df_off["SUBREF_STATE"] == 1]
                sig_on_rows = df_on_sig.index.to_numpy()
                ref_on_rows = df_on_ref.index.to_numpy()
                sig_off_rows = df_off_sig.index.to_numpy()
                ref_off_rows = df_off_ref.index.to_numpy()

                # Define how large of a gap between rows we will tolerate to consider
                # a row as part of a cycle.
                # Thinking about it, we should use the SUBREF_STATE=0 as delimiter rather
                # than this.
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

                # Define the calibrated data array, and
                # variables to store weights and exposure times.
                # @TODO: using TDIM7 is fragile if the headers ever change.
                #  nchan should be gotten from data length
                # e.g. len(self._hdu[1].data[:]["DATA"][row])
                # where row is a row number associated with this scan number
                df = self._sdf[sdfi].index(bintable=bintable)
                nchan = int(df["TDIM7"][0][1:-1].split(",")[0])
                ta = np.empty((len(sig_on_groups)), dtype=object)
                ta_avg = np.zeros(nchan, dtype="d")
                wt_avg = 0.0  # A single value for now, but it should be an array once we implement vector TSYS.
                tsys_wt = 0.0
                tsys_avg = 0.0
                exposure = 0.0

                # print("GROUPS ", ref_on_groups, sig_on_groups, ref_off_groups, sig_off_groups)
                # Loop over cycles, calibrating each independently.
                groups_zip = zip(ref_on_groups, sig_on_groups, ref_off_groups, sig_off_groups)

                for i, (rgon, sgon, rgoff, sgoff) in enumerate(groups_zip):
                    # Do it the dysh way.
                    calrows = {"ON": rgon, "OFF": rgoff}
                    tprows = np.sort(np.hstack((rgon, rgoff)))
                    reftp.append(
                        TPScan(
                            self._sdf[sdfi], scan, "BOTH", "BOTH", tprows, calrows, bintable, kwargs_opts["calibrate"]
                        )
                    )
                    calrows = {"ON": sgon, "OFF": sgoff}
                    tprows = np.sort(np.hstack((sgon, sgoff)))
                    sigtp.append(
                        TPScan(
                            self._sdf[sdfi], scan, "BOTH", "BOTH", tprows, calrows, bintable, kwargs_opts["calibrate"]
                        )
                    )
                sb = SubBeamNodScan(sigtp, reftp, method=method, calibrate=True, weights=w)
                scanblock.append(sb)
        elif method == "scan":
            for sdfi in range(len(self._sdf)):
                # Process the whole scan as a single block.
                # This is less accurate, but might be needed if
                # the scan was aborted and there are not enough
                # sig/ref cycles to do a per cycle calibration.

                tpon = self.gettp(
                    scan,
                    sig=None,
                    cal=None,
                    bintable=bintable,
                    fdnum=fdnum,
                    plnum=plnum,
                    ifnum=ifnum,
                    subref=-1,
                    weight=w,
                    calibrate=docal,
                )
                sigtp.append(tpon[0])
                tpoff = self.gettp(
                    scan,
                    sig=None,
                    cal=None,
                    bintable=bintable,
                    fdnum=fdnum,
                    plnum=plnum,
                    ifnum=ifnum,
                    subref=1,
                    weight=w,
                    calibrate=docal,
                )
                reftp.append(tpoff[0])
                # in order to reproduce gbtidl tsys, we need to do a normal
                # total power scan
                ftp = self.gettp(
                    scan,
                    sig=None,
                    cal=None,
                    bintable=bintable,
                    fdnum=fdnum,
                    plnum=plnum,
                    ifnum=ifnum,
                    weight=w,
                    calibrate=docal,
                )  # .timeaverage(weights=w)
                fulltp.append(ftp[0])
            sb = SubBeamNodScan(sigtp, reftp, fulltp, method=method, calibrate=True, weights=w)
            scanblock.append(sb)
        if len(scanblock) == 0:
            raise Exception("Didn't find any scans matching the input selection criteria.")
        return scanblock

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
        bintable : int
            the index for BINTABLE containing the scans
        fitsindex: int
            the index of the FITS file contained in this GBTFITSLoad.  Default:0

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
        dfon = self.select("_OBSTYPE", "PSWITCHON", df)
        dfoff = self.select("_OBSTYPE", "PSWITCHOFF", df)
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
        dfon = self.select("CAL", "T", df)
        dfoff = self.select("CAL", "F", df)
        if ifnum is not None:
            dfon = self.select("IFNUM", ifnum, dfon)
            dfoff = self.select("IFNUM", ifnum, dfoff)
        if plnum is not None:
            dfon = self.select("PLNUM", plnum, dfon)
            dfoff = self.select("PLNUM", plnum, dfoff)
        if fdnum is not None:
            dfon = self.select("FDNUM", fdnum, dfon)
            dfoff = self.select("FDNUM", fdnum, dfoff)
        if subref is not None:
            dfon = self.select("SUBREF_STATE", subref, dfon)
            dfoff = self.select("SUBREF_STATE", subref, dfoff)
        s["ON"] = list(dfon.index)
        s["OFF"] = list(dfoff.index)
        return s

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
        scans = self.onoff_scan_list(scans, ifnum, plnum, bintable, fitsindex=fitsindex)
        # scans is now a dict of "ON" "OFF
        for key in scans:
            rows[key] = self.scan_rows(scans[key], ifnum, plnum, bintable, fitsindex=fitsindex)
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
        rows = list(df.index)
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
        scanidx = self.index[self.index["SCAN"].isin(scans)]
        bt = self.udata("BINTABLE")
        for j in bt:
            df = scanidx[scanidx["BINTABLE"] == j]
            rows.append(list(df.index))
        return rows

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
        # print("Using rows",rows)
        hdu0 = self._hdu[0].copy()
        outhdu = fits.HDUList(hdu0)
        # get the bintables rows as new bintables.
        for i in range(len(rows)):
            ob = self._bintable_from_rows(rows[i], i)
            # print(f"bintable {i} #rows {len(rows[i])} data length {len(ob.data)}")
            if len(ob.data) > 0:
                outhdu.append(ob)
        # print(outhdu.info())
        # write it out!
        outhdu.update_extend()  # possibly unneeded
        outhdu.writeto(fileobj, output_verify=output_verify, overwrite=overwrite, checksum=checksum)
