"""
The classes that define various types of Scan and their calibration methods.
"""

import warnings
from collections import UserList
from copy import deepcopy

import astropy.units as u
import numpy as np
from astropy import constants as ac
from astropy.io.fits import BinTableHDU, Column
from astropy.table import Table, vstack
from scipy import ndimage

from dysh.spectra import core

from ..coordinates import Observatory
from ..util import uniq
from .core import (
    average,
    fft_shift,
    find_non_blanks,
    mean_tsys,
    sq_weighted_avg,
    tsys_weight,
)
from .spectrum import Spectrum

# from typing import Literal


# from astropy.coordinates.spectral_coordinate import NoVelocityWarning


# @todo change this to an ABC class. https://docs.python.org/3/library/abc.html
# and create a SpectrumAverageMixin for timeaverage, polaverage, finalspectrum
class ScanMixin:
    """This class describes the common interface to all Scan classes.
    A Scan represents one IF, one feed, and one or more polarizations.
    Derived classes *must* implement :meth:`calibrate`.
    """

    @property
    def scan(self):
        """
        The scan number

        Returns
        -------
        int
            The scan number of the integrations in the Scan object
        """
        return self._scan

    @property
    def nchan(self):
        """
        The number of channels in this scan

        Returns
        -------
        int
            The number of channels in this scan

        """
        return self._nchan

    @property
    def nrows(self):
        """The number of rows in this Scan

        Returns
        -------
        int
            The number of rows in this Scan
        """
        return self._nrows

    @property
    def npol(self):
        """
        The number of polarizations in this Scan

        Returns
        -------
        int
            The number of polarizations in this Scan

        """
        return self._npol

    @property
    def ifnum(self):
        """The IF number

        Returns
        -------
        int
            The index of the Intermediate Frequency
        """
        return self._ifnum

    @property
    def fdnum(self):
        """The feed number

        Returns
        -------
        int
            The index of the Feed
        """
        return self._fdnum

    @property
    def is_calibrated(self):
        """
        Have the data been calibrated?

        Returns
        -------
        bool
            True if the data have been calibrated, False if not.

        """
        return self._calibrated is not None

    @property
    def meta(self):
        """
        The metadata of this Scan. The metadata is a list of dictionaries, the length of which is
        equal to the number of calibrated integrations in the Scan.

        Returns
        -------
        dict
            Dictionary containing the metadata of this Scan

        """
        return self._meta

    def _meta_as_table(self):
        """get the metadata as an astropy Table"""
        return Table(self._meta)

    def _make_meta(self, rowindices):
        """
        Create the metadata for a Scan.  The metadata is a list of dictionaries, the length of which is
        equal to the number of calibrated integrations in the Scan.

        Parameters
        ----------
        rowindices : list of int
            The list of indices into the parent SDFITS index (DataFrame). These typically point to the
            indices for only the, e.g. the SIG data, since the resultant metadata will be used to make the calibrated spectra,
            which SIGREF data result in N/2 spectra.

        Returns
        -------
        None

        """
        # FITS can't handle NaN in a header, so just drop any column where NaN appears
        # Ideally this should be done on the individual record level in say, Spectrum.make_spectrum.
        df = self._sdfits.index(bintable=self._bintable_index).iloc[rowindices].dropna(axis=1, how="any")
        self._meta = df.to_dict("records")  # returns dict(s) with key = row number.
        for i in range(len(self._meta)):
            if "CUNIT1" not in self._meta[i]:
                self._meta[i][
                    "CUNIT1"
                ] = "Hz"  # @todo this is in gbtfits.hdu[0].header['TUNIT11'] but is it always TUNIT11?
            self._meta[i]["CUNIT2"] = "deg"  # is this always true?
            self._meta[i]["CUNIT3"] = "deg"  # is this always true?
            restfrq = self._meta[i]["RESTFREQ"]
            rfq = restfrq * u.Unit(self._meta[i]["CUNIT1"])
            restfreq = rfq.to("Hz").value
            self._meta[i]["RESTFRQ"] = restfreq  # WCS wants no E

    def _add_calibration_meta(self):
        """Add metadata that are computed after calibration."""
        if not self.is_calibrated:
            raise Exception("Data have to be calibrated first to add calibration metadata")
        for i in range(len(self._meta)):
            self._meta[i]["TSYS"] = self._tsys[i]
            self._meta[i]["EXPOSURE"] = self._exposure[i]
            self._meta[i]["NAXIS1"] = len(self._calibrated[i])
            self._meta[i]["TSYS"] = self._tsys[i]
            self._meta[i]["EXPOSURE"] = self.exposure[i]

    @property
    def pols(self):
        """The polarization number(s)

        Returns
        -------
        list
            The list of integer polarization number(s)
        """
        return self._pols

    def _set_if_fd(self, df):
        """Set the IF and FD numbers from the input dataframe and
        raise an error of there are more than one

        Parameters
        ----------
        df : ~pandas.DataFrame
            The DataFrame describing the selected data
        """
        self._ifnum = uniq(df["IFNUM"])
        self._fdnum = uniq(df["FDNUM"])
        if len(self._ifnum) > 1:
            raise Exception(f"Only one IFNUM is allowed per Scan, found {self.ifnum}")
        if len(self._fdnum) > 1:
            raise Exception(f"Only one FDNUM is allowed per Scan, found {self.fdnum}")
        self._ifnum = self._ifnum[0]
        self._fdnum = self._fdnum[0]

    def calibrate(self, **kwargs):
        """Calibrate the Scan data"""
        pass

    def timeaverage(self, weights=None):
        r"""Compute the time-averaged spectrum for this scan.

        Parameters
        ----------
        weights: str
            'tsys' or None.  If 'tsys' the weight will be calculated as:

             :math:`w = t_{exp} \times \delta\nu/T_{sys}^2`

            Default: 'tsys'
        Returns
        -------
        spectrum : :class:`~spectra.spectrum.Spectrum`
            The time-averaged spectrum
        """
        pass

    def polaverage(self, weights=None):
        """Average all polarizations in this Scan"""
        pass

    def finalspectrum(self, weights=None):
        """Average all times and polarizations in this Scan"""
        pass

    def make_bintable(self):
        """
        Create a :class:`~astropy.io.fits.BinaryTableHDU` from the calibrated data of this Scan.


        Returns
        -------
        b : :class:`~astropy.io.fits.BinaryTableHDU`
           A FITS binary table HDU, suitable for writing out or appending to a `~astropy.io.fits.HDUList`.

        """
        # Creating a Table from the metadata dictionary and instantiating
        # the BinTableHDU with that takes care of all the tricky
        # numpy/pandas dtypes to FITS single character format string. No need
        # to figure it out oneself (which I spent far too much time trying to do before I
        # discovered this!)
        if self._calibrated is None:
            raise Exception("Data must be calibrated before writing.")
        cd = BinTableHDU(data=self._meta_as_table(), name="SINGLE DISH").columns
        form = f"{np.shape(self._calibrated)[1]}E"
        cd.add_col(Column(name="DATA", format=form, array=self._calibrated))
        b = BinTableHDU.from_columns(cd, name="SINGLE DISH")
        return b

    def write(self, fileobj, output_verify="exception", overwrite=False, checksum=False):
        """
        Write an SDFITS format file (FITS binary table HDU) of the calibrated data in this Scan

        Parameters
        ----------
        fileobj : str, file-like or `pathlib.Path`
            File to write to.  If a file object, must be opened in a
            writeable mode.
        multifile: bool, optional
            If True, write to multiple files if and only if there are multiple SDFITS files in this GBTFITSLoad.
            Otherwise, write to a single SDFITS file.
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

        Returns
        -------
        None.

        """
        self.make_bintable().writeto(name=fileobj, output_verify=output_verify, overwrite=overwrite, checksum=checksum)

    def __len__(self):
        return self._nrows


class ScanBlock(UserList, ScanMixin):
    def __init__(self, *args):
        super().__init__(*args)
        self._nrows = 0
        self._npol = 0
        self._timeaveraged = []
        self._polaveraged = []
        self._finalspectrum = []

    def calibrate(self, **kwargs):
        """Calibrate all scans in this ScanBlock"""
        for scan in self.data:
            scan.calibrate(**kwargs)

    def timeaverage(self, weights="tsys", mode="old"):
        r"""Compute the time-averaged spectrum for all scans in this ScanBlock.

        Parameters
        ----------
        weights: str
            'tsys' or None.  If 'tsys' the weight will be calculated as:

             :math:`w = t_{exp} \times \delta\nu/T_{sys}^2`

            Default: 'tsys'
        Returns
        -------
        timeaverage: list of `~spectra.spectrum.Spectrum`
            List of all the time-averaged spectra
        """
        # warnings.simplefilter("ignore", NoVelocityWarning)
        if mode == "old":
            # average of the averages
            self._timeaveraged = []
            for scan in self.data:
                self._timeaveraged.append(scan.timeaverage(weights))
            if weights == "tsys":
                # There may be multiple integrations, so need to
                # average the Tsys weights
                w = np.array([np.nanmean(k._tsys_weight) for k in self.data])
                if len(np.shape(w)) > 1:  # remove empty axes
                    w = w.squeeze()
            else:
                w = weights
            timeavg = np.array([k.data for k in self._timeaveraged])
            # print(
            #    f"TAsh {np.shape(timeavg)} len(data) = {len(self.data)} weights={w}"
            # )  # " tsysW={self.data[0]._tsys_weight}")

            # Weight the average of the timeaverages by the weights.
            avgdata = average(timeavg, axis=0, weights=w)
            avgspec = np.mean(self._timeaveraged)
            avgspec.meta = self._timeaveraged[0].meta
            avgspec.meta["TSYS"] = np.average(a=[k.meta["TSYS"] for k in self._timeaveraged], axis=0, weights=w)
            avgspec.meta["EXPOSURE"] = np.sum([k.meta["EXPOSURE"] for k in self._timeaveraged])
            # observer = self._timeaveraged[0].observer # nope this has to be a location ugh. see @todo in Spectrum constructor
            # hardcode to GBT for now

            return Spectrum.make_spectrum(
                avgdata * avgspec.flux.unit, meta=avgspec.meta, observer_location=Observatory["GBT"]
            )
        elif mode == "new":
            # average of the integrations
            allcal = np.all([d._calibrate for d in self.data])
            if not allcal:
                raise Exception("Data must be calibrated before time averaging.")
            c = np.concatenate([d._calibrated for d in self.data])
            if weights == "tsys":
                w = np.concatenate([d._tsys_weight for d in self.data])
                # if len(np.shape(w)) > 1:  # remove empty axes
                #    w = w.squeeze()
            else:
                w = None
            timeavg = average(c, weights=w)
            avgspec = self.data[0].calibrated(0)
            avgspec.meta["TSYS"] = np.nanmean([d.tsys for d in self.data])
            avgspec.meta["EXPOSURE"] = np.sum([d.exposure for d in self.data])
            return Spectrum.make_spectrum(
                timeavg * avgspec.flux.unit, meta=avgspec.meta, observer_location=Observatory["GBT"]
            )
        else:
            raise Exception(f"unrecognized mode {mode}")

    def polaverage(self, weights="tsys"):
        # @todo rewrite this to return a spectrum as timeaverage does now.
        r"""Average all polarizations in all scans in this ScanBlock

        Parameters
        ----------
        weights: str
            'tsys' or None.  If 'tsys' the weight will be calculated as:

             :math:`w = t_{exp} \times \delta\nu/T_{sys}^2`

            Default: 'tsys'
        Returns
        -------
        polaverage: list of `~spectra.spectrum.Spectrum`
            List of all the polarization-averaged spectra
        """
        self._polaveraged = []
        for scan in self.data:
            self._polaveraged.append(scan.polaverage(weights))
        return self._polaveraged

    def finalspectrum(self, weights="tsys"):
        r"""Average all times and polarizations in all scans this ScanBlock

        Parameters
        ----------
        weights: str
            'tsys' or None.  If 'tsys' the weight will be calculated as:

             :math:`w = t_{exp} \times \delta\nu/T_{sys}^2`

            Default: 'tsys'
        Returns
        -------
        finalspectra: list of `~spectra.spectrum.Spectrum`
            List of all the time- and polarization-averaged spectra
        """
        self._finalspectrum = []
        for scan in self.data:
            self._finalspectrum.append(scan.finalspectrum(weights))
        return self._finalspectrum

    def write(self, fileobj, output_verify="exception", overwrite=False, checksum=False):
        """
        Write an SDFITS format file (FITS binary table HDU) of the calibrated data in this Scan

        Parameters
        ----------
        fileobj : str, file-like or `pathlib.Path`
            File to write to.  If a file object, must be opened in a
            writeable mode.
        multifile: bool, optional
            If True, write to multiple files if and only if there are multiple SDFITS files in this GBTFITSLoad.
            Otherwise, write to a single SDFITS file.
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

        Returns
        -------
        None.

        """
        s0 = self.data[0]
        # If there is only one scan, delegate to its write method
        if len(self.data) == 1:
            s0.write(fileobj, output_verify, overwrite, checksum)
            return
        # Meta are the keys of the first scan's bintable except for DATA.
        # We can use this to compare with the subsequent Scan
        # keywords without have to create their bintables first.
        defaultkeys = set(s0._meta[0].keys())
        datashape = np.shape(s0._calibrated)
        tablelist = [s0._meta_as_table()]
        for scan in self.data[1:]:
            # check data shapes are the same
            thisshape = np.shape(scan._calibrated)
            if thisshape != datashape:
                # @todo Variable length arrays? https://docs.astropy.org/en/stable/io/fits/usage/unfamiliar.html#variable-length-array-tables
                # or write to separate bintables.
                raise Exception(
                    f"Data shapes of scans are not equal {thisshape}!={datashape}. Can't combine Scans into single"
                    " BinTableHDU"
                )
            # check that the header keywords are the same
            diff = set(scan._meta[0].keys()) - defaultkeys
            if len(diff) > 0:
                raise Exception(
                    f"Scan header keywords are not the same. These keywords were not present in all Scans: {diff}."
                    " Can't combine Scans into single BinTableHDU"
                )
            tablelist.append(scan._meta_as_table())
        # now do the same trick as in Scan.write() of adding "DATA" to the coldefs
        # astropy Tables can be concatenated with vstack thankfully.
        table = vstack(tablelist, join_type="exact")
        cd = BinTableHDU(table, name="SINGLE DISH").columns
        data = np.concatenate([c._calibrated for c in self.data])
        form = f"{np.shape(data)[1]}E"
        cd.add_col(Column(name="DATA", format=form, array=data))
        b = BinTableHDU.from_columns(cd, name="SINGLE DISH")
        b.writeto(name=fileobj, output_verify=output_verify, overwrite=overwrite, checksum=checksum)


class TPScan(ScanMixin):
    """GBT specific version of Total Power Scan

    Parameters
    ----------
    gbtfits : `~fits.gbtfitsload.GBFITSLoad`
        input GBFITSLoad object
    scan: int
        scan number
    sigstate : bool
        Select the signal state used to form the data.  True means select sig='T', False to select sig='F'.
        None means select both.  See table below for explanation.
    calstate : bool
        Select the calibration state used to form the data.  True means select cal='T', False to select cal='F'.
        None means select both. See table below for explanation.
    scanrows : list-like
        the list of rows in `sdfits` corresponding to sigstate integrations
    calrows : dict
        dictionary containing with keys 'ON' and 'OFF' containing list of rows in `sdfits` corresponding to cal=T (ON) and cal=F (OFF) integrations for `scan`
    bintable : int
        the index for BINTABLE in `sdfits` containing the scans
    calibrate: bool
        whether or not to calibrate the data.  If `True`, the data will be (calon - caloff)*0.5, otherwise it will be SDFITS row data. Default:True
    smoothref: int
        the number of channels in the reference to boxcar smooth prior to calibration

    Notes
    -----
    How the total power and system temperature are calculated, depending on signal and reference state parameters:


     ======   =====   ===================================================================       ==================================
     CAL      SIG     RESULT                                                                    TSYS
     ======   =====   ===================================================================       ==================================
     None     None    data = 0.5* (REFCALON + REFCALOFF), regardless of sig state               use all CAL states, all SIG states
     None     True    data = 0.5* (REFCALON + REFCALOFF), where sig = 'T'                       use all CAL states, SIG='T'
     None     False   data = 0.5* (REFCALON + REFCALOFF), where sig = 'F'                       use all CAL states, SIG='F'
     True     None    data = REFCALON, regardless of sig state                                  use all CAL states, all SIG states
     False    None    data = REFCALOFF, regardless of sig state                                 use all CAL states, all SIG states
     True     True    data = REFCALON, where sig='T'                                            use all CAL states, SIG='T'
     True     False   data = REFCALON, where sig='F'                                            use all CAL states, SIG='F'
     False    True    data = REFCALOFF  where sig='T'                                           use all CAL states, SIG='T'
     False    False   data = REFCALOFF, where sig='F'                                           use all CAL states, SIG='F'
     ======   =====   ===================================================================       ==================================

    where `REFCALON` = integrations with `cal=T` and  `REFCALOFF` = integrations with `cal=F`.

    """

    def __init__(
        self,
        gbtfits,
        scan,
        sigstate,
        calstate,
        scanrows,
        calrows,
        bintable,
        calibrate=True,
        smoothref=1,
        observer_location=Observatory["GBT"],
    ):
        self._sdfits = gbtfits  # parent class
        self._scan = scan
        self._sigstate = sigstate
        self._calstate = calstate
        self._scanrows = scanrows
        self._smoothref = smoothref
        if self._smoothref > 1:
            warnings.warn(f"TP smoothref={self._smoothref} not implemented yet")

        # print("BINTABLE = ", bintable)
        # @todo deal with data that crosses bintables
        if bintable is None:
            self._bintable_index = self._sdfits._find_bintable_and_row(self._scanrows[0])[0]
        else:
            self._bintable_index = bintable
        self._observer_location = observer_location
        df = self._sdfits._index
        df = df.iloc[scanrows]
        self._index = df
        self._set_if_fd(df)
        self._pols = uniq(df["PLNUM"])
        self._nint = 0
        self._npol = len(self._pols)
        self._timeaveraged = None
        self._polaveraged = None
        self._nrows = len(scanrows)
        self._tsys = None
        if False:
            self._npol = gbtfits.npol(self._bintable_index)  # @todo deal with bintable
            self._nint = gbtfits.nintegrations(self._bintable_index)
        self._calrows = calrows
        # all cal=T states where sig=sigstate
        self._refonrows = sorted(list(set(self._calrows["ON"]).intersection(set(self._scanrows))))
        # all cal=F states where sig=sigstate
        self._refoffrows = sorted(list(set(self._calrows["OFF"]).intersection(set(self._scanrows))))
        self._refcalon = gbtfits.rawspectra(self._bintable_index)[self._refonrows]
        self._refcaloff = gbtfits.rawspectra(self._bintable_index)[self._refoffrows]
        # now remove blanked integrations
        # seems like this should be done for all Scan classes!
        # PS: yes.
        nb1 = find_non_blanks(self._refcalon)
        nb2 = find_non_blanks(self._refcaloff)
        goodrows = np.intersect1d(nb1, nb2)
        # Tell the user about blank integration(s) that will be ignored.
        if len(goodrows) != len(self._refcalon):
            nblanks = len(self._refcalon) - len(goodrows)
            print(f"Ignoring {nblanks} blanked integration(s).")
        self._refcalon = self._refcalon[goodrows]
        self._refcaloff = self._refcaloff[goodrows]
        self._refonrows = [self._refonrows[i] for i in goodrows]
        self._refoffrows = [self._refoffrows[i] for i in goodrows]
        self._nchan = len(self._refcalon[0])
        self._calibrate = calibrate
        self._data = None
        if self._calibrate:
            self.calibrate()
        self.calc_tsys()

    def calibrate(self):
        """Calibrate the data according to the CAL/SIG table above"""
        # the way the data are formed depend only on cal state
        # since we have downselected based on sig state in the constructor
        if self.calstate is None:
            # print("data = 0.5*(ON+OFF)")
            self._data = (0.5 * (self._refcalon + self._refcaloff)).astype(float)
        elif self.calstate:
            # print("data = ON")
            self._data = self._refcalon.astype(float)
        elif self.calstate == False:
            # print("data = OFF")
            self._data = self._refcaloff.astype(float)
        else:
            raise Exception(f"Unrecognized cal state {self.calstate}")  # should never happen

    @property
    def sigstate(self):
        """The requested signal state

        Returns
        -------
        bool
            True if signal state is on ('T' in the SDFITS header), False otherwise ('F')
        """
        return self._sigstate

    @property
    def calstate(self):
        """The requested calibration state

        Returns
        -------
        bool
            True if calibration state is on ('T' in the SDFITS header), False otherwise ('F')
        """
        return self._calstate

    @property
    def tsys(self):
        """The system temperature array.

        Returns
        -------
        tsys : `~numpy.ndarray`
            System temperature values in K
        """
        return self._tsys

    def calc_tsys(self, **kwargs):
        """
        Calculate the system temperature array, according to table above.
        """
        kwargs_opts = {"verbose": False}
        kwargs_opts.update(kwargs)

        if False:
            if self.calstate is None:
                tcal = list(self._sdfits.index(bintable=self._bintable_index).iloc[self._refonrows]["TCAL"])
                nspect = len(tcal)
                calon = self._refcalcon
                caloff = self._refcaloff
            elif self.calstate:
                tcal = list(self._sdfits.index(bintable=self._bintable_index).iloc[self._refonrows]["TCAL"])
                nspect = len(tcal)
                calon = self._refcalon
            elif self.calstate == False:
                pass
        self._tcal = list(self._sdfits.index(bintable=self._bintable_index).iloc[self._refonrows]["TCAL"])
        nspect = len(self._tcal)
        self._tsys = np.empty(nspect, dtype=float)  # should be same as len(calon)
        # allcal = self._refonrows.copy()
        # allcal.extend(self._refoffrows)
        # tcal = list(self._sdfits.index(self._bintable_index).iloc[sorted(allcal)]["TCAL"])
        # @todo this loop could be replaced with clever numpy
        if len(self._tcal) != nspect:
            raise Exception(f"TCAL length {len(tcal)} and number of spectra {nspect} don't match")
        for i in range(nspect):
            # tsys = mean_tsys(calon=calon[i], caloff=caloff[i], tcal=tcal[i])
            tsys = mean_tsys(calon=self._refcalon[i], caloff=self._refcaloff[i], tcal=self._tcal[i])
            self._tsys[i] = tsys

    @property
    def exposure(self):
        """Get the array of exposure (integration) times.  The value depends on the cal state:

            =====  ======================================
            CAL    EXPOSURE
            =====  ======================================
            None   :math:`t_{EXP,REFON} + t_{EXP,REFOFF}`
            True   :math:`t_{EXP,REFON}`
            False  :math:`t_{EXP,REFOFF}`
            =====  ======================================

        Returns
        -------
        exposure : `~numpy.ndarray`
            The exposure time in units of the EXPOSURE keyword in the SDFITS header
        """
        if self.calstate is None:
            exp_ref_on = self._sdfits.index(bintable=self._bintable_index).iloc[self._refonrows]["EXPOSURE"].to_numpy()
            exp_ref_off = (
                self._sdfits.index(bintable=self._bintable_index).iloc[self._refoffrows]["EXPOSURE"].to_numpy()
            )
        elif self.calstate:
            exp_ref_on = self._sdfits.index(bintable=self._bintable_index).iloc[self._refonrows]["EXPOSURE"].to_numpy()
            exp_ref_off = 0
        elif self.calstate == False:
            exp_ref_on = 0
            exp_ref_off = (
                self._sdfits.index(bintable=self._bintable_index).iloc[self._refoffrows]["EXPOSURE"].to_numpy()
            )

        exposure = exp_ref_on + exp_ref_off
        return exposure

    @property
    def delta_freq(self):
        r"""Get the array of channel frequency width. The value depends on the cal state:


        =====  ================================================================
        CAL     :math:`\Delta\nu`
        =====  ================================================================
        None    :math:`0.5 * ( \Delta\nu_{REFON}+ \Delta\nu_{REFOFF} )`
        True    :math:`\Delta\nu_{REFON}`
        False   :math:`\Delta\nu_{REFOFF}`
        =====  ================================================================


        Returns
        -------
        delta_freq: `~numpy.ndarray`
            The channel frequency width in units of the CDELT1 keyword in the SDFITS header
        """
        df_ref_on = self._sdfits.index(bintable=self._bintable_index).iloc[self._refonrows]["CDELT1"].to_numpy()
        df_ref_off = self._sdfits.index(bintable=self._bintable_index).iloc[self._refoffrows]["CDELT1"].to_numpy()
        if self.calstate is None:
            delta_freq = 0.5 * (df_ref_on + df_ref_off)
        elif self.calstate:
            delta_freq = df_ref_on
        elif self.calstate == False:
            delta_freq = df_ref_off
        return delta_freq

    @property
    def _tsys_weight(self):
        r"""The system temperature weighting array computed from current
        :math:`T_{sys}, t_{exp}`, and `\delta\nu`. See :meth:`tsys_weight`
        """
        return tsys_weight(self.exposure, self.delta_freq, self.tsys)

    def tpmeta(self, i):
        ser = self._sdfits.index(bintable=self._bintable_index).iloc[self._scanrows[i]]
        # meta = self._sdfits.index(bintable=self._bintable_index).iloc[self._scanrows[i]].dropna().to_dict()
        meta = ser.dropna().to_dict()
        meta["TSYS"] = self._tsys[i]
        meta["EXPOSURE"] = self.exposure[i]
        meta["NAXIS1"] = len(self._data[i])
        if "CUNIT1" not in meta:
            meta["CUNIT1"] = "Hz"  # @todo this is in gbtfits.hdu[0].header['TUNIT11'] but is it always TUNIT11?
        meta["CUNIT2"] = "deg"  # is this always true?
        meta["CUNIT3"] = "deg"  # is this always true?
        restfrq = meta["RESTFREQ"]
        rfq = restfrq * u.Unit(meta["CUNIT1"])
        restfreq = rfq.to("Hz").value
        meta["RESTFRQ"] = restfreq  # WCS wants no E
        return meta

    def total_power(self, i):
        """Return the total power spectrum

        Parameters
        ----------
        i : int
            The index into the data array

        Returns
        -------
        spectrum : `~spectra.spectrum.Spectrum`
        """
        if not self._calibrate:
            raise Exception("You must calibrate first to get a total power spectrum")
        # print(len(self._scanrows), i)
        ser = self._sdfits.index(bintable=self._bintable_index).iloc[self._scanrows[i]]
        # meta = self._sdfits.index(bintable=self._bintable_index).iloc[self._scanrows[i]].dropna().to_dict()
        meta = ser.dropna().to_dict()
        meta["TSYS"] = self._tsys[i]
        meta["EXPOSURE"] = self.exposure[i]
        meta["NAXIS1"] = len(self._data[i])
        if "CUNIT1" not in meta:
            meta["CUNIT1"] = "Hz"  # @todo this is in gbtfits.hdu[0].header['TUNIT11'] but is it always TUNIT11?
        meta["CUNIT2"] = "deg"  # is this always true?
        meta["CUNIT3"] = "deg"  # is this always true?
        restfrq = meta["RESTFREQ"]
        rfq = restfrq * u.Unit(meta["CUNIT1"])
        restfreq = rfq.to("Hz").value
        meta["RESTFRQ"] = restfreq  # WCS wants no E
        return Spectrum.make_spectrum(self._data[i] * u.ct, meta, observer_location=self._observer_location)

    def timeaverage(self, weights="tsys"):
        r"""Compute the time-averaged spectrum for this set of scans.

        Parameters
        ----------
        weights: str
            'tsys' or None.  If 'tsys' the weight will be calculated as:

             :math:`w = t_{exp} \times \delta\nu/T_{sys}^2`

            Default: 'tsys'
        Returns
        -------
        spectrum : :class:`~spectra.spectrum.Spectrum`
            The time-averaged spectrum
        """
        if self._npol > 1:
            raise Exception("Can't yet time average multiple polarizations")
        self._timeaveraged = deepcopy(self.total_power(0))
        if weights == "tsys":
            w = self._tsys_weight
        else:
            w = np.ones_like(self._tsys_weight)
        non_blanks = find_non_blanks(self._data)[0]
        self._timeaveraged._data = average(self._data, axis=0, weights=w)
        self._timeaveraged.meta["MEANTSYS"] = np.mean(self._tsys[non_blanks])
        self._timeaveraged.meta["WTTSYS"] = sq_weighted_avg(self._tsys[non_blanks], axis=0, weights=w[non_blanks])
        self._timeaveraged.meta["TSYS"] = self._timeaveraged.meta["WTTSYS"]
        self._timeaveraged.meta["EXPOSURE"] = self.exposure[non_blanks].sum()
        return self._timeaveraged


#        @todo   'scans' should become 'scan'
class PSScan(ScanMixin):
    """GBT specific version of Position Switch Scan. A position switch scan object has
    one IF, one feed, and one or more polarizations.

    Parameters
    ----------
    gbtfits : `~fits.gbtfitsload.GBFITSLoad`
        input GBFITSLoad object
    scans : dict
        dictionary with keys 'ON' and 'OFF' containing unique list of ON (signal) and OFF (reference) scan numbers NOTE: there should be one ON and one OFF, a pair
    scanrows : dict
        dictionary with keys 'ON' and 'OFF' containing the list of rows in `sdfits` corresponding to ON (signal) and OFF (reference) integrations
    calrows : dict
        dictionary containing with keys 'ON' and 'OFF' containing list of rows in `sdfits` corresponding to cal=T (ON) and cal=F (OFF) integrations.
    bintable : int
        the index for BINTABLE in `sdfits` containing the scans
    calibrate: bool
        whether or not to calibrate the data.  If true, data will be calibrated as TSYS*(ON-OFF)/OFF. Default: True
    smoothref: int
        the number of channels in the reference to boxcar smooth prior to calibration
    observer_location : `~astropy.coordinates.EarthLocation`
        Location of the observatory. See `~dysh.coordinates.Observatory`.
        This will be transformed to `~astropy.coordinates.ITRS` using the time of
        observation DATE-OBS or MJD-OBS in
        the SDFITS header.  The default is the location of the GBT.
    """

    def __init__(
        self,
        gbtfits,
        scans,
        scanrows,
        calrows,
        bintable,
        calibrate=True,
        smoothref=1,
        observer_location=Observatory["GBT"],
    ):
        # The rows of the original bintable corresponding to ON (sig) and OFF (reg)
        self._sdfits = gbtfits  # parent class
        self._scans = scans
        self._scan = scans["ON"]
        self._scanrows = scanrows
        self._nrows = len(self._scanrows["ON"])
        self._smoothref = smoothref
        # print(f"PJT len(scanrows ON) {len(self._scanrows['ON'])}")
        # print(f"PJT len(scanrows OFF) {len(self._scanrows['OFF'])}")
        # print("PJT scans", scans)
        # print("PJT scanrows", scanrows)
        # print("PJT calrows", calrows)
        # print(f"len(scanrows ON) {len(self._scanrows['ON'])}")
        # print(f"len(scanrows OFF) {len(self._scanrows['OFF'])}")

        # calrows perhaps not needed as input since we can get it from gbtfits object?
        # calrows['ON'] are rows with noise diode was on, regardless of sig or ref
        # calrows['OFF'] are rows with noise diode was off, regardless of sig or ref
        self._calrows = calrows
        # print("BINTABLE = ", bintable)
        # @todo deal with data that crosses bintables
        if bintable is None:
            self._bintable_index = gbtfits._find_bintable_and_row(self._scanrows["ON"][0])[0]
        else:
            self._bintable_index = bintable
        # print(f"bintable index is {self._bintable_index}")
        self._observer_location = observer_location
        # df = selection.iloc[scanrows["ON"]]
        df = self._sdfits._index.iloc[scanrows["ON"]]
        self._set_if_fd(df)
        self._pols = uniq(df["PLNUM"])
        self._npol = len(self._pols)
        if False:
            self._nint = gbtfits.nintegrations(self._bintable_index)
        # so quick with slicing!
        self._sigonrows = sorted(list(set(self._calrows["ON"]).intersection(set(self._scanrows["ON"]))))
        self._sigoffrows = sorted(list(set(self._calrows["OFF"]).intersection(set(self._scanrows["ON"]))))
        self._refonrows = sorted(list(set(self._calrows["ON"]).intersection(set(self._scanrows["OFF"]))))
        self._refoffrows = sorted(list(set(self._calrows["OFF"]).intersection(set(self._scanrows["OFF"]))))
        self._sigcalon = gbtfits.rawspectra(self._bintable_index)[self._sigonrows]
        self._nchan = len(self._sigcalon[0])
        self._sigcaloff = gbtfits.rawspectra(self._bintable_index)[self._sigoffrows]
        self._refcalon = gbtfits.rawspectra(self._bintable_index)[self._refonrows]
        self._refcaloff = gbtfits.rawspectra(self._bintable_index)[self._refoffrows]
        self._tsys = None
        self._exposure = None
        self._calibrated = None
        self._calibrate = calibrate
        self._make_meta(self._sigonrows)
        if self._calibrate:
            self.calibrate()

    @property
    def scans(self):
        """The dictionary of the ON and OFF scan numbers in the PSScan.

        Returns
        -------
        scans : dict
            The scan number dictionary
        """
        return self._scans

    @property
    def tsys(self):
        """The system temperature array. This will be `None` until calibration is done.

        Returns
        -------
        tsys : `~numpy.ndarray`
            System temperature values in K
        """
        return self._tsys

    # @todo something clever
    # self._calibrated_spectrum = Spectrum(self._calibrated,...) [assuming same spectral axis]
    def calibrated(self, i):
        """Return the calibrated Spectrum.

        Parameters
        ----------
        i : int
            The index into the calibrated array

        Returns
        -------
        spectrum : `~spectra.spectrum.Spectrum`
        """
        return Spectrum.make_spectrum(
            self._calibrated[i] * u.K, meta=self.meta[i], observer_location=self._observer_location
        )

    def calibrate(self, **kwargs):
        """
        Position switch calibration, following equations 1 and 2 in the GBTIDL calibration manual
        """
        kwargs_opts = {"verbose": False}
        kwargs_opts.update(kwargs)
        if self._smoothref > 1 and kwargs_opts["verbose"]:
            print(f"PS smoothref={self._smoothref}")

        self._status = 1
        nspect = self.nrows // 2
        self._calibrated = np.empty((nspect, self._nchan), dtype="d")
        self._tsys = np.empty(nspect, dtype="d")
        self._exposure = np.empty(nspect, dtype="d")
        # print("REFONROWS ", self._refonrows)
        tcal = list(self._sdfits.index(bintable=self._bintable_index).iloc[self._refonrows]["TCAL"])
        # @todo  this loop could be replaced with clever numpy
        if len(tcal) != nspect:
            raise Exception(f"TCAL length {len(tcal)} and number of spectra {nspect} don't match")
        for i in range(nspect):
            tsys = mean_tsys(calon=self._refcalon[i], caloff=self._refcaloff[i], tcal=tcal[i])
            sig = 0.5 * (self._sigcalon[i] + self._sigcaloff[i])
            ref = 0.5 * (self._refcalon[i] + self._refcaloff[i])
            if self._smoothref > 1:
                ref = core.smooth(ref, "boxcar", self._smoothref)
            self._calibrated[i] = tsys * (sig - ref) / ref
            self._tsys[i] = tsys
            self._exposure[i] = self.exposure[i]
        # print("Calibrated %d spectra" % nspect)
        self._add_calibration_meta()

    # tip o' the hat to Pedro S. for exposure and delta_freq
    @property
    def exposure(self):
        """Get the array of exposure (integration) times

        exposure = [ 0.5*(exp_ref_on + exp_ref_off) + 0.5*(exp_sig_on + exp_sig_off) ] / 2

        Returns
        -------
        exposure : ~numpy.ndarray
            The exposure time in units of the EXPOSURE keyword in the SDFITS header
        """
        exp_ref_on = self._sdfits.index(bintable=self._bintable_index).iloc[self._refonrows]["EXPOSURE"].to_numpy()
        exp_ref_off = self._sdfits.index(bintable=self._bintable_index).iloc[self._refoffrows]["EXPOSURE"].to_numpy()
        exp_sig_on = self._sdfits.index(bintable=self._bintable_index).iloc[self._sigonrows]["EXPOSURE"].to_numpy()
        exp_sig_off = self._sdfits.index(bintable=self._bintable_index).iloc[self._sigoffrows]["EXPOSURE"].to_numpy()
        exp_ref = exp_ref_on + exp_ref_off
        exp_sig = exp_sig_on + exp_sig_off
        if self._smoothref > 1:
            nsmooth = self._smoothref
        else:
            nsmooth = 1.0
        exposure = exp_sig * exp_ref * nsmooth / (exp_sig + exp_ref * nsmooth)
        return exposure

    @property
    def delta_freq(self):
        """Get the array of channel frequency width

        df = [ 0.5*(df_ref_on + df_ref_off) + 0.5*(df_sig_on + df_sig_off) ] / 2

        Returns
        -------
             delta_freq: ~numpy.ndarray
                 The channel frequency width in units of the CDELT1 keyword in the SDFITS header
        """
        df_ref_on = self._sdfits.index(bintable=self._bintable_index).iloc[self._refonrows]["CDELT1"].to_numpy()
        df_ref_off = self._sdfits.index(bintable=self._bintable_index).iloc[self._refoffrows]["CDELT1"].to_numpy()
        df_sig_on = self._sdfits.index(bintable=self._bintable_index).iloc[self._sigonrows]["CDELT1"].to_numpy()
        df_sig_off = self._sdfits.index(bintable=self._bintable_index).iloc[self._sigoffrows]["CDELT1"].to_numpy()
        df_ref = 0.5 * (df_ref_on + df_ref_off)
        df_sig = 0.5 * (df_sig_on + df_sig_off)
        delta_freq = 0.5 * (df_ref + df_sig)
        return delta_freq

    @property
    def _tsys_weight(self):
        r"""The system temperature weighting array computed from current
        :math`T_{sys}`, :math:`t_{int}`, and :math:`\delta\nu`. See :meth:`tsys_weight`
        """
        return tsys_weight(self.exposure, self.delta_freq, self.tsys)

    def timeaverage(self, weights="tsys"):
        r"""Compute the time-averaged spectrum for this set of scans.

        Parameters
        ----------
        weights: str
            'tsys' or None.  If 'tsys' the weight will be calculated as:

             :math:`w = t_{exp} \times \delta\nu/T_{sys}^2`

            Default: 'tsys'
        Returns
        -------
        spectrum : :class:`~spectra.spectrum.Spectrum`
            The time-averaged spectrum
        """
        if self._calibrated is None or len(self._calibrated) == 0:
            raise Exception("You can't time average before calibration.")
        if self._npol > 1:
            raise Exception("Can't yet time average multiple polarizations")
        self._timeaveraged = deepcopy(self.calibrated(0))
        data = self._calibrated
        if weights == "tsys":
            w = self._tsys_weight
        else:
            w = np.ones_like(self._tsys_weight)
        self._timeaveraged._data = average(data, axis=0, weights=w)
        non_blanks = find_non_blanks(data)
        self._timeaveraged.meta["MEANTSYS"] = np.mean(self._tsys[non_blanks])
        self._timeaveraged.meta["WTTSYS"] = sq_weighted_avg(self._tsys[non_blanks], axis=0, weights=w[non_blanks])
        self._timeaveraged.meta["EXPOSURE"] = np.sum(self._exposure[non_blanks])
        self._timeaveraged.meta["TSYS"] = self._timeaveraged.meta["WTTSYS"]
        return self._timeaveraged


class FSScan(ScanMixin):
    """GBT specific version of Frequency Switch Scan

    Parameters
    ----------

    gbtfits : `~fit.gbtfitsload.GBFITSLoad`
        Input GBFITSLoad object.
    scan : int
        Scan number that contains integrations with a series of sig/ref and calon/caloff states.
    sigrows :dict
        Dictionary containing with keys 'ON' and 'OFF' containing list of rows in `sdfits`
        corresponding to sig=T (ON) and sig=F (OFF) integrations.
    calrows : dict
        Dictionary containing with keys 'ON' and 'OFF' containing list of rows in `sdfits`
        corresponding to cal=T (ON) and cal=F (OFF) integrations.
    bintable : int
        The index for BINTABLE in `sdfits` containing the scans.
    calibrate : bool
        Whether or not to calibrate the data.  If true, data will be calibrated as TSYS*(ON-OFF)/OFF.
        Default: True
    fold : bool
        Whether or not to fold the spectrum. Default: True
    shift_method : str
        Method to use when shifting the spectra for folding. One of 'fft' or 'interpolate'.
        'fft' uses a phase shift in the time domain. 'interpolate' interpolates the signal. Default: 'fft'
    use_sig : bool
        Whether to use the sig as the sig, or the ref as the sig. Default: True
    smoothref: int
        The number of channels in the reference to boxcar smooth prior to calibration.
    observer_location : `~astropy.coordinates.EarthLocation`
        Location of the observatory. See `~dysh.coordinates.Observatory`.
        This will be transformed to `~astropy.coordinates.ITRS` using the time of
        observation DATE-OBS or MJD-OBS in the SDFITS header.  The default is the location of the GBT.
    """

    def __init__(
        self,
        gbtfits,
        scan,
        sigrows,
        calrows,
        bintable,
        calibrate=True,
        fold=True,
        shift_method="fft",
        use_sig=True,
        smoothref=1,
        observer_location=Observatory["GBT"],
        debug=False,
    ):
        # The rows of the original bintable corresponding to ON (sig) and OFF (reg)
        self._sdfits = gbtfits  # parent class
        self._scan = scan  # for FS everything is an "ON"
        self._sigrows = sigrows  # dict with "ON" and "OFF"
        self._calrows = calrows  # dict with "ON" and "OFF"
        self._folded = False
        self._use_sig = use_sig
        self._smoothref = smoothref
        if self._smoothref > 1:
            print(f"FS smoothref={self._smoothref} not implemented yet")

        self._sigonrows = sorted(list(set(self._calrows["ON"]).intersection(set(self._sigrows["ON"]))))
        self._sigoffrows = sorted(list(set(self._calrows["OFF"]).intersection(set(self._sigrows["ON"]))))
        self._refonrows = sorted(list(set(self._calrows["ON"]).intersection(set(self._sigrows["OFF"]))))
        self._refoffrows = sorted(list(set(self._calrows["OFF"]).intersection(set(self._sigrows["OFF"]))))

        self._debug = debug

        if self._debug:
            print("---------------------------------------------------")
            print("FSSCAN: ")
            print("SigOff", self._sigoffrows)
            print("SigOn", self._sigonrows)
            print("RefOff", self._refoffrows)
            print("RegOn", self._refonrows)

        nsigrows = len(self._sigonrows) + len(self._sigoffrows)
        nrefrows = len(self._refonrows) + len(self._refoffrows)
        if nsigrows != nrefrows:
            raise Exception("Number of sig rows does not match ref rows. Dangerous to proceed")
        if self._debug:
            print("sigonrows", nsigrows, self._sigonrows)
        self._nrows = nsigrows

        a_scanrow = self._sigonrows[0]

        # print("BINTABLE = ", bintable)
        # @todo deal with data that crosses bintables
        if bintable is None:
            self._bintable_index = gbtfits._find_bintable_and_row(a_scanrow)[0]
        else:
            self._bintable_index = bintable
        if self._debug:
            print(f"bintable index is {self._bintable_index}")
        self._observer_location = observer_location
        self._scanrows = list(set(self._calrows["ON"])) + list(set(self._calrows["OFF"]))

        df = self._sdfits._index.iloc[self._scanrows]
        if self._debug:
            print("len(df) = ", len(df))
        self._set_if_fd(df)
        self._pols = uniq(df["PLNUM"])
        if self._debug:
            print(f"FSSCAN #pol = {self._pols}")
        self._npol = len(self._pols)
        if False:
            self._nint = gbtfits.nintegrations(self._bintable_index)
        # @todo use gbtfits.velocity_convention(veldef,velframe)
        # so quick with slicing!

        self._sigcalon = gbtfits.rawspectra(self._bintable_index)[self._sigonrows]
        self._sigcaloff = gbtfits.rawspectra(self._bintable_index)[self._sigoffrows]
        self._refcalon = gbtfits.rawspectra(self._bintable_index)[self._refonrows]
        self._refcaloff = gbtfits.rawspectra(self._bintable_index)[self._refoffrows]
        self._nchan = len(self._sigcalon[0])
        self._tsys = None
        self._exposure = None
        self._calibrated = None
        self._calibrate = calibrate
        self._make_meta(self._sigonrows)
        if self._calibrate:
            self.calibrate(fold=fold, shift_method=shift_method)
        if self._debug:
            print("---------------------------------------------------")

    @property
    def folded(self):
        """
        Has the FSscan been folded?

        Returns
        -------
        boolean
            True if the signal and reference integrations have been folded. False if not.
        """
        return self._folded

    @property
    def tsys(self):
        """The system temperature array. This will be `None` until calibration is done.

        Returns
        -------
        tsys : `~numpy.ndarray`
            System temperature values in K
        """
        return self._tsys

    def calibrated(self, i):
        """Return the calibrated Spectrum of this FSscan

        Parameters
        ----------
        i : int
            The index into the calibrated array

        Returns
        -------
        spectrum : `~spectra.spectrum.Spectrum`
        """
        return Spectrum.make_spectrum(
            self._calibrated[i] * u.K, meta=self.meta[i], observer_location=self._observer_location
        )

    def calibrate(self, **kwargs):
        """
        Frequency switch calibration, following equations ...
        fold=True or fold=False is required
        """
        if self._debug:
            print(f'FOLD={kwargs["fold"]}')
            print(f'METHOD={kwargs["shift_method"]}')

        # some helper functions, courtesy proto_getfs.py
        def channel_to_frequency(crval1, crpix1, cdelt1, vframe, nchan, nint, ndim=1):
            """ """

            # Compute the correction factor.
            beta = (vframe * u.m / u.s) / ac.c
            vcorr = np.sqrt((1.0 + beta) / (1.0 - beta))

            # The +1 is to start counting from 1.
            indx = np.arange(nchan) + 1
            if ndim == 1:
                freq = crval1 + cdelt1 * (indx - crpix1)
                freq *= vcorr
            elif ndim == 2:
                indx = np.tile(indx, (nint, 1))
                freq = crval1[:, np.newaxis] + cdelt1[:, np.newaxis] * (indx - crpix1[:, np.newaxis])
                freq *= vcorr[:, np.newaxis]

            return freq

        def index_frequency(df):
            """
            Create a frequency axis from an index.
            This assumes all entries in the index have the same number of channels.
            """
            # Could you do this with gbtfits.getspec(row).spectral_axis?

            ndim = len(df.shape)
            nint = df.shape[0]

            if ndim == 1:
                nchan = np.array([int(df["TDIM7"][1:-1].split(",")[0])])
            else:
                nchan = np.array([int(df["TDIM7"].iloc[i][1:-1].split(",")[0]) for i in range(len(df))])

            crval1 = df["CRVAL1"]
            crpix1 = df["CRPIX1"]
            cdelt1 = df["CDELT1"]
            vframe = df["VFRAME"]  # Use the velocity frame requested by the user.

            if ndim == 2:
                crval1 = crval1.to_numpy()
                crpix1 = crpix1.to_numpy()
                cdelt1 = cdelt1.to_numpy()
                vframe = vframe.to_numpy()

            freq = channel_to_frequency(crval1, crpix1, cdelt1, vframe, nchan[0], nint, ndim=ndim)

            # Apply units.
            try:
                cunit1 = u.Unit(df["CUNIT1"])
                if ndim == 2:
                    cunit1 = cunit[0]
            except KeyError:
                cunit1 = u.Hz

            return freq * cunit1

        def do_total_power(no_cal, cal, tcal):
            """ """

        def vec_mean_tsys(on, off, tcal):
            """
            mean_tsys implements this, albeit only in 1D
            """
            pass

        def do_sig_ref(sig, ref, tsys, smooth=False):
            """
            smooth=True would implement smoothing the reference (or something)
            """
            return (sig - ref) / ref * tsys

        def do_fold(sig, ref, sig_freq, ref_freq, remove_wrap=False, shift_method="fft"):
            """ """
            chan_shift = (sig_freq[0] - ref_freq[0]) / np.abs(np.diff(sig_freq)).mean()
            # print("do_fold: ",sig_freq[0], ref_freq[0],chan_shift)
            ref_shift = do_shift(ref, chan_shift, remove_wrap=remove_wrap, method=shift_method)
            # @todo weights
            avg = (sig + ref_shift) / 2
            return avg

        def do_shift(data, offset, remove_wrap=False, method="fft"):
            """
            Shift the data of a numpy array using roll/shift

            @todo   use the fancier GBTIDL fft based shift
            """

            ishift = int(np.round(offset))  # Integer shift.
            fshift = offset - ishift  # Fractional shift.

            # print("FOLD:  ishift=%d fshift=%g" % (ishift, fshift))
            data2 = np.roll(data, ishift, axis=0)
            if remove_wrap:
                if ishift < 0:
                    data2[ishift:] = np.nan
                else:
                    data2[:ishift] = np.nan
            # Now the fractional shift, each row separate since ndimage.shift() cannot deal with np.nan
            if method == "fft":
                # Set `pad=False` to avoid edge effects.
                # This needs to be sorted out.
                data2 = fft_shift(data2, fshift, pad=False)
            elif method == "interpolate":
                data2 = ndimage.shift(data2, [fshift])
            return data2

        kwargs_opts = {"verbose": False}
        kwargs_opts.update(kwargs)
        _fold = kwargs.get("fold", False)
        _mode = 1  # 1: keep the sig    else: keep the ref     (not externally supported)
        nspect = self.nrows // 2
        self._calibrated = np.empty((nspect, self._nchan), dtype="d")
        self._tsys = np.empty(nspect, dtype="d")
        self._exposure = np.empty(nspect, dtype="d")
        #
        sig_freq = self._sigcalon[0]
        df_sig = self._sdfits.index(bintable=self._bintable_index).iloc[self._sigonrows]
        df_ref = self._sdfits.index(bintable=self._bintable_index).iloc[self._refonrows]
        if self._debug:
            print("df_sig", type(df_sig), len(df_sig))
        sig_freq = index_frequency(df_sig)
        ref_freq = index_frequency(df_ref)
        chan_shift = abs(sig_freq[0, 0] - ref_freq[0, 0]) / np.abs(np.diff(sig_freq)).mean()
        if self._debug:
            print("FS: shift=%g  nchan=%d" % (chan_shift, self._nchan))

        #  tcal is the same for REF and SIG, and the same for all integrations actually.
        tcal = list(self._sdfits.index(bintable=self._bintable_index).iloc[self._sigonrows]["TCAL"])
        if self._debug:
            print("TCAL:", len(tcal), tcal[0])
        if len(tcal) != nspect:
            raise Exception(f"TCAL length {len(tcal)} and number of spectra {nspect} don't match")
        # @todo   the nspect loop could be replaced with clever numpy?
        for i in range(nspect):
            tsys_sig = mean_tsys(calon=self._sigcalon[i], caloff=self._sigcaloff[i], tcal=tcal[i])
            tsys_ref = mean_tsys(calon=self._refcalon[i], caloff=self._refcaloff[i], tcal=tcal[i])
            if i == 0 and self._debug:
                print("Tsys(sig/ref)[0]=", tsys_sig, tsys_ref)
            tp_sig = 0.5 * (self._sigcalon[i] + self._sigcaloff[i])
            tp_ref = 0.5 * (self._refcalon[i] + self._refcaloff[i])
            #
            cal_sig = do_sig_ref(tp_sig, tp_ref, tsys_ref)
            cal_ref = do_sig_ref(tp_ref, tp_sig, tsys_sig)
            #
            if _fold:
                cal_sig_fold = do_fold(cal_sig, cal_ref, sig_freq[i], ref_freq[i], shift_method=kwargs["shift_method"])
                cal_ref_fold = do_fold(cal_ref, cal_sig, ref_freq[i], sig_freq[i], shift_method=kwargs["shift_method"])
                self._folded = True
                if self._use_sig:
                    self._calibrated[i] = cal_sig_fold
                    self._tsys[i] = tsys_ref
                else:
                    self._calibrated[i] = cal_ref_fold
                    self._tsys[i] = tsys_sig
                self._exposure[i] = 2 * self.exposure[i]  # @todo
            else:
                if self._use_sig:
                    self._calibrated[i] = cal_sig
                    self._tsys[i] = tsys_ref
                else:
                    self._calibrated[i] = cal_ref
                    self._tsys[i] = tsys_sig
                self._exposure[i] = self.exposure[i]

        self._add_calibration_meta()
        # print("Calibrated %d spectra with fold=%s and use_sig=%s" % (nspect, repr(_fold), repr(self._use_sig)))

    # tip o' the hat to Pedro S. for exposure and delta_freq
    @property
    def exposure(self):
        """Get the array of exposure (integration) times for FSscan

        exposure = [ 0.5*(exp_ref_on + exp_ref_off) + 0.5*(exp_sig_on + exp_sig_off) ] / 2

        Returns
        -------
        exposure : ~numpy.ndarray
            The exposure time in units of the EXPOSURE keyword in the SDFITS header
        """
        exp_ref_on = self._sdfits.index(bintable=self._bintable_index).iloc[self._refonrows]["EXPOSURE"].to_numpy()
        exp_ref_off = self._sdfits.index(bintable=self._bintable_index).iloc[self._refoffrows]["EXPOSURE"].to_numpy()
        exp_sig_on = self._sdfits.index(bintable=self._bintable_index).iloc[self._sigonrows]["EXPOSURE"].to_numpy()
        exp_sig_off = self._sdfits.index(bintable=self._bintable_index).iloc[self._sigoffrows]["EXPOSURE"].to_numpy()
        exp_ref = exp_ref_on + exp_ref_off
        exp_sig = exp_sig_on + exp_sig_off
        if self._smoothref > 1:
            nsmooth = self._smoothref
        else:
            nsmooth = 1.0
        exposure = exp_sig * exp_ref * nsmooth / (exp_sig + exp_ref * nsmooth)
        return exposure

    @property
    def delta_freq(self):
        """Get the array of channel frequency width

        df = [ 0.5*(df_ref_on + df_ref_off) + 0.5*(df_sig_on + df_sig_off) ] / 2

        Returns
        -------
             delta_freq: ~numpy.ndarray
                 The channel frequency width in units of the CDELT1 keyword in the SDFITS header
        """
        df_ref_on = self._sdfits.index(bintable=self._bintable_index).iloc[self._refonrows]["CDELT1"].to_numpy()
        df_ref_off = self._sdfits.index(bintable=self._bintable_index).iloc[self._refoffrows]["CDELT1"].to_numpy()
        df_sig_on = self._sdfits.index(bintable=self._bintable_index).iloc[self._sigonrows]["CDELT1"].to_numpy()
        df_sig_off = self._sdfits.index(bintable=self._bintable_index).iloc[self._sigoffrows]["CDELT1"].to_numpy()
        df_ref = 0.5 * (df_ref_on + df_ref_off)
        df_sig = 0.5 * (df_sig_on + df_sig_off)
        delta_freq = 0.5 * (df_ref + df_sig)
        return delta_freq

    @property
    def _tsys_weight(self):
        r"""The system temperature weighting array computed from current
        :math`T_{sys}`, :math:`t_{int}`, and :math:`\delta\nu`. See :meth:`tsys_weight`
        """
        return tsys_weight(self.exposure, self.delta_freq, self.tsys)

    def timeaverage(self, weights="tsys"):
        r"""Compute the time-averaged spectrum for this set of FSscans.

        Parameters
        ----------
        weights: str
            'tsys' or None.  If 'tsys' the weight will be calculated as:

             :math:`w = t_{exp} \times \delta\nu/T_{sys}^2`

            Default: 'tsys'
        Returns
        -------
        spectrum : :class:`~spectra.spectrum.Spectrum`
            The time-averaged spectrum
        """
        if self._calibrated is None or len(self._calibrated) == 0:
            raise Exception("You can't time average before calibration.")
        if self._npol > 1:
            raise Exception("Can't yet time average multiple polarizations %d" % self._npol)
        self._timeaveraged = deepcopy(self.calibrated(0))
        data = self._calibrated
        if weights == "tsys":
            w = self._tsys_weight
        else:
            w = np.ones_like(self._tsys_weight)
        self._timeaveraged._data = average(data, axis=0, weights=w)
        non_blanks = find_non_blanks(data)
        self._timeaveraged.meta["MEANTSYS"] = np.mean(self._tsys[non_blanks])
        self._timeaveraged.meta["WTTSYS"] = sq_weighted_avg(self._tsys[non_blanks], axis=0, weights=w[non_blanks])
        self._timeaveraged.meta["EXPOSURE"] = np.sum(self._exposure[non_blanks])
        self._timeaveraged.meta["TSYS"] = self._timeaveraged.meta["WTTSYS"]
        return self._timeaveraged


class SubBeamNodScan(ScanMixin):
    r"""
    Parameters
    ----------
    sigtp:  list of ~spectra.scan.TPScan
        Signal total power scans
    reftp:  list ~spectra.scan.TPScan
        Reference total power scans
    fulltp:  ~spectra.scan.TPScan
        A full (sig+ref) total power scans, used only for method='scan'
    method: str
        Method to use when processing. One of 'cycle' or 'scan'.  'cycle' is more accurate and averages data in each SUBREF_STATE cycle. 'scan' reproduces GBTIDL's snodka function which has been shown to be less accurate.  Default:'cycle'
    calibrate: bool
        Whether or not to calibrate the data.
    smoothref: int
        the number of channels in the reference to boxcar smooth prior to calibration
    observer_location : `~astropy.coordinates.EarthLocation`
        Location of the observatory. See `~dysh.coordinates.Observatory`.
        This will be transformed to `~astropy.coordinates.ITRS` using the time of
        observation DATE-OBS or MJD-OBS in
        the SDFITS header.  The default is the location of the GBT.
    weights: str
        Weighting scheme to use when averaging the signal and reference scans
        'tsys' or None.  If 'tsys' the weight will be calculated as:

         :math:`w = t_{exp} \times \delta\nu/T_{sys}^2`

        Default: 'tsys'
    """

    def __init__(
        self,
        sigtp,
        reftp,
        fulltp=None,
        method="cycle",
        calibrate=True,
        smoothref=1,
        observer_location=Observatory["GBT"],
        **kwargs,
    ):
        kwargs_opts = {
            "timeaverage": False,
            "weights": "tsys",  # or None or ndarray
            "debug": False,
        }
        kwargs_opts.update(kwargs)
        w = kwargs_opts["weights"]
        if len(reftp) != len(sigtp):
            raise ValueError(
                f"Reference and signal total power arrays are different lengths: {len(reftp)} != {len(sigtp)}"
            )
        self._scan = sigtp[0]._scan
        self._sigtp = sigtp
        self._reftp = reftp
        self._fulltp = fulltp
        self._ifnum = self._sigtp[0].ifnum
        self._fdnum = self._sigtp[0].fdnum
        self._npol = self._sigtp[0].npol
        self._pols = self._sigtp[0]._pols
        self._nchan = len(reftp[0]._data[0])
        self._nint = 0
        self._method = method.lower()
        if self._method not in ["cycle", "scan"]:
            raise ValueError(f"Method {self._method} unrecognized. Must be one of 'cycle' or 'scan'")
        self._smoothref = smoothref
        if self._smoothref > 1:
            print(f"SubBeamNodS smoothref={self._smoothref} not implemented yet")
        self._observer_location = observer_location
        self._calibrated = None
        if calibrate:
            self.calibrate(weights=w)

    def calibrate(self, **kwargs):
        """Calibrate the SUbBeamNodScan data"""
        nspect = len(self._reftp)
        self._tsys = np.empty(nspect, dtype=float)
        self._exposure = np.empty(nspect, dtype=float)
        self._delta_freq = np.empty(nspect, dtype=float)
        self._calibrated = np.empty((nspect, self._nchan), dtype=float)
        if self._method == "cycle":
            for i in range(nspect):
                ref_avg = self._reftp[i].timeaverage(weights=kwargs["weights"])
                sig_avg = self._sigtp[i].timeaverage(weights=kwargs["weights"])
                # Combine sig and ref.
                ta = ((sig_avg - ref_avg) / ref_avg).flux.value * ref_avg.meta["WTTSYS"]
                self._tsys[i] = ref_avg.meta["WTTSYS"]
                self._exposure[i] = sig_avg.meta["EXPOSURE"]
                self._delta_freq[i] = sig_avg.meta["CDELT1"]
                self._calibrated[i] = ta

        elif self._method == "scan":
            # Process the whole scan as a single block.
            # This is less accurate, but might be needed if
            # the scan was aborted and there are not enough
            # sig/ref cycles to do a per cycle calibration.
            for i in range(len(self._reftp)):
                on = self._sigtp[i].timeaverage(weights=kwargs["weights"]).data
                off = self._reftp[i].timeaverage(weights=kwargs["weights"]).data
                fulltpavg = self._fulltp[i].timeaverage(weights=kwargs["weights"])
                tsys = fulltpavg.meta["TSYS"]
                # data is a Spectrum.  not consistent with other Scan classes
                # where _calibrated is a numpy array.
                data = tsys * (on - off) / off
                # data.meta["MEANTSYS"] = 0.5 * np.mean((on.meta["TSYS"] + off.meta["TSYS"]))
                # data.meta["WTTSYS"] = tsys
                # data.meta["TSYS"] = data.meta["WTTSYS"]
                # self._tsys.append(ref_avg.meta["WTTSYS"])
                self._calibrated[i] = data
        else:
            raise ValueError(f"Method {self._method} unrecognized. Must be one of 'cycle' or 'scan'")
        # self._add_calibration_meta()

    def calibrated(self, i):
        meta = deepcopy(self._sigtp[i].timeaverage().meta)  # use self._sigtp.meta? instead?
        meta["TSYS"] = self._tsys[i]
        meta["EXPOSURE"] = self._exposure[i]
        naxis1 = len(self._calibrated[i])
        meta["NAXIS1"] = naxis1
        if "CUNIT1" not in meta:
            meta["CUNIT1"] = "Hz"  # @todo this is in gbtfits.hdu[0].header['TUNIT11'] but is it always TUNIT11?
        meta["CUNIT2"] = "deg"  # is this always true?
        meta["CUNIT3"] = "deg"  # is this always true?
        restfrq = meta["RESTFREQ"]
        rfq = restfrq * u.Unit(meta["CUNIT1"])
        restfreq = rfq.to("Hz").value
        meta["RESTFRQ"] = restfreq  # WCS wants no E
        return Spectrum.make_spectrum(self._calibrated[i] * u.K, meta, observer_location=self._observer_location)

    @property
    def exposure(self):
        return self._exposure

    @property
    def delta_freq(self):
        return self._delta_freq

    @property
    def tsys(self):
        return self._tsys

    @property
    def _tsys_weight(self):
        r"""The system temperature weighting array computed from current
        :math:`T_{sys}, t_{exp}`, and `\delta\nu`. See :meth:`tsys_weight`
        """
        return tsys_weight(self.exposure, self.delta_freq, self.tsys)

    def timeaverage(self, weights="tsys"):
        if self._calibrated is None or len(self._calibrated) == 0:
            raise Exception("You can't time average before calibration.")
        if self._npol > 1:
            raise Exception(f"Can't yet time average multiple polarizations {self._npol}")
        self._timeaveraged = deepcopy(self.calibrated(0))
        data = self._calibrated
        nchan = len(data[0])
        if weights == "tsys":
            w = self._tsys_weight
        else:
            w = None
        if self._method == "scan":
            w = None  # don't double tsys weight(?)
            self._timeaveraged = deepcopy(self.calibrated(0))
            self._timeaveraged._data = average(data, axis=0, weights=w)
            self._timeaveraged.meta["MEANTSYS"] = np.mean(self._tsys)
            # data.meta["MEANTSYS"] = 0.5 * np.mean((on.meta["TSYS"] + off.meta["TSYS"]))
            self._timeaveraged.meta["WTTSYS"] = self._timeaveraged.meta["WTTSYS"]
            self._timeaveraged.meta["TSYS"] = self._timeaveraged.meta["WTTSYS"]
            self._timeaveraged.meta["EXPOSURE"] = np.sum(self.exposure)
        if self._method == "cycle":  # or weights = "gbtidl"
            # GBTIDL method of weighting subbeamnod data
            ta_avg = np.zeros(nchan, dtype="d")
            wt_avg = 0.0  # A single value for now, but it should be an array once we implement vector TSYS.
            tsys_wt = 0.0
            tsys_avg = 0.0
            for i in range(len(data)):
                wt_avg += self.tsys[i] ** -2.0
                tsys_wt_ = tsys_weight(self.exposure[i], self.delta_freq[i], self.tsys[i])
                tsys_wt += tsys_wt_
                ta_avg[:] += data[i] * self.tsys[i] ** -2.0
            wt1 = self.tsys**-2.0
            wt2 = tsys_weight(self.exposure, self.delta_freq, self.tsys)
            ta_avg /= wt_avg
            tsys_avg /= tsys_wt
            self._timeaveraged._data = ta_avg
            self._timeaveraged.meta["MEANTSYS"] = np.mean(self._tsys)
            self._timeaveraged.meta["WTTSYS"] = tsys_avg  # sq_weighted_avg(self._tsys, axis=0, weights=w)
            self._timeaveraged.meta["TSYS"] = self._timeaveraged.meta["WTTSYS"]
            self._timeaveraged.meta["EXPOSURE"] = np.sum(self.exposure)

        return self._timeaveraged
