from collections import UserList
from copy import deepcopy

import astropy.units as u
import numpy as np
from astropy.wcs import WCS

from ..util import uniq
from . import (
    average,
    find_non_blanks,
    mean_tsys,
    sq_weighted_avg,
    tsys_weight,
    veldef_to_convention,
)
from .spectrum import Spectrum


class ScanMixin(object):
    @property
    def status(self):
        """Status flag, will be used later for undo"""
        return self._status

    @property
    def nchan(self):
        return self._nchan

    @property
    def nrows(self):
        """The number of rows in this Scan"""
        return self._nrows

    @property
    def npol(self):
        """The number of polarizations in this Scan"""
        return self._npol

    def nif(self):
        """The number of IFs in this Scan"""
        return self._nif

    def nfeed(self):
        """The number of feeds in this Scan"""
        return self._nfeed

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

    def __len__(self):
        return self._nrows


class TPScan(object):
    r"""
    Holds a total power scan.

    Parameters
    ----------
    sdfits : ~SDFITSLoad
        Input SDFITSLoad object (or derivative).
    scan: int
        Scan number.
    sigstate : str
        One of 'SIG' or 'REF' to indicate if this is the signal or reference scan.
    scanrows : list-like
        The list of rows in `sdfits` corresponding to sig_state integrations.
    bintable : int
        The index for BINTABLE in `sdfits` containing the scans.
    """

    def __init__(self, sdfits, scan, sigstate, calstate, scanrows, bintable):
        self._sdfits = sdfits  # parent class
        self._scan = scan
        self._sigstate = sigstate
        self._calstate = calstate
        self._scanrows = scanrows
        self._bintable_index = bintable
        self._data = self._sdfits.rawspectra(bintable)[scanrows]  # all cal states
        self._status = 0  # @TODO make these an enumeration, possibly dict
        #                           # ex1:
        self._nint = 0  # 11
        self._npol = 0  #  2
        self._timeaveraged = None  #  2
        self._polaveraged = None  #  1
        self._nrows = len(scanrows)
        self._tsys = None
        print(f"TPSCAN nrows = {self.nrows}")


class ScanBlock(UserList, ScanMixin):
    def __init__(self, *args):
        super().__init__(*args)
        self._status = None
        self._nrows = 0
        self._npol = 0
        self._timeaveraged = []
        self._polaveraged = []
        self._finalspectrum = []

    def calibrate(self, **kwargs):
        """Calibrate all scans in this ScanBlock"""
        for scan in self.data:
            scan.calibrate(**kwargs)

    def timeaverage(self, weights="tsys"):
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
        for scan in self.data:
            self._timeaveraged.append(scan.timeaverage(weights))
        return self._timeaveraged

    def polaverage(self, weights="tsys"):
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
        for scan in self.data:
            self._finalspectrum.append(scan.finalspectrum(weights))
        return self._finalspectrum


class TPScan(ScanMixin):
    """GBT specific version of Total Power Scan

    Parameters
    ----------
    gbtfits : `~fits.gbtfitsload.GBFITSLoad`
        input GBFITSLoad object
    scan: int
        scan number
    sigstate : str
        one of 'SIG' or 'REF' to indicate if this is the signal or reference scan or 'BOTH' if it contains both
    calstate : str
        one of 'ON' or 'OFF' to indicate the calibration state of this scan, or 'BOTH' if it contains both
    scanrows : list-like
        the list of rows in `sdfits` corresponding to sig_state integrations
    calrows : dict
        dictionary containing with keys 'ON' and 'OFF' containing list of rows in `sdfits` corresponding to cal=T (ON) and cal=F (OFF) integrations for `scan`
    bintable : int
        the index for BINTABLE in `sdfits` containing the scans
    calibrate: bool
        whether or not to calibrate the data.  If `True`, the data will be (calon - caloff)*0.5, otherwise it will be SDFITS row data. Default:True
    """

    # @TODO get rid of calrows and calc tsys in gettp and pass it in.
    def __init__(self, gbtfits, scan, sigstate, calstate, scanrows, calrows, bintable, calibrate=True):
        self._sdfits = gbtfits  # parent class
        self._scan = scan
        self._sigstate = sigstate  # ignored?
        self._calstate = calstate  # ignored?
        self._scanrows = scanrows
        # print("BINTABLE = ", bintable)
        # @TODO deal with data that crosses bintables
        if bintable is None:
            self._bintable_index = self._sdfits._find_bintable_and_row(self._scanrows[0])[0]
        else:
            self._bintable_index = bintable
        self._data = self._sdfits.rawspectra(self._bintable_index)[scanrows]  # all cal states
        df = self._sdfits._index
        df = df.iloc[scanrows]
        self._index = df
        self._feeds = uniq(df["FDNUM"])
        self._pols = uniq(df["PLNUM"])
        self._ifs = uniq(df["IFNUM"])
        self._status = 0  # @TODO make these an enumeration, possibly dict
        self._nint = 0
        self._npol = len(self._pols)
        self._nfeed = len(self._feeds)
        self._nif = len(self._ifs)
        self._timeaveraged = None
        self._polaveraged = None
        self._nrows = len(scanrows)
        self._tsys = None
        if False:
            self._npol = gbtfits.npol(self._bintable_index)  # TODO deal with bintable
            self._nint = gbtfits.nintegrations(self._bintable_index)
        self._calrows = calrows
        self._refonrows = self._calrows["ON"]
        self._refoffrows = self._calrows["OFF"]
        self._refcalon = gbtfits.rawspectra(self._bintable_index)[self._refonrows]
        self._refcaloff = gbtfits.rawspectra(self._bintable_index)[self._refoffrows]
        self._nchan = len(self._refcalon[0])
        self._calibrate = calibrate
        if self._calibrate:
            self._data = (0.5 * (self._refcalon + self._refcaloff)).astype(float)
        # print(f"# scanrows {len(self._scanrows)}, # calrows ON {len(self._calrows['ON'])}  # calrows OFF {len(self._calrows['OFF'])}")
        self.calc_tsys()

    @property
    def sigstate(self):
        return self._sigstate

    @property
    def calstate(self):
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
        Calculate the system temperature array
        """
        kwargs_opts = {"verbose": False}
        kwargs_opts.update(kwargs)

        self._status = 1

        tcal = list(self._sdfits.index(bintable=self._bintable_index).iloc[self._refonrows]["TCAL"])
        nspect = len(tcal)
        self._tsys = np.empty(nspect, dtype=float)  # should be same as len(calon)
        # allcal = self._refonrows.copy()
        # allcal.extend(self._refoffrows)
        # tcal = list(self._sdfits.index(self._bintable_index).iloc[sorted(allcal)]["TCAL"])
        # @Todo  this loop could be replaces with clever numpy
        if len(tcal) != nspect:
            raise Exception(f"TCAL length {len(tcal)} and number of spectra {nspect} don't match")
        for i in range(nspect):
            tsys = mean_tsys(calon=self._refcalon[i], caloff=self._refcaloff[i], tcal=tcal[i])
            self._tsys[i] = tsys

    @property
    def exposure(self):
        """Get the array of exposure (integration) times

            exposure =  0.5*(exp_ref_on + exp_ref_off)

        Note we only have access to the refon and refoff row indices so
        can't use sig here.  This is probably incorrect

        Returns
        -------
        exposure : `~numpy.ndarray`
            The exposure time in units of the EXPOSURE keyword in the SDFITS header
        """
        exp_ref_on = self._sdfits.index(bintable=self._bintable_index).iloc[self._refonrows]["EXPOSURE"].to_numpy()
        exp_ref_off = self._sdfits.index(bintable=self._bintable_index).iloc[self._refoffrows]["EXPOSURE"].to_numpy()
        exposure = exp_ref_on + exp_ref_off
        return exposure

    @property
    def delta_freq(self):
        """Get the array of channel frequency width

           df =  0.5*(df_ref_on + df_ref_off)

        Note we only have access to the refon and refoff row indices so
        can't use sig here.  This is probably incorrect

        Returns
        -------
        delta_freq: `~numpy.ndarray`
            The channel frequency width in units of the CDELT1 keyword in the SDFITS header
        """
        df_ref_on = self._sdfits.index(bintable=self._bintable_index).iloc[self._refonrows]["CDELT1"].to_numpy()
        df_ref_off = self._sdfits.index(bintable=self._bintable_index).iloc[self._refoffrows]["CDELT1"].to_numpy()
        delta_freq = 0.5 * (df_ref_on + df_ref_off)
        return delta_freq

    @property
    def _tsys_weight(self):
        r"""The system temperature weighting array computed from current
        :math:`T_{sys}, t_{exp}`, and `\delta\nu`. See :meth:`tsys_weight`
        """
        return tsys_weight(self.exposure, self.delta_freq, self.tsys)

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
        # print(len(self._scanrows), i)
        meta = dict(self._sdfits.index(bintable=self._bintable_index).iloc[self._scanrows[i]])
        meta["TSYS"] = self._tsys[i]
        meta["EXPOSURE"] = self.exposure[i]
        naxis1 = len(self._data[i])
        meta["CTYPE1"]
        ctype2 = meta["CTYPE2"]
        ctype3 = meta["CTYPE3"]
        crval1 = meta["CRVAL1"]
        crval2 = meta["CRVAL2"]
        crval3 = meta["CRVAL3"]
        crpix1 = meta["CRPIX1"]
        cdelt1 = meta["CDELT1"]
        restfrq = meta["RESTFREQ"]
        if "CUNIT1" in meta:
            cunit1 = meta["CUNIT1"]
        else:
            cunit1 = "Hz"  # @TODO this is in gbtfits.hdu[0].header['TUNIT11'] but is it always TUNIT11?
        rfq = restfrq * u.Unit(cunit1)
        restfreq = rfq.to("Hz").value

        # @TODO WCS is expensive.  Figure how to calculate spectral_axis instead.
        wcs = WCS(
            header={
                "CDELT1": cdelt1,
                "CRVAL1": crval1,
                "CUNIT1": cunit1,
                "CTYPE1": "FREQ",
                "CRPIX1": crpix1,
                "RESTFRQ": restfreq,
                "CTYPE2": ctype2,
                "CRVAL2": crval2,
                "CRPIX2": 1,
                "CTYPE3": ctype3,
                "CRVAL3": crval3,
                "CRPIX3": 1,
                "CUNIT2": "deg",
                "CUNIT3": "deg",
                "NAXIS1": naxis1,
                "NAXIS2": 1,
                "NAXIS3": 1,
            },
        )
        vc = veldef_to_convention(meta["VELDEF"])
        s = Spectrum(self._data[i] * u.ct, wcs=wcs, meta=meta, velocity_convention=vc)
        return s

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
        if self._nif > 1:
            raise Exception("Can't yet time average multiple IFs")
        if self._nfeed > 1:
            raise Exception("Can't yet time average multiple feeds")
        self._timeaveraged = deepcopy(self.total_power(0))
        if weights == "tsys":
            w = self._tsys_weight
        else:
            w = np.ones_like(self._tsys_weight)
        non_blanks = find_non_blanks(self._data)
        self._timeaveraged._data = average(self._data, axis=0, weights=w)
        self._timeaveraged.meta["MEANTSYS"] = np.mean(self._tsys[non_blanks])
        self._timeaveraged.meta["WTTSYS"] = sq_weighted_avg(self._tsys[non_blanks], axis=0, weights=w[non_blanks])
        self._timeaveraged.meta["TSYS"] = self._timeaveraged.meta["WTTSYS"]
        self._timeaveraged.meta["EXPOSURE"] = self.exposure[non_blanks].sum()
        return self._timeaveraged


class PSScan(ScanMixin):
    """GBT specific version of Position Switch Scan

    Parameters
    ----------

    gbtfits : `~fit.gbtfitsload.GBFITSLoad`
        input GBFITSLoad object
    scans : dict
        dictionary with keys 'ON' and 'OFF' containing unique list of ON (signal) and OFF (reference) scan numbers
    scanrows : dict
        dictionary with keys 'ON' and 'OFF' containing the list of rows in `sdfits` corresponding to ON (signal) and OFF (reference) integrations
    calrows : dict
        dictionary containing with keys 'ON' and 'OFF' containing list of rows in `sdfits` corresponding to cal=T (ON) and cal=F (OFF) integrations.
    bintable : int
        the index for BINTABLE in `sdfits` containing the scans
    calibrate: bool
        whether or not to calibrate the data.  If true, data will be calibrated as TSYS*(ON-OFF)/OFF. Default: True
    """

    def __init__(self, gbtfits, scans, scanrows, calrows, bintable, calibrate=True):
        # The rows of the original bintable corresponding to ON (sig) and OFF (reg)
        self._sdfits = gbtfits  # parent class
        self._scans = scans
        self._scanrows = scanrows
        self._nrows = len(self._scanrows["ON"])
        # print(f"len(scanrows ON) {len(self._scanrows['ON'])}")
        # print(f"len(scanrows OFF) {len(self._scanrows['OFF'])}")

        # calrows perhaps not needed as input since we can get it from gbtfits object?
        # calrows['ON'] are rows with noise diode was on, regardless of sig or ref
        # calrows['OFF'] are rows with noise diode was off, regardless of sig or ref
        self._calrows = calrows
        # print("BINTABLE = ", bintable)
        # @TODO deal with data that crosses bintables
        if bintable is None:
            self._bintable_index = gbtfits._find_bintable_and_row(self._scanrows["ON"][0])[0]
        else:
            self._bintable_index = bintable
        df = self._sdfits._index
        df = df.iloc[scanrows["ON"]]
        self._feeds = uniq(df["FDNUM"])
        self._pols = uniq(df["PLNUM"])
        self._ifs = uniq(df["IFNUM"])
        self._npol = len(self._pols)
        self._nfeed = len(self._feeds)
        self._nif = len(self._ifs)
        if False:
            self._nint = gbtfits.nintegrations(self._bintable_index)
        # todo use gbtfits.velocity_convention(veldef,velframe)
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
        if self._calibrate:
            self.calibrate()

    @property
    def tsys(self):
        """The system temperature array. This will be `None` until calibration is done.

        Returns
        -------
        tsys : `~numpy.ndarray`
            System temperature values in K
        """
        return self._tsys

    # TODO something clever
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
        meta = dict(self._sdfits.index(bintable=self._bintable_index).iloc[self._scanrows["ON"][i]])
        meta["TSYS"] = self._tsys[i]
        meta["EXPOSURE"] = self._exposure[i]
        naxis1 = len(self._calibrated[i])
        meta["CTYPE1"]
        ctype2 = meta["CTYPE2"]
        ctype3 = meta["CTYPE3"]
        crval1 = meta["CRVAL1"]
        crval2 = meta["CRVAL2"]
        crval3 = meta["CRVAL3"]
        crpix1 = meta["CRPIX1"]
        cdelt1 = meta["CDELT1"]
        restfrq = meta["RESTFREQ"]
        if "CUNIT1" in meta:
            cunit1 = meta["CUNIT1"]
        else:
            cunit1 = "Hz"  # @TODO this is in gbtfits.hdu[0].header['TUNIT11'] but is it always TUNIT11?
        rfq = restfrq * u.Unit(cunit1)
        restfreq = rfq.to("Hz").value

        # @TODO WCS is expensive.  Figure how to calculate spectral_axis instead.
        wcs = WCS(
            header={
                "CDELT1": cdelt1,
                "CRVAL1": crval1,
                "CUNIT1": cunit1,
                "CTYPE1": "FREQ",
                "CRPIX1": crpix1,
                "RESTFRQ": restfreq,
                "CTYPE2": ctype2,
                "CRVAL2": crval2,
                "CRPIX2": 1,
                "CTYPE3": ctype3,
                "CRVAL3": crval3,
                "CRPIX3": 1,
                "CUNIT2": "deg",
                "CUNIT3": "deg",
                "NAXIS1": naxis1,
                "NAXIS2": 1,
                "NAXIS3": 1,
            },
        )
        vc = veldef_to_convention(meta["VELDEF"])

        return Spectrum(self._calibrated[i] * u.K, wcs=wcs, meta=meta, velocity_convention=vc)

    def calibrate(self, **kwargs):
        """
        Position switch calibration, following equations 1 and 2 in the GBTIDL calibration manual
        """
        kwargs_opts = {"verbose": False}
        kwargs_opts.update(kwargs)

        self._status = 1
        nspect = self.nrows // 2
        self._calibrated = np.empty((nspect, self._nchan), dtype="d")
        self._tsys = np.empty(nspect, dtype="d")
        self._exposure = np.empty(nspect, dtype="d")

        tcal = list(self._sdfits.index(bintable=self._bintable_index).iloc[self._refonrows]["TCAL"])
        # @Todo  this loop could be replaced with clever numpy
        if len(tcal) != nspect:
            raise Exception(f"TCAL length {len(tcal)} and number of spectra {nspect} don't match")
        for i in range(nspect):
            tsys = mean_tsys(calon=self._refcalon[i], caloff=self._refcaloff[i], tcal=tcal[i])
            sig = 0.5 * (self._sigcalon[i] + self._sigcaloff[i])
            ref = 0.5 * (self._refcalon[i] + self._refcaloff[i])
            self._calibrated[i] = tsys * (sig - ref) / ref
            self._tsys[i] = tsys
            self._exposure[i] = self.exposure[i]

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
        # exposure = 0.5*(exp_ref + exp_sig)
        # exposure = exp_ref + exp_sig
        nsmooth = 1.0  # In case we start smoothing the reference spectra.
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
        if self._nif > 1:
            raise Exception("Can't yet time average multiple IFs")
        if self._nfeed > 1:
            raise Exception("Can't yet time average multiple feeds")
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


class SubBeamNodScan(ScanMixin):  # SBNodScan?
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
    weights: str
        Weighting scheme to use when averaging the signal and reference scans
        'tsys' or None.  If 'tsys' the weight will be calculated as:

         :math:`w = t_{exp} \times \delta\nu/T_{sys}^2`

        Default: 'tsys'
    """

    def __init__(self, sigtp, reftp, fulltp=None, method="cycle", calibrate=True, **kwargs):
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
        self._sigtp = sigtp
        self._reftp = reftp
        self._fulltp = fulltp
        self._nchan = len(reftp[0]._data[0])
        self._npol = 1
        self._nint = 0
        self._nif = 1
        self._nfeed = 1
        self._method = method.lower()
        if self._method not in ["cycle", "scan"]:
            raise ValueError(f"Method {self._method} unrecognized. Must be one of 'cycle' or 'scan'")
        self._calibrated = None
        if calibrate:
            self.calibrate(weights=w)

    def calibrate(self, **kwargs):
        """Calibrate the Scan data"""
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

    def calibrated(self, i):
        meta = deepcopy(self._sigtp[i].timeaverage().meta)
        meta["TSYS"] = self._tsys[i]
        naxis1 = len(self._calibrated[i])
        meta["CTYPE1"]
        ctype2 = meta["CTYPE2"]
        ctype3 = meta["CTYPE3"]
        crval1 = meta["CRVAL1"]
        crval2 = meta["CRVAL2"]
        crval3 = meta["CRVAL3"]
        crpix1 = meta["CRPIX1"]
        cdelt1 = meta["CDELT1"]
        restfrq = meta["RESTFREQ"]
        if "CUNIT1" in meta:
            cunit1 = meta["CUNIT1"]
        else:
            cunit1 = "Hz"  # @TODO this is in gbtfits.hdu[0].header['TUNIT11'] but is it always TUNIT11?
        rfq = restfrq * u.Unit(cunit1)
        restfreq = rfq.to("Hz").value

        # @TODO WCS is expensive.  Figure how to calculate spectral_axis instead.
        wcs = WCS(
            header={
                "CDELT1": cdelt1,
                "CRVAL1": crval1,
                "CUNIT1": cunit1,
                "CTYPE1": "FREQ",
                "CRPIX1": crpix1,
                "RESTFRQ": restfreq,
                "CTYPE2": ctype2,
                "CRVAL2": crval2,
                "CRPIX2": 1,
                "CTYPE3": ctype3,
                "CRVAL3": crval3,
                "CRPIX3": 1,
                "CUNIT2": "deg",
                "CUNIT3": "deg",
                "NAXIS1": naxis1,
                "NAXIS2": 1,
                "NAXIS3": 1,
            },
        )
        vc = veldef_to_convention(meta["VELDEF"])

        return Spectrum(self._calibrated[i] * u.K, wcs=wcs, meta=meta, velocity_convention=vc)

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
            raise Exception("Can't yet time average multiple polarizations")
        if self._nif > 1:
            raise Exception("Can't yet time average multiple IFs")
        if self._nfeed > 1:
            raise Exception("Can't yet time average multiple feeds")
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
