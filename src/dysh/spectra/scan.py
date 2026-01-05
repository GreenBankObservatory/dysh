"""
The classes that define various types of Scan and their calibration methods.
"""

import warnings
from abc import abstractmethod
from collections import UserList
from copy import deepcopy

import astropy.units as u
import numpy as np
from astropy import constants as ac
from astropy.io.fits import BinTableHDU, Column
from astropy.table import Table, vstack
from astropy.time import Time
from astropy.utils.masked import Masked

from dysh.spectra import core

from ..coordinates import Observatory
from ..log import HistoricalBase, log_call_to_history, logger
from ..plot import scanplot as sp
from ..util import isot_to_mjd, minimum_string_match
from ..util.docstring_manip import copy_docstring
from ..util.gaincorrection import GBTGainCorrection
from .core import (
    available_smooth_methods,
    find_non_blanks,
    find_nonblank_ints,
    mean_tsys,
    smooth,
    sq_weighted_avg,
    tsys_weight,
)
from .spectrum import Spectrum, average_spectra
from .vane import VaneSpectrum


class SpectralAverageMixin:
    @log_call_to_history
    def smooth(self, method="hanning", width=1, decimate=0):
        """
        Smooth or convolve the underlying calibrated data array, optionally decimating the data.

        A number of methods from astropy.convolution can be selected
        with the `method` keyword.

        Default smoothing is hanning.

        Note: Any previously computed/removed baseline will remain unchanged.

        Parameters
        ----------.
        method : {'hanning', 'boxcar', 'gaussian'}, optional
            Smoothing method. Valid are: 'hanning', 'boxcar' and
            'gaussian'. Minimum match applies.
            The default is 'hanning'.
        width : int, optional
            Effective width of the convolving kernel.  Should ideally be an
            odd number.
            For 'hanning' this should be 1, with a 0.25,0.5,0.25 kernel.
            For 'boxcar' an even value triggers an odd one with half the
            signal at the edges, and will thus not reproduce GBTIDL.
            For 'gaussian' this is the FWHM of the final beam. We normally
            assume the input beam has FWHM=1, pending resolution on cases
            where CDELT1 is not the same as FREQRES.
            The default is 1.
        decimate : int, optional
            Decimation factor of the spectrum by returning every decimate channel.
            -1:   no decimation
            0:    use the width parameter
            >1:   user supplied decimation (use with caution)
            The default is 0, meaning decimation is by `width`

        Returns
        -------
        None

        Raises
        ------
        ValueError
            If no valid smoothing method is given.
        """

        valid_methods = available_smooth_methods()
        this_method = minimum_string_match(method, valid_methods)
        if width < 1:
            raise ValueError(f"`width` ({width}) must be >=1.")

        if this_method is None:
            raise ValueError(f"Unrecognized method ({method}). Valid methods are {valid_methods}")
        if decimate == 0:
            # Take the default decimation by `width`.
            decimate = int(abs(width))
            if not float(width).is_integer():
                logger.info(f"Adjusting decimation factor to be a natural number. Will decimate by {decimate}")
        clen, _nchan = self._calibrated.shape
        sdata = []
        smask = []
        meta = []
        for i in range(clen):
            c = self._calibrated[i]
            newdata, newmeta = smooth(
                data=c,
                method=method,
                width=width,
                ndecimate=decimate,
                kernel=None,
                meta=self.meta[i],
            )
            if hasattr(newdata, "mask"):
                smask.append(newdata.mask)
            sdata.append(newdata)
            meta.append(newmeta)
        self._calibrated = np.ma.masked_array(sdata, smask)
        self._meta = meta
        # If decimation occurs we must recompute delta_freq.
        if decimate > -1:
            self._set_delta_freq_from_meta()

    def _set_delta_freq_from_meta(self):
        """After decimation reset the delta_freq variable from the recomputed metadata."""
        self._delta_freq = np.array([x["CDELT1"] for x in self._meta])

    @abstractmethod
    def _calc_delta_freq(self):
        """

        Calculate the channel frequency spacing.

        Returns
        -------
        None.

        """
        pass

    @log_call_to_history
    def timeaverage(self, weights: str | np.ndarray = "tsys", use_wcs=True) -> Spectrum:
        r"""
        Compute the time-averaged spectrum. For a Scan this will average all the integrations in the Scan,
        according to the given weights.
        For a `~dysh.spectra.scan.ScanBlock`, it will average all the integrations in all Scans contained in the `~dysh.spectra.scan.ScanBlock`.

        Parameters
        ----------
        weights: None, str or `~numpy.ndarray`
            If None, the channel weights will be equal and set to unity.

            If 'tsys' the channel weights will be calculated as:

             :math:`w = t_{exp} \times \delta\nu/T_{sys}^2`

            If an array, it must have shape `(Nint,)` or `(Nint,nchan)` where `Nint` is the number
            of integrations in the Scan or ScanBlock and `nchan` is the number of channels in the Scan.

        use_wcs : bool
            Create a WCS object for the resulting `~dysh.spectra.spectrum.Spectrum`.
            Creating a WCS object adds computation time to the creation of a `~dysh.spectra.spectrum.Spectrum` object.
            There may be cases where the WCS is not needed, so setting this boolean to False will save computation.

        Returns
        -------
        spectrum : `~dysh.spectra.spectrum.Spectrum`
            The time-averaged spectrum. The weights array of the Spectrum will have shape `(nchan,)` and will be set to the sum of the input weights.

        .. note::

           Data that are masked will have values set to zero.  This is a feature of `numpy.ma.average`. Data mask fill value is NaN (np.nan)

        """
        pass

    @property
    def exposure(self) -> np.ndarray:
        """
        The array of exposure (integration) times. How the exposure is calculated
        varies for different derived classes.  See :meth:`_calc_exposure`.

        Returns
        -------
        exposure : `~numpy.ndarray`
            The exposure time in units of the EXPOSURE keyword in the SDFITS header.
        """
        return self._exposure

    @property
    def duration(self) -> np.ndarray:
        """
        The array of duration times. How the duration is calculated
        varies for different derived classes.  See :meth:`_calc_exposure`.
        Duration includes the blanking time, exposure does not, so `duration>=exposure`.

        Returns
        -------
        duration : `~numpy.ndarray`
            The duration time in units of the DURATION keyword in the SDFITS header.
        """
        return self._duration

    @property
    def delta_freq(self) -> np.ndarray:
        """
        The array of channel frequency width.
        How the channel width is calculated varies for different derived classes. See :meth:`_calc_delta_freq`.

        Returns
        -------
        delta_freq: `~numpy.ndarray`
            The channel frequency width in units of the CDELT1 keyword in the SDFITS header.

        """
        return self._delta_freq

    @property
    def tsys(self):
        """
        The system temperature array.

        Returns
        -------
        tsys : `~numpy.ndarray`
            System temperature values in K.
        """
        return self._tsys

    @property
    def tsys_weight(self) -> np.ndarray:
        r"""
        The system temperature weighting array computed from current
        :math:`T_{sys}`, :math:`t_{int}`, and :math:`\delta\nu`. See :meth:`tsys_weight`
        """
        return tsys_weight(self.exposure, self.delta_freq, self.tsys)

    @property
    def weights(self) -> np.ndarray:
        """
        The weights for each integration after an averaging operation.  If `Scan.timeaverage()` has not been
        called, the weights will be unity.  The weights array can have shape `(nint,)` or `(nint,nchan)` depending on
        how `timeaverage()` was called.
        """
        return self._weights


class ScanBase(HistoricalBase, SpectralAverageMixin):
    """This class describes the common interface to all Scan classes.
    A Scan represents one scan number, one IF, one feed, and one polarization.
    Derived classes *must* implement :meth:`calibrate`.
    """

    def __init__(
        self,
        sdfits,
        smoothref,
        apply_flags,
        observer_location,
        fdnum=-1,
        ifnum=-1,
        plnum=-1,
        tsys=None,
        tcal=None,
        ap_eff=None,
        surface_error=None,
        zenith_opacity=None,
        channel=None,
    ):
        HistoricalBase.__init__(self)
        self._fdnum = fdnum
        self._ifnum = ifnum
        self._plnum = plnum
        self._nchan = -1
        self._scan = -1
        self._nrows = -1
        self._bintable_index = -1
        self._pols = ""  # currently unused. Will contain polarization stokes string.
        self._nint = -1
        self._sdfits = sdfits
        self._nocal = False
        self._meta = {}
        self._tscale_fac = np.array([1.0])
        self._tscale = "ta"
        self._tsys = tsys
        self._tcal = tcal
        self._ap_eff = ap_eff
        self._ap_eff_array = None
        self._surface_error = surface_error
        self._surface_error_array = None
        self._zenith_opacity = zenith_opacity
        self._exposure = None
        self._duration = None
        self._calibrated = None
        self._smoothref = smoothref
        self._apply_flags = apply_flags
        self._observer_location = observer_location
        self._tscale_to_unit = {"ta": u.K, "ta*": u.K, "flux": u.Jy, "raw": u.ct, "counts": u.ct, "count": u.ct}
        if channel is not None:
            self._channel_slice = slice(channel[0], channel[1])
        else:
            self._channel_slice = slice(0, None)
        # @todo Baseline fitting of scanblock. See issue (RFE) #607 https://github.com/GreenBankObservatory/dysh/issues/607
        self._baseline_model = None
        self._subtracted = False  # This is False if and only if baseline_model is None so we technically don't need a separate boolean.
        self._plotter = None
        self._check_gain_factors(self._ap_eff, self._surface_error)

    def _validate_defaults(self):
        _required = {
            "bintable_index": self._bintable_index,
            "ifnum": self._ifnum,
            "fdnum": self._fdnum,
            "nchan": self._nchan,
            "nrows": self._nrows,
            "plnum": self._plnum,
            "scan": self._scan,
        }
        unset = []
        for k, v in _required.items():
            if v == -1:
                unset.append(k)
        if len(unset) > 0:
            raise Exception(
                f"The following required Scan attributes were not set by the derived class {self.__class__.__name__}:"
                f" {unset}"
            )
        for k in ["ifnum", "plnum", "fdnum"]:
            v = _required[k]
            if not isinstance(v, (int, np.integer)):
                raise ValueError(f"{self.__class__.__name__}: {k} must be an integer but got {type(v)}")

    @classmethod
    def _check_tscale(self, tscale):
        """
        Check that the requested brightness scale is valid.
        This allows us to not import `GBTGainCorretion` into `GBTFITSLoad`.

        Parameters
        ----------
        tscale : str
            Strings representing valid options for scaling spectral data, specifically
                - 'Raw' : raw value, e.g., count
                - 'Ta'  : Antenna Temperature
                - 'Ta*' : Antenna temperature corrected to above the atmosphere
                - 'Flux'  : flux density in Jansky
            This parameter is case-insensitive.

        Raises
        ------
        ValueError
            If the scale is unrecognized.
        """
        if not GBTGainCorrection.is_valid_scale(tscale):
            raise ValueError(
                f"Unrecognized brightness scale {tscale}. Valid options are {GBTGainCorrection.valid_scales} (case-insensitive)."
            )

    @classmethod
    def _check_gain_factors(self, ap_eff, surface_error):
        if ap_eff is not None and surface_error is not None:
            raise ValueError("Only one of ap_eff or surface_error should be specified")

    def _finish_initialization(
        self, calibrate, calibrate_kwargs, meta_rows, tscale, zenith_opacity, tsys=None, tcal=None
    ):
        if len(meta_rows) == 0:
            raise Exception(
                f"In Scan {self.scan}, no data left to calibrate. Check blank integrations, flags, and selection."
            )
        self._calibrate = calibrate

        self._nint = len(meta_rows)
        self._make_meta(meta_rows)
        self._tscale_fac = np.ones(self._nint)
        self._init_tsys(tsys)
        self._init_tcal(tcal)
        self._weights = np.ones(self.nint)
        self._calc_exposure()
        self._calc_delta_freq()
        if self._calibrate:
            if calibrate_kwargs is not None:
                self.calibrate(**calibrate_kwargs)
            else:
                self.calibrate()
            self._add_calibration_meta()
        if zenith_opacity is not None or self._ap_eff is not None:
            self.scale(tscale, zenith_opacity)
        self._update_scale_meta()
        self._validate_defaults()

    @abstractmethod
    def _calc_exposure(self):
        """Method to compute the specific exposure array for the given Scan type"""
        raise NotImplementedError(f"Exposure calculation for {self.__class__.__name__} needs to be implemented.")
        # actually you won't even be able to instantiate the class if the method is not implemented.

    @abstractmethod
    def _calc_delta_freq(self):
        """Method to compute the specific channel width array for the given Scan type"""
        raise NotImplementedError(
            f"Delta Freq (channel width) calculation for {self.__class__.__name__} needs to be implemented."
        )

    def getspec(self, i: int, use_wcs: bool = True) -> Spectrum:  ##SCANBASE
        """Return the i-th calibrated Spectrum from this Scan.

        Parameters
        ----------
        i : int
            The index into the calibrated array
        use_wcs : bool
            Create a WCS object for the resulting `~dysh.spectra.spectrum.Spectrum`.
            Creating a WCS object adds computation time to the creation of a `~dysh.spectra.spectrum.Spectrum` object.
            There may be cases where the WCS is not needed, so setting this boolean to False will save computation.

        Returns
        -------
        spectrum : `~dysh.spectra.spectrum.Spectrum`
        """
        s = Spectrum.make_spectrum(
            Masked(
                self._calibrated[i] * self._tscale_to_unit[self.tscale.lower()],
                self._calibrated[i].mask,
            ),
            meta=self.meta[i],
            observer_location=self._observer_location,
            use_wcs=use_wcs,
        )
        s.merge_commentary(self)
        s._baseline_model = self._baseline_model
        s._subtracted = self._subtracted
        return s

    @property
    def is_scaled(self) -> bool:
        r"""Is this Scan scaled to something other than antenna temperature :math:`T_A`.

        Returns
        -------
        is_scaled : bool
            True if scale is e.g :math:`T_A^*` or Flux (:math:`S_\nu`).
        """
        return self._tscale.lower() != "ta"

    @property
    def tscale(self) -> str:
        """
        The descriptive brightness unit of the data. One of
            - 'Raw' : raw value, e.g., count
            - 'Ta'  : Antenna Temperature in K
            - 'Ta*' : Antenna temperature corrected to above the atmosphere in K
            - 'Flux': flux density in Jansky

        Returns
        -------
        tscale : str
            Brightness unit string.

        """
        return self._tscale[0].upper() + self._tscale[1:]

    @property
    def tscale_fac(self) -> np.ndarray:
        """
        The factor(s) by which the data have been scale from antenna temperature to corrected antenna temperature
        or flux density.

        Returns
        -------
        tscale_fac : `~numpy.ndarray`
            An array of floats, one per integration in the scan.

        """
        return self._tscale_fac

    @property
    def tunit(self) -> u.Unit:
        """The brightness unit (temperature or flux density)  of this Scan's data

        Returns
        -------
        tunit : `~astropy.units.Unit`
            The brightness unit
        """
        return self._tscale_to_unit[self._tscale.lower()]

    def _scaleby(self, factor):
        """Scale the calibrated data array by a factor. This is an NxM * N multiplication

        Parameters
        ----------
            factor - `~numpy.ndarray` or float

            The factor to scale the spectral data by

        Returns
        -------
            None
        """
        # Array multiplication is slightly faster than a loop
        self._calibrated = self._calibrated * (np.array([factor]).T)

    @log_call_to_history
    def scale(self, tscale: str, zenith_opacity: float | None = None):
        """
        Scale the data to the given brightness scale (temperature or flux density) using the zenith opacity.
        If data are already scaled, they will be unscaled first.

        Parameters
        ----------
        tscale : str
            Strings representing valid options for scaling spectral data, specifically
                - 'Ta'  : Antenna Temperature in K
                - 'Ta*' : Antenna temperature corrected to above the atmosphere in K
                - 'Flux': flux density in Jansky
            This parameter is case-insensitive.

        zenith_opacity : float, optional
            The zenith opacity. Required if `tscale` is 'Ta*' or 'Flux'.

        Returns
        -------
        None

        Raises
        ------
        TypeError
            If scaling to the `tscale` unit is not applicable to the scan type, e.g., a total power scan.
        ValueError
            If `tscale` is unrecognized or `zenith_opacity` is negative.

        """
        if self.__class__ == TPScan:
            raise TypeError("Total power data cannot be directly scaled to temperature or flux density.")
        self._check_tscale(tscale)
        s = tscale.lower()
        if s == self._tscale.lower():
            # requested scale is the current scale, nothing to be done.
            return
        if zenith_opacity is not None and zenith_opacity < 0:
            raise ValueError("Zenith opacity cannot be negative.")
        if s != "ta" and zenith_opacity is None:
            # and self._ap_eff is not None:
            raise ValueError("Zenith opacity must be provided when scaling to Ta* or Flux.")
        ntscale = self._tscale_to_unit[s].to_string()
        # unscale the data if it was already scaled.
        if self.is_scaled:
            self._scaleby(1.0 / self._tscale_fac)

        # if scaling back to antenna temperature, reset the scale factor to one and return.
        if s == "ta":
            self._tscale_fac = np.full(self._nint, 1.0)
            self._tscale = s
            self._set_all_meta("TSCALFAC", 1.0)
            self._set_all_meta("TSCALE", self.tscale)
            self._set_all_meta("BUNIT", ntscale)
            self._set_all_meta("TUNIT7", ntscale)
            return

        gc = GBTGainCorrection()
        elev = np.array([x["ELEVATIO"] for x in self._meta]) * u.degree
        freq = np.array([x["CRVAL1"] for x in self._meta]) * u.Hz
        date = Time([x["DATE"] for x in self._meta], scale="utc", format="isot")
        factor = gc.scale_ta_to(
            tscale, freq, elev, date, zenith_opacity, zd=False, surface_error=self._surface_error, ap_eff=self._ap_eff
        )
        self._scaleby(factor)
        self._tscale_fac = factor
        self._tscale = s
        self._set_all_meta("BUNIT", ntscale)
        self._set_all_meta("TUNIT7", ntscale)
        self._set_all_meta("TSCALE", self._tscale)
        for i in range(len(self._meta)):
            self._meta[i]["TSCALFAC"] = self.tscale_fac[i]

    def _set_all_meta(self, key, value):
        for i in range(len(self._meta)):
            self._meta[i][key] = value

    def _get_all_meta(self, key):
        return [x[key] for x in self._meta]

    def _check_model(self, model, c0, sa, tol):
        # make sure flux units match
        if model.return_units != c0.unit:
            raise ValueError(f"Units of model {model.return_units} and calibrated data {c0.unit} must be the same.")
        # Warn if domain of model doesn't encompass domain of spectral axis.
        # Sort domain and spectral axis to make them both ascending.
        domain = sorted(model.domain) * model.input_units
        ssa = sorted(sa.value) * sa.unit
        cdelt = ssa[1] - ssa[0]
        toldelt = abs(tol * cdelt)
        diff0 = ssa[0] - domain[0]
        diff1 = domain[-1] - ssa[-1]
        if (diff0 < -toldelt) or (diff1 > toldelt):
            raise ValueError(f"Baseline model would extrapolate on spectral axis by more than {tol} channels.")

    @log_call_to_history
    def subtract_baseline(self, model, tol: int = 1, force: bool = False):
        """
        Subtract a (previously computed) baseline model from every integration in this Scan.

        Parameters
        ----------
        model : `~astropy.modeling.Model`
            The baseline model to subtract. This is typically a `~specutils.utils.quantity_model.QuantityModel`
            derived by removing a baseline from a similar spectrum.
        tol : int, optional
            The maximum number of channels on either end of the spectrum to extrapolate the baseline model,
            if the spectral domain of the baseline model is smaller than the spectral axis of the Scan.
            For instance, if `tol=1`, then
            one channel on the low frequency and one channel on the high frequency end are allowed to be extrapolated.
            The default is 1.
        force : bool, optional
            Force subtraction of the input baseline model, even if another baseline model has been previously subtracted.
            Note: The previous baseline model will **not** be undone (added back in) before subtraction of the input baseline model.
            The default is False.

        Raises
        ------
        ValueError
            If the data are not yet calibrated or the tolerance `tol` is exceeded.

        Returns
        -------
        None

        """
        if self._calibrated is None:
            raise ValueError("Data must be calibrated before a baseline can be subtracted.")
        if self._subtracted:
            if not force:
                warnings.warn(
                    "A baseline model has already been subtracted from this scan. Use 'force=True' to force removal of another model.",
                    stacklevel=2,
                )
                return
        if tol < 0:
            raise ValueError("tol must be non-negative.")
        c0 = self.getspec(0)
        sa = c0.spectral_axis
        self._check_model(model, c0, sa, tol)
        self._calibrated -= model(sa).value
        self._subtracted = True
        self._baseline_model = model

    def undo_baseline(self):
        """
        Undo the applied (subtracted) baseline. The subtracted baseline
        will be added back to the data. The `baseline_model` attribute is set to None.
        """
        if self._baseline_model is None:
            return
        sa = self.getspec(0).spectral_axis
        self._calibrated += self._baseline_model(sa).value
        self._baseline_model = None
        self._subtracted = False

    @property
    def baseline_model(self):
        """Returns the subtracted baseline model or None if it has not yet been computed."""
        return self._baseline_model

    @property
    def calibrated(self) -> np.ndarray:
        """Returns the calibrated integrations in the Scan as a numpy array."""
        return self._calibrated

    @property
    def subtracted(self) -> bool:
        """Has a baseline model been subtracted?

        Returns
        -------
        True if a baseline model has been subtracted, False otherwise
        """
        return self._subtracted

    @property
    def scan(self) -> int:
        """
        The scan number

        Returns
        -------
        int
            The scan number of the integrations in the Scan object
        """
        return self._scan

    @property
    def ap_eff(self) -> np.ndarray:
        """
        The aperture efficiencies for the integrations in this Scan

        Returns
        -------
        ap_eff : `~numpy.ndarray`
            The aperture efficiencies, an array of floats between 0 and 1, one value per integration.
        """
        # compute it the first time if not computed.
        if self._ap_eff_array is None:
            if self._ap_eff is not None:
                self._ap_eff_array = np.full(self.nint, fill_value=self._ap_eff, dtype=float)
            else:
                gc = GBTGainCorrection()
                elev = np.array([x["ELEVATIO"] for x in self._meta]) * u.degree
                freq = np.array([x["CRVAL1"] for x in self._meta]) * u.Hz
                date = Time([x["DATE"] for x in self._meta], scale="utc", format="isot")
                self._ap_eff_array = gc.aperture_efficiency(
                    specval=freq, angle=elev, date=date, zd=False, surface_error=self.surface_error
                )

        return self._ap_eff_array

    @property
    def surface_error(self) -> u.quantity.Quantity:
        """
        The dish surface errors for the integrations in this Scan

        Returns
        -------
        surface_error : `~astropy.units.quantity.Quantity`
            The surface_errors with dimension length, one value per integration.
        """
        # compute it the first time if not computed.
        if self._surface_error_array is None:
            if self._surface_error is not None:
                self._surface_error_array = (
                    np.full(self.nint, fill_value=self._surface_error.to(u.micron).value, dtype=float) * u.micron
                )
            else:
                gc = GBTGainCorrection()
                date = Time([x["DATE"] for x in self._meta], scale="utc", format="isot")
                # While any Quantity with dimensions length is fine; be nice and convert to
                # traditional micron units
                self._surface_error_array = gc._surface_error_array(date).to(u.micron)
        return self._surface_error_array

    @property
    def zenith_opacity(self) -> float | None:
        """
        The zenith opacity of this Scan, used to calculated aperture efficiency.

        Returns
        -------
        zenith_opacity: float or None
            A non-negative float. If None, no zenith opacity was set for this scan (e.g. TPScan); the data cannot be scaled.
        """
        return self._zenith_opacity

    @property
    def nchan(self) -> int:
        """
        The number of channels in this scan.

        Returns
        -------
        int
            The number of channels in this scan.

        """
        return self._nchan

    @property
    def nint(self) -> int:
        """
        The number of integrations in this scan.

        Returns
        -------
        int
            The number of integrations in this scan.
        """
        return self._nint

    @property
    def nrows(self) -> int:
        """
        The number of rows in this scan.

        Returns
        -------
        int
            The number of rows in this scan.
        """
        return self._nrows

    @property
    def ifnum(self) -> int:
        """
        The intermediate frequency (IF) number.

        Returns
        -------
        int
            The index of the IF.
        """
        return self._ifnum

    @property
    def fdnum(self) -> int:
        """
        The feed number.

        Returns
        -------
        int
            The index of the feed.
        """
        return self._fdnum

    @property
    def plnum(self) -> int:
        """
        The polarization number.

        Returns
        -------
        int
            The polarization number.
        """
        return self._plnum

    # @property
    # @todo implement this using Evans crval4_to_pol
    # def pols(self):
    #     """The polarization description

    #     Returns
    #    -------
    #      list
    #         The list of integer polarization number(s)
    #     """
    #    return self._pols

    @property
    def is_calibrated(self) -> bool:
        """
        Have the data been calibrated?

        Returns
        -------
        bool
            True if the data have been calibrated, False if not.

        """
        return self._calibrated is not None

    @property
    def meta(self) -> dict:
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
        # remember, self._meta is a list of dicts,
        # and despite its name it is not Table metadata, it
        # it Table Columns values
        d = {}
        if len(self.history) != 0:
            d["HISTORY"] = self.history
        if len(self.comments) != 0:
            d["COMMENT"] = self.comments
        return Table(self._meta, meta=d)

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
        df = self._sdfits.index(bintable=self._bintable_index).iloc[rowindices]
        self._meta = df.to_dict("records")  # returns dict(s) with key = row number.
        for i in range(len(self._meta)):
            if "CUNIT1" not in self._meta[i]:
                self._meta[i]["CUNIT1"] = (
                    "Hz"  # @todo this is in gbtfits.hdu[0].header['TUNIT11'] but is it always TUNIT11?
                )
            self._meta[i]["CUNIT2"] = "deg"  # is this always true?
            self._meta[i]["CUNIT3"] = "deg"  # is this always true?
            restfrq = self._meta[i]["RESTFREQ"]
            rfq = restfrq * u.Unit(self._meta[i]["CUNIT1"])
            restfreq = rfq.to("Hz").value
            self._meta[i]["RESTFRQ"] = restfreq  # WCS wants no E
            self._meta[i]["BUNIT"] = self._tscale_to_unit[self.tscale.lower()].to_string()
            self._meta[i]["TUNIT7"] = self._meta[i]["BUNIT"]
            self._meta[i]["TSCALE"] = self.tscale
            self._meta[i]["CRPIX1"] -= self._channel_slice.start  # adjustment for user trimmed channels

    def _add_calibration_meta(self):
        """Add metadata that are computed after calibration."""
        if not self.is_calibrated:
            raise Exception("Data have to be calibrated first to add calibration metadata")
        for i in range(len(self._meta)):
            self._meta[i]["TSYS"] = self._tsys[i]
            self._meta[i]["TCAL"] = self._tcal[i]
            # I'm not really sure we should be setting NAXIS1 as this is a FITS reserve keyword.
            # For not leave it in as some tests depend on it.
            # @todo ask PJT
            self._meta[i]["NAXIS1"] = len(self._calibrated[i])
            self._meta[i]["TSYS"] = self._tsys[i]
            self._meta[i]["EXPOSURE"] = self.exposure[i]
            self._meta[i]["DURATION"] = self.duration[i]

    def _update_scale_meta(self):
        """Update metadata that described how integrations were scaled to Ta, Ta*, or Flux"""
        a = self.ap_eff
        s = self.surface_error.value
        seunit = str(self.surface_error.unit)
        for i in range(len(self._meta)):
            self._meta[i]["BUNIT"] = self._tscale_to_unit[self.tscale.lower()].to_string()
            self._meta[i]["TSCALE"] = self.tscale
            self._meta[i]["TSCALFAC"] = self.tscale_fac[i]
            self._meta[i]["AP_EFF"] = a[i]
            self._meta[i]["SURF_ERR"] = s[i]
            self._meta[i]["SE_UNIT"] = seunit

    @abstractmethod
    def calibrate(self, **kwargs):  ## SCANBASE
        """Calibrate the Scan data"""
        pass

    def get_vane_tcal(self):
        """Get the tcal value for the vane."""
        dateobs = self._get_all_meta("DATE-OBS")
        mjd = isot_to_mjd(dateobs)
        obsfreq = np.mean(self._get_all_meta("OBSFREQ"))
        elevation = self._get_all_meta("ELEVATIO")
        tcal = self._vane._get_tcal(
            obsfreq * u.Hz,
            mjd,
            elevation * u.deg,
            zenith_opacity=self._zenith_opacity,
        )
        return tcal

    @log_call_to_history
    @copy_docstring(SpectralAverageMixin.timeaverage)
    def timeaverage(self, weights="tsys", use_wcs=True):  ## SCANBASE
        if self._calibrated is None or len(self._calibrated) == 0:
            raise Exception("You can't time average before calibration.")
        self._timeaveraged = deepcopy(self.getspec(0, use_wcs=use_wcs))
        data = self._calibrated
        w = None
        if isinstance(weights, np.ndarray):
            ws = weights.shape
            if ws != len(self) and ws != self.tsys_weight.shape and ws != (len(self), self.nchan):
                raise ValueError(
                    f"Bad shape for weight array: {ws}. Was expecting {len(self)}, {self.tsys_weight.shape}, or ({len(self)}, {self.nchan})."
                )
            if w is None:
                w = weights
        elif weights == "tsys":
            w = self.tsys_weight
        elif weights is None:
            w = np.ones_like(self.tsys_weight)
        else:
            raise ValueError("Unrecognized weights: must be 'tsys', None, or an array of numbers")
        data_avg, sum_of_weights = np.ma.average(data, axis=0, weights=w, returned=True)
        self._timeaveraged._data = data_avg
        self._timeaveraged.mask = data_avg.mask
        self._timeaveraged._data.set_fill_value(np.nan)
        non_blanks = find_non_blanks(data)
        if w.shape == (len(self), self.nchan):
            w_collapsed = np.average(w, axis=1)
        else:
            w_collapsed = w
        self._timeaveraged.meta["MEANTSYS"] = np.mean(self._tsys[non_blanks])
        self._timeaveraged.meta["WTTSYS"] = sq_weighted_avg(
            self._tsys[non_blanks], axis=0, weights=w_collapsed[non_blanks]
        )
        self._timeaveraged.meta["EXPOSURE"] = np.sum(self._exposure[non_blanks])
        self._timeaveraged.meta["DURATION"] = np.sum(self._duration[non_blanks])
        self._timeaveraged.meta["TSYS"] = self._timeaveraged.meta["WTTSYS"]
        self._timeaveraged.meta["AP_EFF"] = sq_weighted_avg(
            self.ap_eff[non_blanks], axis=0, weights=w_collapsed[non_blanks]
        )
        self._timeaveraged.meta["SURF_ERR"] = sq_weighted_avg(
            self.surface_error[non_blanks].value, axis=0, weights=w_collapsed[non_blanks]
        )
        if self.zenith_opacity is not None:
            self._timeaveraged.meta["TAU_Z"] = self.zenith_opacity
        self._timeaveraged._weights = sum_of_weights
        self._weights = w
        return self._timeaveraged

    def _make_bintable(self, flags: bool) -> BinTableHDU:
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
        # Table metadata aren't preserved in BinTableHDU, so we
        # have to grab them here and add them
        # data_table = self._meta_as_table()
        # table_meta = data_table.meta
        cd = BinTableHDU(data=self._meta_as_table(), name="SINGLE DISH").columns
        form = f"{np.shape(self._calibrated)[1]}E"
        cd.add_col(Column(name="DATA", format=form, array=self._calibrated))
        logger.debug(f"Writing {len(self._calibrated)} rows for output from ScanBase.")
        # re-arrange so DATA is column 7
        cd1 = cd[:6] + cd[-1] + cd[6:-1]
        if flags:
            flags = self._calibrated.mask.astype(np.uint8)
            flagform = f"{np.shape(flags)[1]}B"
            cd1.add_col(Column(name="FLAGS", format=flagform, array=flags))
        b = BinTableHDU.from_columns(cd1, name="SINGLE DISH")
        return b

    def write(self, fileobj, flags=True, output_verify="exception", overwrite=False, checksum=False):
        """
        Write an SDFITS format file (FITS binary table HDU) of the calibrated data in this ScanBase

        Parameters
        ----------
        fileobj : str, file-like or `pathlib.Path`
            File to write to.  If a file object, must be opened in a
            writeable mode.
        flags: bool, optional
            If True, write the applied flags to a `FLAGS` column in the binary table.
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
        self._make_bintable(flags=flags).writeto(
            name=fileobj, output_verify=output_verify, overwrite=overwrite, checksum=checksum
        )

    def plot(self, **kwargs):
        if self._plotter is None:
            self._plotter = sp.ScanPlot(self, **kwargs)
        self._plotter.plot(**kwargs)
        return self._plotter

    def __len__(self):
        return self._nint

    def _init_tsys(self, tsys=None):
        """
        Initialize the array of system temperature values.
        This assumes that the input system temperature is either
        a scalar or a (N,) array with N equal to the number of integrations.
        """

        # Noise diode firing and no user provided tsys.
        if not self._nocal and tsys is None:
            self._tsys = np.full(self._nint, np.nan, dtype=float)
        # User provided tsys or TSYS column.
        elif tsys is not None:
            self._tsys = np.ones(self._nint, dtype=float) * tsys[: self._nint]

    def _init_tcal(self, tcal=None):
        """
        Initialize the array of noise diode temperatures.
        """

        self._tcal = np.empty((self._nint), dtype=float)
        self._tcal[:] = tcal


class ScanBlock(UserList, HistoricalBase, SpectralAverageMixin):
    @log_call_to_history
    def __init__(self, *args):
        UserList.__init__(self, *args)
        HistoricalBase.__init__(self)
        self._nrows = 0
        self._npol = 0  # always 1?
        self._nfeed = 0  # always 1?
        self._nif = 0  # always 1?
        self._timeaveraged = []
        self._plotter = None

    def _scanblock_property(self, prop: str, desc: str):
        """
        Utility method to return a property which should have only one value

        Parameters
        ----------
        prop : str
            The property name
        desc : str
            A descriptive string for the warning message if needed

        Raises
        ------
        AttributeError
            If the ScanBlock doesn't have the property

        Returns
        -------
        value : Any
            The property value.
        """
        if not hasattr(self.data[0], prop):
            raise AttributeError("'ScanBlock' object has no attribute '{prop}'")
        _prop = set([getattr(scan, prop) for scan in self.data])
        if len(_prop) > 1:
            logger.warning(f"The Scans in this ScanBlock have differing {desc} {_prop}")
            return list(_prop)
        return list(_prop)[0]  # noqa: RUF015

    def _aggregate_scan_property(self, prop: str) -> np.ndarray:
        """
        Utility method to collect values of a particular property of all
        Scans in this ScanBlock

        Parameters
        ----------
        prop : str
            The property name.

        Raises
        ------
        AttributeError
            If a Scan of the ScanBlock does not have the property.

        Returns
        -------
        value : `~numpy.ndarray`
            The values of the Scan property
        """
        # check the first scan
        if not hasattr(self.data[0], prop):
            raise AttributeError(f"Contained scans have no attribute '{prop}'")
        # There is no guarantee that each scan in the scanblock has the same number
        # of integrations, so we can't create a conventional numpy array of [nscan,nintegration]
        value = np.empty(len(self.data), dtype=np.ndarray)
        for i, scan in enumerate(self.data):
            # Will raise AttributeError if scan does not have prop
            value[i] = getattr(scan, prop)
            i += 1
        return value

    @property
    def tsys(self):
        """
        The system temperatures for all scans in this ScanBlock

        Returns
        -------
        tsys :  `~numpy.ndarray`
            The system temperatures for all scans in this ScanBlock
        """
        return self._aggregate_scan_property("tsys")

    @property
    def delta_freq(self):
        """
        Returns
        -------
        delta_freq : `~numpy.ndarray`
            The channel frequency spacings for all scans in this ScanBlock
        """
        return self._aggregate_scan_property("delta_freq")

    @property
    def exposure(self):
        """
        Returns
        -------
        exposure : `~numpy.ndarray`
            The exposure times for all scans in this ScanBlock
        """
        return self._aggregate_scan_property("exposure")

    @property
    def duration(self):
        """
        Returns
        -------
        duration : `~numpy.ndarray`
            The duration times for all scans in this ScanBlock
        """
        return self._aggregate_scan_property("duration")

    @property
    def nchan(self):
        """
        The number of channels in the first Scan in this ScanBlock.
        (Assumes number of channels in each enclosed Scan are the same).

        Returns
        -------
        int
            The number of channels in this ScanBlock.

        """
        return self.data[0].nchan

    @property
    def nint(self):
        """The total number of integrations in this Scanblock

        Returns
        -------
        int
            The number of integerationsin this ScanBlock.

        """
        return np.sum([i.nint for i in self])

    @property
    def weights(self) -> list:
        """
        The weights associated with the Scans and integrations in this ScanBlock

        Returns
        -------
        list
            A list containing the weight arrays for each scan. Because Scans can have different numbers of integrations,
            this is a list instead of a numpy array.
        """
        weights = []
        for scan in self.data:
            weights.append(scan.weights)
        return weights

    @log_call_to_history
    def calibrate(self, **kwargs):
        """Calibrate all scans in this ScanBlock"""
        for scan in self.data:
            scan.calibrate(**kwargs)

    @log_call_to_history
    def smooth(self, method="hanning", width=1, decimate=0, kernel=None):  # ScanBlock
        """
        Smooth all scans in this ScanBlock.
        Smooth or convolve the  calibrated data arrays in contained Scans, optionally decimating the data.
        A number of methods from astropy.convolution can be selected
        with the `method` keyword.

        Default smoothing is hanning.

        Note: Any previously computed/removed baseline will remain unchanged.

        Parameters
        ----------.
        method : string, optional
            Smoothing method. Valid are: 'hanning', 'boxcar' and
            'gaussian'. Minimum match applies.
            The default is 'hanning'.
        width : int, optional
            Effective width of the convolving kernel.  Should ideally be an
            odd number.
            For 'hanning' this should be 1, with a 0.25,0.5,0.25 kernel.
            For 'boxcar' an even value triggers an odd one with half the
            signal at the edges, and will thus not reproduce GBTIDL.
            For 'gaussian' this is the FWHM of the final beam. We normally
            assume the input beam has FWHM=1, pending resolution on cases
            where CDELT1 is not the same as FREQRES.
            The default is 1.
        decimate : int, optional
            Decimation factor of the spectrum by returning every decimate channel.
            -1:   no decimation
            0:    use the width parameter
            >1:   user supplied decimation (use with caution)
            The default is 0, meaning decimation is by `width`
        kernel : `~numpy.ndarray`, optional
            A numpy array which is the kernel by which the signal is convolved.
            Use with caution, as it is assumed the kernel is normalized to
            one, and is symmetric. Since width is ill-defined here, the user
            should supply an appropriate number manually.
            NOTE: not implemented yet.
            The default is None.
        Raises
        ------
        Exception
            If no valid smoothing method is given.

        Returns
        -------
        None

        """
        for scan in self.data:
            scan.smooth(method, width, decimate, kernel)

    @log_call_to_history
    @copy_docstring(SpectralAverageMixin.timeaverage)
    def timeaverage(self, weights="tsys", use_wcs: bool = True):  ## SCANBLOCK
        # average of the averages
        ary = False
        if isinstance(weights, np.ndarray):
            ws = weights.shape
            ary = True
            if ws != (self.nint,) and ws != (self.nint, self.nchan):
                raise ValueError(f"Bad shape for weight array: {ws}. Was expecting {len(self)},)  or ({self.nint},).")
        self._timeaveraged = []
        index = 0
        if ary:
            for scan in self.data:
                w = weights[index : index + scan.nint]
                index += scan.nint
                self._timeaveraged.append(scan.timeaverage(w, use_wcs=use_wcs))
            s = average_spectra(self._timeaveraged, weights="spectral")
        else:
            for scan in self.data:
                self._timeaveraged.append(scan.timeaverage(weights, use_wcs=use_wcs))
            s = average_spectra(self._timeaveraged, weights=weights)
        s.merge_commentary(self)
        return s

    @log_call_to_history
    def scale(self, tscale, zenith_opacity):
        """
        Scale all the data in this `ScanBlock` to the given brightness scale and zenith opacity. If data are already
        scaled, they will be unscaled first.

        Parameters
        ----------
        tscale : str
            Strings representing valid options for scaling spectral data, specifically
                - 'Ta'  : Antenna Temperature
                - 'Ta*' : Antenna temperature corrected to above the atmosphere
                - 'Flux'  : flux density in Jansky
            This parameter is case-insensitive.

        zenith_opacity : float
            The zenith opacity

        Returns
        -------
        None

        Raises
        ------
        TypeError
            If scaling to temperature is not applicable to the scan type, e.g., a total power scan.
        ValueError
            if `tscale` is unrecognized or `zenith_opacity` is negative.

        """
        for scan in self.data:
            scan.scale(tscale, zenith_opacity)

    @property
    def tscale(self):
        """
        The descriptive brightness unit of the data.  One of
            - 'Raw' : raw value, e.g., count
            - 'Ta'  : Antenna Temperature
            - 'Ta*' : Antenna temperature corrected to above the atmosphere
            - 'Flux'  : flux density in Jansky

        Returns
        -------
        tscale : str
            brightness unit string
        """
        return self._scanblock_property("tscale", "brightness scales")

    @property
    def tscale_fac(self):
        """
        The factor(s) by which the data have been scaled from antenna temperature to corrected antenna temperature
        or flux density.

        Returns
        -------
        tscale_fac : `~numpy.ndarray`
            An array of floats, one per integration in the ScanBlock.

        """
        return np.squeeze(np.array([scan.tscale_fac for scan in self.data]))

    @property
    def tunit(self):
        """The brightness unit (temperature or flux density) of this ScanBlocks's data

        Returns
        -------
        tunit : `~astropy.units.Unit`
            The brightness unit
        """
        return self._scanblock_property("tunit", "brightness units")

    # possible @todo:  We could have a baseline() method with same signature as Spectrum.baseline, which would compute
    # timeaverage for each Scan in a ScanBlock, and for each Scan calculate and remove that baseline from t
    # the integrations in that Scan.

    @log_call_to_history
    def subtract_baseline(self, model, tol: int = 1, force: bool = False):
        """
        Subtract a (previously computed) baseline model from every integration of every Scan in this ScanBlock.

        Parameters
        ----------
        model : `~astropy.modeling.Model`
            The baseline model to subtract. This is typically a `~specutils.utils.quantity_model.QuantityModel`
            derived by removing a baseline from a similar spectrum.
        tol : int, optional
            The maximum number of channels on either end of the spectrum to extrapolate the baseline model,
            if the spectral domain of the baseline model is smaller than the spectral axis of the Scan. For instance, if `tol=1`, then
            one channel on the low frequency and one channel on the high frequency end are allowed to be extrapolated.
            The default is 1.
        force : bool, optional
            Force subtraction of the input baseline model, even if another baseline model has been previously subtracted.
            Note: The previous baseline model will **not** be undone (added back in) before subtraction of the input baseline model.
            The default is False.

        Raises
        ------
        ValueError
            If the data are not yet calibrated or the tolerance `tol` is exceeded.

        Returns
        -------
        None
        """
        for scan in self.data:
            scan.subtract_baseline(model, tol, force)

    @log_call_to_history
    def undo_baseline(self):
        """
        For all Scans in this ScanBlock, undo the applied (subtracted) baseline. The subtracted baseline
        will be added back to the data. The Scan's baseline_model` attribute is set to None.
        """
        for scan in self.data:
            scan.undo_baseline()

    def write(self, fileobj, flags=True, output_verify="exception", overwrite=False, checksum=False):
        """
        Write an SDFITS format file (FITS binary table HDU) of the calibrated data in this ScanBlock

        Parameters
        ----------
        fileobj : str, file-like or `pathlib.Path`
            File to write to.  If a file object, must be opened in a
            writeable mode.
        flags: bool, optional
            If True, write the applied flags to a `FLAGS` column in the binary table.
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
        # Meta are the keys of the first scan's bintable except for DATA.
        # We can use this to compare with the subsequent Scan
        # keywords without have to create their bintables first.
        defaultkeys = set(s0._meta[0].keys())
        datashape = np.shape(s0._calibrated)
        tablelist = [s0._meta_as_table()]
        nrows = 0
        for scan in self.data:  # [1:]:
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
            if nrows > 0:
                tablelist.append(scan._meta_as_table())
            nrows = nrows + thisshape[0]
        # now do the same trick as in Scan.write() of adding "DATA" to the coldefs
        # astropy Tables can be concatenated with vstack thankfully.
        table = vstack(tablelist, join_type="exact")
        # need to preserve table.meta because it gets lost in created of "cd" ColDefs
        table_meta = table.meta
        cd = BinTableHDU(table, name="SINGLE DISH").columns
        data = np.concatenate([c._calibrated.filled(np.nan) for c in self.data])
        form = f"{np.shape(data)[1]}E"
        cd.add_col(Column(name="DATA", format=form, array=data))
        # re-arrange so DATA is column 7
        cd1 = cd[:6] + cd[-1] + cd[6:-1]
        if flags:
            flags = np.concatenate([c._calibrated.mask * 1 for c in self.data]).astype(np.uint8)
            flagform = f"{np.shape(flags)[1]}B"
            cd1.add_col(Column(name="FLAGS", format=flagform, array=flags))
        b = BinTableHDU.from_columns(cd1, name="SINGLE DISH")

        # preserve any meta
        for k, v in table_meta.items():
            if k == "HISTORY" or k == "COMMENT":
                continue  # will deal with these later
            b.header[k] = v
        for h in self._history:
            b.header["HISTORY"] = h
        for c in self._comments:
            b.header["COMMENT"] = c

        logger.debug(f"Saving {nrows} scans in {fileobj} from ScanBlock.")
        b.writeto(name=fileobj, output_verify=output_verify, overwrite=overwrite, checksum=checksum)

    def plot(self, **kwargs):
        if self._plotter is None:
            self._plotter = sp.ScanPlot(self, **kwargs)
        self._plotter.plot(**kwargs)
        return self._plotter


class TPScan(ScanBase):
    """GBT specific version of Total Power Scan

    Parameters
    ----------
    gbtfits : `~dysh.fits.gbtfitsload.GBTFITSLoad`
        input GBTFITSLoad object
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
    fdnum: int
        The feed number
    ifnum : int
        The IF number
    plnum : int
        The polarization number
    bintable : int
        the index for BINTABLE in `sdfits` containing the scans
    calibrate: bool
        whether or not to calibrate the data.  If `True`, the data will be (calon + caloff)*0.5, otherwise it will be SDFITS row data. Default:True
    smoothref: int
        the number of channels in the reference to boxcar smooth prior to calibration
    apply_flags : boolean, optional.  If True, apply flags before calibration.
    observer_location : `~astropy.coordinates.EarthLocation`
        Location of the observatory. See `~dysh.coordinates.Observatory`.
        This will be transformed to `~astropy.coordinates.ITRS` using the time of
        observation DATE-OBS or MJD-OBS in
        the SDFITS header.  The default is the location of the GBT.
    tscale : str, optional
        The brightess unit scale for the output scan, must be one of (case-insensitive)
            - 'Raw' : raw value, e.g., count
            - 'Ta'  : Antenna Temperature
            - 'Ta*' : Antenna temperature corrected to above the atmosphere
            - 'Flux'  : flux density in Jansky
        Default: 'Raw'
    channel: list or None
        An inclusive list of `[firstchan, lastchan]` to use in the calibration. The channel list is zero-based. If provided,
        only data channels in the inclusive range `[firstchan,lastchan]` will be used. If a reference spectrum has been given, it will also be
        trimmed to `[firstchan,lastchan]`. If channels have already been selected through
        :meth:`~dysh.fits.gbtfitsload.GBTFITSLoad.select_channel`, a ValueError will be raised.
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
        fdnum,
        ifnum,
        plnum,
        bintable,
        calibrate=True,
        smoothref=1,
        apply_flags=False,
        tsys=None,
        tcal=None,
        observer_location=Observatory["GBT"],
        tscale="Raw",  # @todo why is this even an exposed parameter for TPScan?
        channel=None,
    ):
        ScanBase.__init__(
            self, gbtfits, smoothref, apply_flags, observer_location, fdnum, ifnum, plnum, tcal=tcal, channel=channel
        )
        self._sdfits = gbtfits  # parent class
        self._scan = scan
        self._sigstate = sigstate
        self._calstate = calstate
        self._scanrows = scanrows
        self._smoothref = smoothref
        self._apply_flags = apply_flags
        self._observer_location = observer_location
        self._tscale = tscale
        self._tunit = self._tscale_to_unit[self._tscale.lower()]
        if self._smoothref > 1:
            raise NotImplementedError(f"TP smoothref={self._smoothref} not implemented yet")

        # @todo deal with data that crosses bintables
        if bintable is None:
            self._bintable_index = self._sdfits._find_bintable_and_row(self._scanrows[0])[0]
        else:
            self._bintable_index = bintable
        df = self._sdfits._index
        df = df.iloc[scanrows]
        self._index = df
        self._nint = 0
        self._timeaveraged = None
        self._nrows = len(scanrows)
        self._tsys = None
        self._calrows = calrows
        # all cal=T states where sig=sigstate
        self._refonrows = sorted(list(set(self._calrows["ON"]).intersection(set(self._scanrows))))
        # all cal=F states where sig=sigstate
        self._refoffrows = sorted(list(set(self._calrows["OFF"]).intersection(set(self._scanrows))))
        self._refcalon = gbtfits.rawspectra(self._bintable_index, setmask=apply_flags)[
            self._refonrows, self._channel_slice
        ]
        self._refcaloff = gbtfits.rawspectra(self._bintable_index, setmask=apply_flags)[
            self._refoffrows, self._channel_slice
        ]
        nb1 = find_non_blanks(self._refcalon)
        nb2 = find_non_blanks(self._refcaloff)
        goodrows = np.intersect1d(nb1, nb2)
        if len(self._refcalon) == 0:
            # special case for notpcal (when calrows["ON"] is 0)
            self._calstate = False  # set calstate to false so that _calc_exposure() doesn't raise exception
            goodrows = np.intersect1d(nb2, nb2)  # isn't this just nb2.flatten()?
            if len(goodrows) != len(self._refcaloff):
                nblanks = len(self._refcaloff) - len(goodrows)
                logger.info(f"Ignoring {nblanks} blanked integration(s).")
            self._refcalon = None
            self._refcaloff = self._refcaloff[goodrows]
            self._refonrows = []
            self._refoffrows = [
                self._refoffrows[i] for i in goodrows
            ]  # why not self._refoffrows[goodrows] ?? -> because it is a list.
            self._nchan = len(self._refcaloff[0])
        else:
            # Tell the user about blank integration(s) that will be ignored.
            if len(goodrows) != len(self._refcalon):
                nblanks = len(self._refcalon) - len(goodrows)
                logger.info(f"Ignoring {nblanks} blanked integration(s).")
            self._refcalon = self._refcalon[goodrows]
            self._refcaloff = self._refcaloff[goodrows]
            self._refonrows = [self._refonrows[i] for i in goodrows]
            self._refoffrows = [self._refoffrows[i] for i in goodrows]
            self._nrows = len(self._refonrows) + len(self._refoffrows)  # ??
            self._nchan = len(self._refcalon[0])
        self._calc_exposure()
        self._calc_delta_freq()
        self._validate_defaults()
        # Use 'ta' as tscale in this call so that scaling is not attempted.
        self._finish_initialization(calibrate, None, self._refoffrows, "ta", None, tsys=tsys, tcal=tcal)

    def calibrate(self, **kwargs):  ## TPSCAN
        """Calibrate the total power data according to the CAL/SIG table above"""
        # the way the data are formed depend only on cal state
        # since we have downselected based on sig state in the constructor
        if self.calstate is None:
            self._calibrated = (0.5 * (self._refcalon + self._refcaloff)).astype(float)
        elif self.calstate:
            self._calibrated = self._refcalon.astype(float)
        elif self.calstate == False:  # noqa: E712
            self._calibrated = self._refcaloff.astype(float)
        else:
            raise Exception(f"Unrecognized cal state {self.calstate}")  # should never happen
        if np.all(np.isnan(self._tsys)):
            self._calc_tsys()

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

    def _calc_tsys(self, **kwargs):
        """
        Calculate the system temperature array, according to table above.
        """
        # self._tcal = list(self._sdfits.index(bintable=self._bintable_index).iloc[self._refoffrows]["TCAL"])
        if len(self._calrows["ON"]) == 0:
            if np.all(np.isnan(self._tsys)):
                self._tsys = np.ones(self._nint, dtype=float)
        else:
            nspect = len(self._tcal)
            self._tsys = np.empty(nspect, dtype=float)  # should be same as len(calon)
            if len(self._tcal) != nspect:
                raise Exception(f"TCAL length {len(self._tcal)} and number of spectra {nspect} don't match")
            for i in range(nspect):
                tsys = mean_tsys(calon=self._refcalon[i], caloff=self._refcaloff[i], tcal=self._tcal[i])
                self._tsys[i] = tsys

    def _calc_exposure(self):
        """Calculate the exposure time for TPScan.

        The value depends on the cal state:

           =====  ======================================
           CAL    EXPOSURE
           =====  ======================================
           None   :math:`t_{EXP,REFCALON} + t_{EXP,REFCALOFF}`
           True   :math:`t_{EXP,REFCALON}`
           False  :math:`t_{EXP,REFCALOFF}`
           =====  ======================================

        where `REFCALON` = integrations with `cal=T` and  `REFCALOFF` = integrations with `cal=F`.
        """
        if self.calstate is None:
            exp_ref_on = self._sdfits.index(bintable=self._bintable_index).iloc[self._refonrows]["EXPOSURE"].to_numpy()
            exp_ref_off = (
                self._sdfits.index(bintable=self._bintable_index).iloc[self._refoffrows]["EXPOSURE"].to_numpy()
            )
            dur_ref_on = self._sdfits.index(bintable=self._bintable_index).iloc[self._refonrows]["DURATION"].to_numpy()
            dur_ref_off = (
                self._sdfits.index(bintable=self._bintable_index).iloc[self._refoffrows]["DURATION"].to_numpy()
            )

        elif self.calstate:
            exp_ref_on = self._sdfits.index(bintable=self._bintable_index).iloc[self._refonrows]["EXPOSURE"].to_numpy()
            exp_ref_off = 0
            dur_ref_on = self._sdfits.index(bintable=self._bintable_index).iloc[self._refonrows]["DURATION"].to_numpy()
            dur_ref_off = 0
        elif self.calstate == False:  # noqa: E712
            exp_ref_on = 0
            exp_ref_off = (
                self._sdfits.index(bintable=self._bintable_index).iloc[self._refoffrows]["EXPOSURE"].to_numpy()
            )
            dur_ref_on = 0
            dur_ref_off = (
                self._sdfits.index(bintable=self._bintable_index).iloc[self._refoffrows]["DURATION"].to_numpy()
            )

        self._exposure = exp_ref_on + exp_ref_off
        self._duration = dur_ref_on + dur_ref_off

    def _calc_delta_freq(self):  # TPSCAN
        """Calculate the channel width.

        The value depends on the cal state:

        =====  ================================================================
        CAL     :math:`\\Delta\nu`
        =====  ================================================================
        None    :math:`0.5 * ( \\Delta\nu_{REFON}+ \\Delta\nu_{REFOFF} )`
        True    :math:`\\Delta\nu_{REFON}`
        False   :math:`\\Delta\nu_{REFOFF}`
        =====  ================================================================
        """
        df_ref_on = self._sdfits.index(bintable=self._bintable_index).iloc[self._refonrows]["CDELT1"].to_numpy()
        df_ref_off = self._sdfits.index(bintable=self._bintable_index).iloc[self._refoffrows]["CDELT1"].to_numpy()
        if self.calstate is None:
            delta_freq = 0.5 * (df_ref_on + df_ref_off)
        elif self.calstate:
            delta_freq = df_ref_on
        elif self.calstate == False:  # noqa: E712
            delta_freq = df_ref_off
        self._delta_freq = delta_freq

    def total_power(self, i):
        """Return the i-th total power spectrum in this Scan.
        This is a synonym for :meth:`calibrated`

        Parameters
        ----------
        i : int
            The index into the data array

        Returns
        -------
        spectrum : `~dysh.spectra.spectrum.Spectrum`
        """
        return self.getspec(i)


class PSScan(ScanBase):
    """GBT specific version of Position Switch Scan. A position switch scan object has
    one IF, one feed, and one polarization

    Parameters
    ----------
    gbtfits : `~dysh.fits.gbtfitsload.GBTFITSLoad`
        Input GBTFITSLoad object.
    scan : dict
        dictionary with keys 'ON' and 'OFF' containing unique list of ON (signal) and OFF (reference) scan numbers NOTE: there should be one ON and one OFF, a pair.
    scanrows : dict
        dictionary with keys 'ON' and 'OFF' containing the list of rows in `sdfits` corresponding to ON (signal) and OFF (reference) integrations.
    calrows : dict
        dictionary containing with keys 'ON' and 'OFF' containing list of rows in `sdfits` corresponding to cal=T (ON) and cal=F (OFF) integrations.
    fdnum: int
        The feed number.
    ifnum : int
        The intermediate frequency (IF) number.
    plnum : int
        The polarization number.
    bintable : int
        The index for BINTABLE in `sdfits` containing the scans.
    calibrate: bool
        Whether or not to calibrate the data. If true, data will be calibrated as TSYS*(ON-OFF)/OFF. Default: True
    smoothref: int
        If >1 smooth the reference with a boxcar kernel with a width of `smooth_ref` channels. The default is to not smooth the reference.
    apply_flags : boolean, optional
        If True, apply flags before calibration.
    observer_location : `~astropy.coordinates.EarthLocation`
        Location of the observatory. See `~Observatory`.
        This will be transformed to `~astropy.coordinates.ITRS` using the time of
        observation DATE-OBS or MJD-OBS in
        the SDFITS header.  The default is the location of the GBT.
    tscale : str, optional
        The brightess unit scale for the output scan, must be one of (case-insensitive)
                - 'Ta'  : Antenna Temperature
                - 'Ta*' : Antenna temperature corrected to above the atmosphere
                - 'Flux'  : flux density in Jansky
        If 'ta*' or 'flux' the zenith opacity must also be given. Default: 'ta'
    zenith_opacity: float, optional
        The zenith opacity to use in calculating the scale factors for the integrations. Default: None
    refspec : int or `~dysh.spectra.spectrum.Spectrum`, optional
        If given, the Spectrum will be used as the reference rather than using scan data.
    tsys : float or `~numpy.ndarray`
        If given, this is the system temperature in Kelvin. It overrides the values calculated using the noise diodes.
        If not given, and signal and reference are scan numbers, the system temperature will be calculated from the reference
        scan and the noise diode. If not given, and the reference is a `~dysh.spectra.spectrum.Spectrum`, the reference system temperature as given
        in the metadata header will be used. The default is to use the noise diode or the metadata, as appropriate.
        If `vane` is provided, `tsys` will be ignored.
    ap_eff : float or None
        Aperture efficiency to be used when scaling data to brightness temperature of flux. The provided aperture
        efficiency must be a number between 0 and 1.  If None, `dysh` will calculate it as described in
        :meth:`~dysh.util.GBTGainCorrection.aperture_efficiency`. Only one of `ap_eff` or `surface_error`
        can be provided.
    surface_error: Quantity or None
        Surface rms error, in units of length (typically microns), to be used in the Ruze formula when calculating the
        aperture efficiency.  If None, `dysh` will use the known GBT surface error model.  Only one of `ap_eff` or `surface_error`
        can be provided.
    channel: list or None
        An inclusive list of `[firstchan, lastchan]` to use in the calibration. The channel list is zero-based. If provided,
        only data channels in the inclusive range `[firstchan,lastchan]` will be used. If a reference spectrum has been given, it will also be
        trimmed to `[firstchan,lastchan]`. If channels have already been selected through
        :meth:`~dysh.fits.gbtfitsload.GBTFITSLoad.select_channel`, a ValueError will be raised.
    vane : `~dysh.spectra.vane.VaneSpectrum` or None
        Vane calibration spectrum. This will be used to derive the system temperature.
        If provided, `tsys` will be ignored.
    """

    def __init__(
        self,
        gbtfits,
        scan,
        scanrows,
        calrows,
        fdnum,
        ifnum,
        plnum,
        bintable,
        calibrate=True,
        smoothref=1,
        apply_flags=False,
        observer_location=Observatory["GBT"],
        tscale="ta",
        zenith_opacity=0.0,
        refspec=None,
        tsys=None,
        ap_eff=None,
        surface_error=None,
        tcal=None,
        nocal=False,
        channel=None,
        vane: VaneSpectrum | None = None,
    ):
        ScanBase.__init__(
            self,
            gbtfits,
            smoothref,
            apply_flags,
            observer_location,
            fdnum,
            ifnum,
            plnum,
            tsys,
            ap_eff=ap_eff,
            surface_error=surface_error,
            zenith_opacity=zenith_opacity,
            tcal=tcal,
            channel=channel,
        )
        # The rows of the original bintable corresponding to ON (sig) and OFF (reg)
        # self._history = deepcopy(gbtfits._history)
        self._scan = scan["ON"]
        self._sigscan = scan["ON"]
        self._refscan = scan["OFF"]
        self._scanrows = scanrows
        self._nrows = len(self._scanrows["ON"])
        self._refspec = refspec
        if isinstance(self.refspec, Spectrum):
            self._has_refspec = True
            self._refspec = self._refspec[self._channel_slice]
        else:
            self._has_refspec = False
        self._sigspec = None
        self._nocal = nocal
        self._vane = vane

        # calrows perhaps not needed as input since we can get it from gbtfits object?
        # calrows['ON'] are rows with noise diode was on, regardless of sig or ref
        # calrows['OFF'] are rows with noise diode was off, regardless of sig or ref
        self._calrows = calrows
        # @todo deal with data that crosses bintables
        if bintable is None:
            self._bintable_index = gbtfits._find_bintable_and_row(self._scanrows["ON"][0])[0]
        else:
            self._bintable_index = bintable
        # noise diode on, signal position
        self._sigonrows = sorted(list(set(self._calrows["ON"]).intersection(set(self._scanrows["ON"]))))
        # noise diode off, signal position
        self._sigoffrows = sorted(list(set(self._calrows["OFF"]).intersection(set(self._scanrows["ON"]))))
        self._sigcalon = gbtfits.rawspectra(self._bintable_index, setmask=apply_flags)[
            self._sigonrows, self._channel_slice
        ]
        self._sigcaloff = gbtfits.rawspectra(self._bintable_index, setmask=apply_flags)[
            self._sigoffrows, self._channel_slice
        ]

        if self._has_refspec:
            self._refoffrows = None
            self._refoffrows = None
            self._refcalon = None
            self._refcaloff = None
            # Catch blank integrations.
            goodrows = find_nonblank_ints(self._sigcaloff, self._sigcalon)
        else:
            # noise diode on, reference position
            self._refonrows = sorted(list(set(self._calrows["ON"]).intersection(set(self._scanrows["OFF"]))))
            # noise diode off, reference position
            self._refoffrows = sorted(list(set(self._calrows["OFF"]).intersection(set(self._scanrows["OFF"]))))
            self._refcalon = gbtfits.rawspectra(self._bintable_index, setmask=apply_flags)[
                self._refonrows, self._channel_slice
            ]
            self._refcaloff = gbtfits.rawspectra(self._bintable_index, setmask=apply_flags)[
                self._refoffrows, self._channel_slice
            ]

            # Catch blank integrations.
            if not self._nocal:
                goodrows = find_nonblank_ints(self._sigcaloff, self._refcaloff, self._sigcalon, self._refcalon)
                self._refcalon = self._refcalon[goodrows]
                self._refcaloff = self._refcaloff[goodrows]
                self._refonrows = [self._refonrows[i] for i in goodrows]
                self._refoffrows = [self._refoffrows[i] for i in goodrows]
                self._sigcalon = self._sigcalon[goodrows]
                self._sigcaloff = self._sigcaloff[goodrows]
                self._sigonrows = [self._sigonrows[i] for i in goodrows]
                self._sigoffrows = [self._sigoffrows[i] for i in goodrows]
                # Update number of rows after removing blanks.
                nsigrows = len(self._sigonrows) + len(self._sigoffrows)
                self._nrows = nsigrows
            else:
                goodrows = find_nonblank_ints(self._sigcaloff, self._refcaloff)
                self._refcaloff = self._refcaloff[goodrows]
                self._refoffrows = [self._refoffrows[i] for i in goodrows]
                self._sigcaloff = self._sigcaloff[goodrows]
                self._sigoffrows = [self._sigoffrows[i] for i in goodrows]
                # Update number of rows after removing blanks.
                nsigrows = len(self._sigoffrows)
                self._nrows = nsigrows
        nchan = gbtfits.nchan(self._bintable_index)
        if self._channel_slice == slice(0, None) or self._channel_slice == slice(0, nchan):
            self._nchan = nchan
        else:
            self._nchan = len(self._sigcalon[0])
        self._finish_initialization(
            calibrate,
            None,
            meta_rows=self._sigoffrows,
            tscale=tscale,
            zenith_opacity=zenith_opacity,
            tsys=tsys,
            tcal=tcal,
        )

    @property
    def sigscan(self) -> int:
        """The scan number associated with the signal"""
        return self._sigscan

    @property
    def refscan(self) -> int | None:
        """The scan number associated with the reference.

        Returns
        -------
         int or None;
             Integer scan number or None if the reference was a Spectrum object
        """
        return self._refscan

    @property
    def sigspec(self) -> Spectrum | None:
        """The signal Spectrum if one was given at construction.

        Returns
        -------
         Spectrum or None;
             Spectrum object if given as signal or None.
        """
        return self._sigspec

    @property
    def refspec(self) -> Spectrum | None:
        """The reference Spectrum if one was given at construction.

        Returns
        -------
         Spectrum or None;
             Spectrum object if given as reference or None.
        """
        return self._refspec

    def calibrate(self, **kwargs):  ##PSSCAN
        """
        Position switch calibration, following equations 1 and 2 in the GBTIDL calibration manual
        """
        kwargs_opts = {"verbose": False}
        kwargs_opts.update(kwargs)
        if self._smoothref > 1 and kwargs_opts["verbose"]:
            logger.debug(f"PSScan smoothref={self._smoothref}")
        nspect = self._nint
        self._calibrated = np.ma.empty((nspect, self._nchan), dtype="d")

        if self._vane is not None:
            self._tcal[:] = self.get_vane_tcal()

        if self._has_refspec:
            if self._smoothref > 1:
                ref, _meta = core.smooth(self.refspec.data, "boxcar", self._smoothref)
            else:
                ref = self.refspec.data
            for i in range(nspect):
                if not self._nocal:
                    sig = 0.5 * (self._sigcalon[i] + self._sigcaloff[i])
                else:
                    sig = self._sigcaloff[i]
                if self._vane is not None:
                    self._tsys[i] = self._vane._get_tsys(ref, self._tcal[i])
                tsys = self._tsys[i]
                self._calibrated[i] = tsys * (sig - ref) / ref
        else:
            tcal = self._tcal
            if len(tcal) != nspect:
                raise Exception(f"TCAL length {len(tcal)} and number of spectra {nspect} don't match")
            if not self._nocal:
                for i in range(nspect):
                    if not np.isnan(self._tsys[i]):
                        tsys = self._tsys[i]
                    else:
                        tsys = mean_tsys(calon=self._refcalon[i], caloff=self._refcaloff[i], tcal=tcal[i])
                    sig = 0.5 * (self._sigcalon[i] + self._sigcaloff[i])
                    ref = 0.5 * (self._refcalon[i] + self._refcaloff[i])
                    if self._smoothref > 1:
                        ref, _meta = core.smooth(ref, "boxcar", self._smoothref)
                    self._calibrated[i] = tsys * (sig - ref) / ref
                    self._tsys[i] = tsys
            else:
                for i in range(nspect):
                    sig = self._sigcaloff[i]
                    ref = self._refcaloff[i]
                    if self._smoothref > 1:
                        ref, _meta = core.smooth(ref, "boxcar", self._smoothref)
                    if self._vane is not None:
                        self._tsys[i] = self._vane._get_tsys(ref, self._tcal[i])
                    tsys = self._tsys[i]
                    self._calibrated[i] = tsys * (sig - ref) / ref
        logger.debug(f"Calibrated {nspect} PSScan spectra")

    def _calc_exposure(self):
        """The array of exposure (integration) times for PSScan

        exposure = exp_sig * exp_ref * nsmooth/(exp_sig+exp_ref *nsmooth)

        with `nsmooth` the reference spectrum smoothing parameter, `exp_sig` the exposure of the signal
        spectra, and `exp_ref` the exposure of the reference spectra. If the signal and reference spectra were
        observed with two noise diode states, then their exposure times are the sum of the exposure times in each state.

        Returns
        -------
        exposure : ~numpy.ndarray
            The exposure time in units of the EXPOSURE keyword in the SDFITS header
        """
        exp_sig_on = self._sdfits.index(bintable=self._bintable_index).iloc[self._sigonrows]["EXPOSURE"].to_numpy()
        exp_sig_off = self._sdfits.index(bintable=self._bintable_index).iloc[self._sigoffrows]["EXPOSURE"].to_numpy()
        dur_sig_on = self._sdfits.index(bintable=self._bintable_index).iloc[self._sigonrows]["DURATION"].to_numpy()
        dur_sig_off = self._sdfits.index(bintable=self._bintable_index).iloc[self._sigoffrows]["DURATION"].to_numpy()
        if self._has_refspec:
            exp_ref = self.refspec.meta.get("EXPOSURE", None)
            dur_ref = self.refspec.meta.get("DURATION", None)
            if exp_ref is None:
                raise ValueError(
                    "Can't set exposure time for PSScan integrations because reference spectrum has no exposure time in its metadata. Solve with refspec.meta['EXPOSURE']=value."
                )
            if dur_ref is None:
                raise ValueError(
                    "Can't set duration time for PSScan integrations because reference spectrum has no duration time in its metadata. Solve with refspec.meta['DURATION']=value."
                )
        else:
            exp_ref_on = self._sdfits.index(bintable=self._bintable_index).iloc[self._refonrows]["EXPOSURE"].to_numpy()
            exp_ref_off = (
                self._sdfits.index(bintable=self._bintable_index).iloc[self._refoffrows]["EXPOSURE"].to_numpy()
            )
            dur_ref_on = self._sdfits.index(bintable=self._bintable_index).iloc[self._refonrows]["DURATION"].to_numpy()
            dur_ref_off = (
                self._sdfits.index(bintable=self._bintable_index).iloc[self._refoffrows]["DURATION"].to_numpy()
            )

            if not self._nocal:
                exp_ref = exp_ref_on + exp_ref_off
                dur_ref = dur_ref_on + dur_ref_off
            else:
                exp_ref = exp_ref_off
                dur_ref = dur_ref_off
        if not self._nocal:
            exp_sig = exp_sig_on + exp_sig_off
            dur_sig = dur_sig_on + dur_sig_off
        else:
            exp_sig = exp_sig_off
            dur_sig = dur_sig_off
        if self._smoothref > 1:
            nsmooth = self._smoothref
        else:
            nsmooth = 1.0
        self._exposure = exp_sig * exp_ref * nsmooth / (exp_sig + exp_ref * nsmooth)
        self._duration = dur_sig * dur_ref * nsmooth / (dur_sig + dur_ref * nsmooth)

    def _calc_delta_freq(self):
        """calculate the channel width

        If the calibration diode has been fired

        df = [ 0.5*(df_ref_on + df_ref_off) + 0.5*(df_sig_on + df_sig_off) ] / 2

        otherwise

        df = 0.5 * (df_ref_off + df_sig_off)

        where `df_ref_on` and `df_ref_off` are the channel widths of the reference spectra for cal='F' and cal='T',
        respectively and `df_sig_on` and `df_sig_off` are the channel widths of the signal spectra for cal='F' and cal='T',
        respectively.

        """
        df_sig_on = self._sdfits.index(bintable=self._bintable_index).iloc[self._sigonrows]["CDELT1"].to_numpy()
        df_sig_off = self._sdfits.index(bintable=self._bintable_index).iloc[self._sigoffrows]["CDELT1"].to_numpy()
        if self._has_refspec:
            df_ref_on = df_ref_off = np.full_like(self._sigoffrows, self.refspec.meta["CDELT1"])
        else:
            df_ref_on = self._sdfits.index(bintable=self._bintable_index).iloc[self._refonrows]["CDELT1"].to_numpy()
            df_ref_off = self._sdfits.index(bintable=self._bintable_index).iloc[self._refoffrows]["CDELT1"].to_numpy()
        if not self._nocal:
            df_ref = 0.5 * (df_ref_on + df_ref_off)
            df_sig = 0.5 * (df_sig_on + df_sig_off)
        else:
            df_ref = df_ref_off
            df_sig = df_sig_off
        self._delta_freq = 0.5 * (df_ref + df_sig)


class NodScan(ScanBase):
    """GBT specific version of Nodding Scan. A nod scan object has
    one IF, two feeds, and one polarization.

    Parameters
    ----------
    gbtfits : `~dysh.fits.gbtfitsload.GBTFITSLoad`
        input GBTFITSLoad object
    scan : dict
        dictionary with keys 'ON' and 'OFF' containing unique list of ON (signal) and OFF (reference) scan numbers
        NOTE: there should be one ON and one OFF, a pair. There should be at least two beams (the nodding beams)
        which will be resp. on source in each scan.
    beam1: bool
        Is this scan BEAM1 or BEAM2?  BEAM1 is defined as being on source in the first scan of the pair, BEAM2 in the second of the pair
    scanrows : dict
        dictionary with keys 'ON' and 'OFF' containing the list of rows in `sdfits` corresponding to ON (signal) and OFF (reference) integrations
    calrows : dict
        dictionary containing with keys 'ON' and 'OFF' containing list of rows in `sdfits` corresponding to cal=T (ON) and cal=F (OFF) integrations.
    fdnum: int
        The feed number
    ifnum : int
        The IF number
    plnum : int
        The polarization number
    bintable : int
        The index for BINTABLE in `sdfits` containing the scans
    calibrate: bool
        Whether or not to calibrate the data.  If true, data will be calibrated as TSYS*(ON-OFF)/OFF.
        Default: True
    smoothref: int
        The number of channels in the reference to boxcar smooth prior to calibration (if applicable)
    apply_flags : boolean
        If True, apply flags before calibration.
    tsys : float
        User provided value for the system temperature.
        If `vane` is provided, `tsys` will be ignored.
    nocal : bool
        True if the noise diode was not fired. False if it was fired.
    tscale : str, optional
        The brightness scale unit for the output scan, must be one of (case-insensitive)
                - 'Ta'  : Antenna Temperature
                - 'Ta*' : Antenna temperature corrected to above the atmosphere
                - 'Flux'  : flux density in Jansky
        If 'ta*' or 'flux' the zenith opacity must also be given. Default:'ta'
    ap_eff : float or None
        Aperture efficiency to be used when scaling data to brightness temperature of flux. The provided aperture
        efficiency must be a number between 0 and 1.  If None, `dysh` will calculate it as described in
        :meth:`~dysh.util.GBTGainCorrection.aperture_efficiency`. Only one of `ap_eff` or `surface_error`
        can be provided.
    surface_error: Quantity or None
        Surface rms error, in units of length (typically microns), to be used in the Ruze formula when calculating the
        aperture efficiency.  If None, `dysh` will use the known GBT surface error model.  Only one of `ap_eff` or `surface_error`
        can be provided.
    zenith_opacity: float, optional
        The zenith opacity to use in calculating the scale factors for the integrations.  Default:None
    observer_location : `~astropy.coordinates.EarthLocation`
        Location of the observatory. See `~dysh.coordinates.Observatory`.
        This will be transformed to `~astropy.coordinates.ITRS` using the time of
        observation DATE-OBS or MJD-OBS in
        the SDFITS header.  The default is the location of the GBT.
    channel: list or None
        An inclusive list of `[firstchan, lastchan]` to use in the calibration. The channel list is zero-based. If provided,
        only data channels in the inclusive range `[firstchan,lastchan]` will be used. If a reference spectrum has been given, it will also be
        trimmed to `[firstchan,lastchan]`. If channels have already been selected through
        :meth:`~dysh.fits.gbtfitsload.GBTFITSLoad.select_channel`, a ValueError will be raised.
    vane : `~dysh.spectra.vane.VaneSpectrum` or None
        Vane calibration spectrum. This will be used to derive the system temperature.
        If provided, `tsys` will be ignored.
    """

    def __init__(
        self,
        gbtfits,
        scan,
        beam1,
        scanrows,
        calrows,
        fdnum,
        ifnum,
        plnum,
        bintable,
        calibrate=True,
        smoothref=1,
        apply_flags=False,
        tsys=None,
        tcal=None,
        nocal=False,
        tscale="ta",
        ap_eff=None,
        surface_error=None,
        zenith_opacity=None,
        observer_location=Observatory["GBT"],
        channel=None,
        vane: VaneSpectrum | None = None,
    ):
        ScanBase.__init__(
            self,
            gbtfits,
            smoothref,
            apply_flags,
            observer_location,
            fdnum,
            ifnum,
            plnum,
            tsys,
            ap_eff=ap_eff,
            surface_error=surface_error,
            zenith_opacity=zenith_opacity,
            tcal=tcal,
            channel=channel,
        )
        self._scan = scan["ON"]
        self._scanrows = scanrows
        self._nrows = len(self._scanrows["ON"])
        self._beam1 = beam1
        self._nocal = nocal
        self._vane = vane

        # calrows perhaps not needed as input since we can get it from gbtfits object?
        # calrows['ON'] are rows with noise diode was on, regardless of sig or ref
        # calrows['OFF'] are rows with noise diode was off, regardless of sig or ref
        self._calrows = calrows
        # @todo deal with data that crosses bintables
        if bintable is None:
            self._bintable_index = gbtfits._find_bintable_and_row(self._scanrows["ON"][0])[0]
        else:
            self._bintable_index = bintable
        if False:
            self._nint = gbtfits.nintegrations(self._bintable_index)
        # so quick with slicing!
        self._sigonrows = sorted(list(set(self._calrows["ON"]).intersection(set(self._scanrows["ON"]))))
        self._sigoffrows = sorted(list(set(self._calrows["OFF"]).intersection(set(self._scanrows["ON"]))))
        self._refonrows = sorted(list(set(self._calrows["ON"]).intersection(set(self._scanrows["OFF"]))))
        self._refoffrows = sorted(list(set(self._calrows["OFF"]).intersection(set(self._scanrows["OFF"]))))
        if beam1:
            self._sigcalon = gbtfits.rawspectra(self._bintable_index, setmask=apply_flags)[
                self._sigonrows, self._channel_slice
            ]
            self._sigcaloff = gbtfits.rawspectra(self._bintable_index, setmask=apply_flags)[
                self._sigoffrows, self._channel_slice
            ]
            self._refcalon = gbtfits.rawspectra(self._bintable_index, setmask=apply_flags)[
                self._refonrows, self._channel_slice
            ]
            self._refcaloff = gbtfits.rawspectra(self._bintable_index, setmask=apply_flags)[
                self._refoffrows, self._channel_slice
            ]
        else:
            self._sigcalon = gbtfits.rawspectra(self._bintable_index, setmask=apply_flags)[
                self._refonrows, self._channel_slice
            ]
            self._sigcaloff = gbtfits.rawspectra(self._bintable_index, setmask=apply_flags)[
                self._refoffrows, self._channel_slice
            ]
            self._refcalon = gbtfits.rawspectra(self._bintable_index, setmask=apply_flags)[
                self._sigonrows, self._channel_slice
            ]
            self._refcaloff = gbtfits.rawspectra(self._bintable_index, setmask=apply_flags)[
                self._sigoffrows, self._channel_slice
            ]

        # Catch blank integrations.

        if not self._nocal:
            goodrows = find_nonblank_ints(self._sigcaloff, self._refcaloff, self._sigcalon, self._refcalon)
            self._refcalon = self._refcalon[goodrows]
            self._refcaloff = self._refcaloff[goodrows]
            self._refonrows = [self._refonrows[i] for i in goodrows]
            self._refoffrows = [self._refoffrows[i] for i in goodrows]
            self._sigcalon = self._sigcalon[goodrows]
            self._sigcaloff = self._sigcaloff[goodrows]
            self._sigonrows = [self._sigonrows[i] for i in goodrows]
            self._sigoffrows = [self._sigoffrows[i] for i in goodrows]
            # Update number of rows after removing blanks.
            nsigrows = len(self._sigonrows) + len(self._sigoffrows)
            self._nrows = nsigrows
        else:
            goodrows = find_nonblank_ints(self._sigcaloff, self._refcaloff)
            self._refcaloff = self._refcaloff[goodrows]
            self._refoffrows = [self._refoffrows[i] for i in goodrows]
            self._sigcaloff = self._sigcaloff[goodrows]
            self._sigoffrows = [self._sigoffrows[i] for i in goodrows]
            # Update number of rows after removing blanks.
            nsigrows = len(self._sigoffrows)
            self._nrows = nsigrows

        self._nchan = len(self._sigcaloff[0])
        self._finish_initialization(calibrate, None, self._sigoffrows, tscale, zenith_opacity, tsys=tsys, tcal=tcal)

    def calibrate(self, **kwargs):  ##NODSCAN
        """
        NodScan calibration
        """
        kwargs_opts = {"verbose": False}
        kwargs_opts.update(kwargs)
        if self._smoothref > 1 and kwargs_opts["verbose"]:
            logger.debug(f"NodScan smoothref={self._smoothref}")
        nspect = self._nint
        self._calibrated = np.ma.empty((nspect, self._nchan), dtype="d")
        self._calc_exposure()
        if self._vane is not None:
            self._tcal[:] = self.get_vane_tcal()
        tcal = self._tcal
        if len(tcal) != nspect:
            raise Exception(f"TCAL length {len(tcal)} and number of spectra {nspect} don't match")
        if not self._nocal:
            for i in range(nspect):
                if not np.isnan(self._tsys[i]):
                    tsys = self._tsys[i]
                else:
                    tsys = mean_tsys(calon=self._refcalon[i], caloff=self._refcaloff[i], tcal=tcal[i])
                sig = 0.5 * (self._sigcalon[i] + self._sigcaloff[i])
                ref = 0.5 * (self._refcalon[i] + self._refcaloff[i])
                if self._smoothref > 1:
                    ref, _meta = core.smooth(ref, "boxcar", self._smoothref)
                self._calibrated[i] = tsys * (sig - ref) / ref
                self._tsys[i] = tsys
        else:
            for i in range(nspect):
                sig = self._sigcaloff[i]
                ref = self._refcaloff[i]
                if self._smoothref > 1:
                    ref, _meta = core.smooth(ref, "boxcar", self._smoothref)
                if self._vane is not None:
                    self._tsys[i] = self._vane._get_tsys(ref, self._tcal[i])
                tsys = self._tsys[i]
                self._calibrated[i] = tsys * (sig - ref) / ref
        logger.debug(f"Calibrated {nspect} NODScan spectra")

    def _calc_exposure(self):
        """The array of exposure (integration) times for NodScan

        exposure = exp_sig * exp_ref * nsmooth / (exp_sig + exp_ref * nsmooth)

        with `nsmooth` the reference spectrum smoothing parameter, `exp_sig` the exposure of the signal
        spectra, and `exp_ref` the exposure of the reference spectra. If the signal and reference spectra were
        observed with two noise diode states, then their exposure times are the sum of the exposure times in each state.

        Returns
        -------
        exposure : ~numpy.ndarray
            The exposure time in units of the EXPOSURE keyword in the SDFITS header
        """
        exp_ref_on = self._sdfits.index(bintable=self._bintable_index).iloc[self._refonrows]["EXPOSURE"].to_numpy()
        exp_ref_off = self._sdfits.index(bintable=self._bintable_index).iloc[self._refoffrows]["EXPOSURE"].to_numpy()
        exp_sig_on = self._sdfits.index(bintable=self._bintable_index).iloc[self._sigonrows]["EXPOSURE"].to_numpy()
        exp_sig_off = self._sdfits.index(bintable=self._bintable_index).iloc[self._sigoffrows]["EXPOSURE"].to_numpy()
        dur_ref_on = self._sdfits.index(bintable=self._bintable_index).iloc[self._refonrows]["DURATION"].to_numpy()
        dur_ref_off = self._sdfits.index(bintable=self._bintable_index).iloc[self._refoffrows]["DURATION"].to_numpy()
        dur_sig_on = self._sdfits.index(bintable=self._bintable_index).iloc[self._sigonrows]["DURATION"].to_numpy()
        dur_sig_off = self._sdfits.index(bintable=self._bintable_index).iloc[self._sigoffrows]["DURATION"].to_numpy()
        if not self._nocal:
            exp_ref = exp_ref_on + exp_ref_off
            exp_sig = exp_sig_on + exp_sig_off
            dur_ref = dur_ref_on + dur_ref_off
            dur_sig = dur_sig_on + dur_sig_off
        else:
            exp_ref = exp_ref_off
            exp_sig = exp_sig_off
            dur_ref = dur_ref_off
            dur_sig = dur_sig_off
        if self._smoothref > 1:
            nsmooth = self._smoothref
        else:
            nsmooth = 1.0
        self._exposure = exp_sig * exp_ref * nsmooth / (exp_sig + exp_ref * nsmooth)
        self._duration = dur_sig * dur_ref * nsmooth / (dur_sig + dur_ref * nsmooth)

    def _calc_delta_freq(self):
        """Get the array of channel frequency width

        If the calibration diode has been fired

        df = [ 0.5*(df_ref_on + df_ref_off) + 0.5*(df_sig_on + df_sig_off) ] / 2

        otherwise

        df = 0.5 * (df_ref_off + df_sig_off)

        where `df_ref_on` and `df_ref_off` are the channel widths of the reference spectra for cal='F' and cal='T',
        respectively and `df_sig_on` and `df_sig_off` are the channel widths of the signal spectra for cal='F' and cal='T',
        respectively.

        Returns
        -------
             delta_freq: ~numpy.ndarray
                 The channel frequency width in units of the CDELT1 keyword in the SDFITS header
        """
        df_ref_on = self._sdfits.index(bintable=self._bintable_index).iloc[self._refonrows]["CDELT1"].to_numpy()
        df_ref_off = self._sdfits.index(bintable=self._bintable_index).iloc[self._refoffrows]["CDELT1"].to_numpy()
        df_sig_on = self._sdfits.index(bintable=self._bintable_index).iloc[self._sigonrows]["CDELT1"].to_numpy()
        df_sig_off = self._sdfits.index(bintable=self._bintable_index).iloc[self._sigoffrows]["CDELT1"].to_numpy()
        if not self._nocal:
            df_ref = 0.5 * (df_ref_on + df_ref_off)
            df_sig = 0.5 * (df_sig_on + df_sig_off)
        else:
            df_ref = df_ref_off
            df_sig = df_sig_off
        self._delta_freq = 0.5 * (df_ref + df_sig)


class FSScan(ScanBase):
    """GBT specific version of Frequency Switch Scan

    Parameters
    ----------
    gbtfits : `~dysh.fits.gbtfitsload.GBTFITSLoad`
        input GBTFITSLoad object
    scan : int
        Scan number that contains integrations with a series of sig/ref and calon/caloff states.
    sigrows :dict
        Dictionary containing with keys 'ON' and 'OFF' containing list of rows in `sdfits`
        corresponding to sig=T (ON) and sig=F (OFF) integrations.
    calrows : dict
        Dictionary containing with keys 'ON' and 'OFF' containing list of rows in `sdfits`
        corresponding to cal=T (ON) and cal=F (OFF) integrations.
    fdnum: int
        The feed number
    ifnum : int
        The IF number
    plnum : int
        The polarization number
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
    apply_flags : boolean, optional.  If True, apply flags before calibration.
    tscale : str, optional
        The brightness scale unit for the output scan, must be one of (case-insensitive)
                - 'Ta'  : Antenna Temperature
                - 'Ta*' : Antenna temperature corrected to above the atmosphere
                - 'Flux'  : flux density in Jansky
        If 'ta*' or 'flux' the zenith opacity must also be given. Default:'ta'
    ap_eff : float or None
        Aperture efficiency to be used when scaling data to brightness temperature of flux. The provided aperture
        efficiency must be a number between 0 and 1.  If None, `dysh` will calculate it as described in
        :meth:`~dysh.util.GBTGainCorrection.aperture_efficiency`. Only one of `ap_eff` or `surface_error`
        can be provided.
    surface_error: Quantity or None
        Surface rms error, in units of length (typically microns), to be used in the Ruze formula when calculating the
        aperture efficiency.  If None, `dysh` will use the known GBT surface error model.  Only one of `ap_eff` or `surface_error`
        can be provided.
    zenith_opacity: float, optional
        The zenith opacity to use in calculating the scale factors for the integrations.  Default:None
    observer_location : `~astropy.coordinates.EarthLocation`
        Location of the observatory. See `~dysh.coordinates.Observatory`.
        This will be transformed to `~astropy.coordinates.ITRS` using the time of
        observation DATE-OBS or MJD-OBS in the SDFITS header.  The default is the location of the GBT.
    tsys : float or `~numpy.ndarray`
        If given, this is the system temperature in Kelvin. It overrides the values calculated using the noise diodes.
        If not given, and signal and reference are scan numbers, the system temperature will be calculated from the reference
        scan and the noise diode. If not given, and the reference is a `~dysh.spectra.spectrum.Spectrum`, the reference system temperature as given
        in the metadata header will be used. The default is to use the noise diode or the metadata, as appropriate.
    channel: list or None
        An inclusive list of `[firstchan, lastchan]` to use in the calibration. The channel list is zero-based. If provided,
        only data channels in the inclusive range `[firstchan,lastchan]` will be used. If a reference spectrum has been given, it will also be
        trimmed to `[firstchan,lastchan]`. If channels have already been selected through
        :meth:`~dysh.fits.gbtfitsload.GBTFITSLoad.select_channel`, a ValueError will be raised.
        If `vane` is provided, `tsys` will be ignored.
    vane : `~dysh.spectra.vane.VaneSpectrum` or None
        Vane calibration spectrum. This will be used to derive the system temperature.
        If provided, `tsys` will be ignored.
    """

    def __init__(
        self,
        gbtfits,
        scan,
        sigrows,
        calrows,
        fdnum,
        ifnum,
        plnum,
        bintable,
        calibrate=True,
        fold=True,
        shift_method="fft",
        use_sig=True,
        smoothref=1,
        apply_flags=False,
        tscale="ta",
        ap_eff=None,
        surface_error=None,
        zenith_opacity=None,
        observer_location=Observatory["GBT"],
        tsys=None,
        tcal=None,
        nocal: bool = False,
        channel: list | None = None,
        vane: VaneSpectrum | None = None,
    ):
        ScanBase.__init__(
            self,
            gbtfits,
            smoothref,
            apply_flags,
            observer_location,
            fdnum,
            ifnum,
            plnum,
            tsys,
            ap_eff=ap_eff,
            surface_error=surface_error,
            zenith_opacity=zenith_opacity,
            tcal=tcal,
            channel=channel,
        )
        # The rows of the original bintable corresponding to ON (sig) and OFF (reg)
        self._scan = scan  # for FS everything is an "ON"
        self._sigrows = sigrows  # dict with "ON" and "OFF"
        self._calrows = calrows  # dict with "ON" and "OFF"
        self._folded = False
        self._use_sig = use_sig
        self._nocal = nocal
        self._vane = vane
        self._smoothref = smoothref
        self._sigonrows = sorted(list(set(self._calrows["ON"]).intersection(set(self._sigrows["ON"]))))
        self._sigoffrows = sorted(list(set(self._calrows["OFF"]).intersection(set(self._sigrows["ON"]))))
        self._refonrows = sorted(list(set(self._calrows["ON"]).intersection(set(self._sigrows["OFF"]))))
        self._refoffrows = sorted(list(set(self._calrows["OFF"]).intersection(set(self._sigrows["OFF"]))))

        logger.debug("---------------------------------------------------")
        logger.debug("FSSCAN: ")
        logger.debug(f"SigOff {self._sigoffrows}")
        logger.debug(f"SigOn {self._sigonrows}")
        logger.debug(f"RefOff {self._refoffrows}")
        logger.debug(f"RefOn {self._refonrows}")

        nsigrows = len(self._sigonrows) + len(self._sigoffrows)
        nrefrows = len(self._refonrows) + len(self._refoffrows)
        if nsigrows != nrefrows:
            raise Exception("Number of sig rows does not match ref rows. Dangerous to proceed")
        logger.debug(f"sigonrows {nsigrows}, {self._sigonrows}")
        self._nrows = nsigrows

        # @todo deal with data that crosses bintables
        if bintable is None:
            self._bintable_index = gbtfits._find_bintable_and_row(self._sigoffrows[0])[0]
        else:
            self._bintable_index = bintable
        logger.debug(f"bintable index is {self._bintable_index}")

        self._scanrows = list(set(self._calrows["ON"])) + list(set(self._calrows["OFF"]))
        self._sigcalon = gbtfits.rawspectra(self._bintable_index, setmask=apply_flags)[
            self._sigonrows, self._channel_slice
        ]
        self._sigcaloff = gbtfits.rawspectra(self._bintable_index, setmask=apply_flags)[
            self._sigoffrows, self._channel_slice
        ]
        self._refcalon = gbtfits.rawspectra(self._bintable_index, setmask=apply_flags)[
            self._refonrows, self._channel_slice
        ]
        self._refcaloff = gbtfits.rawspectra(self._bintable_index, setmask=apply_flags)[
            self._refoffrows, self._channel_slice
        ]

        # Catch blank integrations.
        if not self._nocal:
            goodrows = find_nonblank_ints(self._sigcaloff, self._refcaloff, self._sigcalon, self._refcalon)
            self._refcalon = self._refcalon[goodrows]
            self._refcaloff = self._refcaloff[goodrows]
            self._refonrows = [self._refonrows[i] for i in goodrows]
            self._refoffrows = [self._refoffrows[i] for i in goodrows]
            self._sigcalon = self._sigcalon[goodrows]
            self._sigcaloff = self._sigcaloff[goodrows]
            self._sigonrows = [self._sigonrows[i] for i in goodrows]
            self._sigoffrows = [self._sigoffrows[i] for i in goodrows]
            nsigrows = len(self._sigonrows) + len(self._sigoffrows)
        else:
            goodrows = find_nonblank_ints(self._sigcaloff, self._refcaloff)
            self._refcaloff = self._refcaloff[goodrows]
            self._refoffrows = [self._refoffrows[i] for i in goodrows]
            self._sigcaloff = self._sigcaloff[goodrows]
            self._sigoffrows = [self._sigoffrows[i] for i in goodrows]
            nsigrows = len(self._sigoffrows)

        # Update number of rows after removing blanks.
        self._nrows = nsigrows
        nchan = gbtfits.nchan(self._bintable_index)
        if self._channel_slice == slice(0, None) or self._channel_slice == slice(0, nchan):
            self._nchan = nchan
        else:
            self._nchan = len(self._sigcalon[0])

        self._finish_initialization(
            calibrate,
            {"fold": fold, "shift_method": shift_method},
            self._sigoffrows,
            tscale,
            zenith_opacity,
            tsys=tsys,
            tcal=tcal,
        )

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

    def calibrate(self, **kwargs):  # FSSCAN
        """
        Frequency switch calibration.

        Parameters
        ----------

        fold : bool
            Fold the spectrum or not. Required keyword.
        """
        # @todo upgrade fold from kwarg to arg
        logger.debug(f"FOLD={kwargs['fold']}")
        logger.debug(f"METHOD={kwargs['shift_method']}")

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
            crpix1 -= self._channel_slice.start
            freq = channel_to_frequency(crval1, crpix1, cdelt1, vframe, nchan[0], nint, ndim=ndim)

            # Apply units.
            try:
                cunit1 = u.Unit(df["CUNIT1"])
                # if ndim == 2:
                #   cunit1 = cunit[0]  #  @todo undefined cunit[]
            except KeyError:
                cunit1 = u.Hz

            return freq * cunit1

        def do_sig_ref(sig, ref, tsys, smooth=False):
            """
            smooth=True would implement smoothing the reference (or something)
            """
            return (sig - ref) / ref * tsys

        def do_fold(sig, ref, sig_freq, ref_freq, remove_wrap=False, shift_method="fft"):
            """ """
            chan_shift = (ref_freq[0] - sig_freq[0]) / np.diff(sig_freq).mean()
            logger.debug(f"do_fold: sig_freq0={sig_freq[0]}, ref_freq0={ref_freq[0]}, chan_shift={chan_shift}")
            ref_shift = core.data_shift(ref, chan_shift, remove_wrap=remove_wrap, method=shift_method)
            # @todo weights
            avg = (sig + ref_shift) / 2
            return avg

        kwargs_opts = {"verbose": False}
        kwargs_opts.update(kwargs)
        _fold = kwargs.get("fold", False)
        nspect = self._nint
        self._calibrated = np.ma.empty((nspect, self._nchan), dtype="d")
        df_sig = self._sdfits.index(bintable=self._bintable_index).iloc[self._sigoffrows]
        df_ref = self._sdfits.index(bintable=self._bintable_index).iloc[self._refoffrows]
        logger.debug(f"df_sig {type(df_sig)} len(df_sig)")
        sig_freq = index_frequency(df_sig)
        ref_freq = index_frequency(df_ref)

        if self._vane is not None:
            self._tcal[:] = self.get_vane_tcal()

        #  tcal is the same for REF and SIG, and the same for all integrations actually.
        tcal = self._tcal
        logger.debug(f"TCAL: {len(tcal)} {tcal[0]}")
        if len(tcal) != nspect:
            raise Exception(f"TCAL length {len(tcal)} and number of spectra {nspect} don't match")

        # @todo   the nspect loop could be replaced with clever numpy?
        if not self._nocal:
            for i in range(nspect):
                if not np.isnan(self._tsys[i]):
                    tsys = self._tsys[i]
                    tsys_ref = tsys
                    tsys_sig = tsys
                else:
                    tsys_sig = mean_tsys(calon=self._sigcalon[i], caloff=self._sigcaloff[i], tcal=tcal[i])
                    tsys_ref = mean_tsys(calon=self._refcalon[i], caloff=self._refcaloff[i], tcal=tcal[i])
                if i == 0:
                    logger.debug(f"Tsys(sig/ref)[0]={tsys_sig} / {tsys_ref}")
                tp_sig = 0.5 * (self._sigcalon[i] + self._sigcaloff[i])
                tp_ref = 0.5 * (self._refcalon[i] + self._refcaloff[i])
                if self._smoothref > 1:
                    if self._use_sig:
                        tp_ref, _meta = core.smooth(tp_ref, "boxcar", self._smoothref)
                    else:
                        tp_sig, _meta = core.smooth(tp_sig, "boxcar", self._smoothref)
                cal_sig = do_sig_ref(tp_sig, tp_ref, tsys_ref)
                cal_ref = do_sig_ref(tp_ref, tp_sig, tsys_sig)
                if _fold:
                    cal_sig_fold = do_fold(
                        cal_sig, cal_ref, sig_freq[i], ref_freq[i], shift_method=kwargs["shift_method"]
                    )
                    cal_ref_fold = do_fold(
                        cal_ref, cal_sig, ref_freq[i], sig_freq[i], shift_method=kwargs["shift_method"]
                    )
                    self._folded = True
                    if self._use_sig:
                        self._calibrated[i] = cal_sig_fold
                        self._tsys[i] = tsys_ref
                    else:
                        self._calibrated[i] = cal_ref_fold
                        self._tsys[i] = tsys_sig

                elif self._use_sig:
                    self._calibrated[i] = cal_sig
                    self._tsys[i] = tsys_ref
                else:
                    self._calibrated[i] = cal_ref
                    self._tsys[i] = tsys_sig
        else:
            for i in range(nspect):
                tp_sig = self._sigcaloff[i]
                tp_ref = self._refcaloff[i]
                if self._smoothref > 1:
                    if self._use_sig:
                        tp_ref, _meta = core.smooth(tp_ref, "boxcar", self._smoothref)
                    else:
                        tp_sig, _meta = core.smooth(tp_sig, "boxcar", self._smoothref)
                if self._vane is not None:
                    if self._use_sig:
                        self._tsys[i] = self._vane._get_tsys(tp_ref, self._tcal[i])
                    else:
                        self._tsys[i] = self._vane._get_tsys(tp_sig, self._tcal[i])
                tsys = self._tsys[i]
                cal_sig = do_sig_ref(tp_sig, tp_ref, tsys)
                cal_ref = do_sig_ref(tp_ref, tp_sig, tsys)
                if _fold:
                    cal_sig_fold = do_fold(
                        cal_sig, cal_ref, sig_freq[i], ref_freq[i], shift_method=kwargs["shift_method"]
                    )
                    cal_ref_fold = do_fold(
                        cal_ref, cal_sig, ref_freq[i], sig_freq[i], shift_method=kwargs["shift_method"]
                    )
                    if self._use_sig:
                        self._calibrated[i] = cal_sig_fold
                    else:
                        self._calibrated[i] = cal_ref_fold
                elif self._use_sig:
                    self._calibrated[i] = cal_sig
                else:
                    self._calibrated[i] = cal_ref

        if _fold:
            self._exposure = 2 * self.exposure
            self._duration = 2 * self.duration
        logger.debug(f"Calibrated {nspect} spectra with fold={_fold} and use_sig={self._use_sig}")

    def _calc_exposure(self):
        """Calculate the array of exposure (integration) times for FSscan

        exposure = exp_sig * exp_ref * nsmooth/(exp_sig+exp_ref *nsmooth)

        with `nsmooth` the reference spectrum smoothing parameter, `exp_sig` the exposure of the signal
        spectra, and `exp_ref` the exposure of the reference spectra. If the signal and reference spectra were
        observed with two noise diode states, then their exposure times are the sum of the exposure times in each state.

        Returns
        -------
        exposure : `~numpy.ndarray`
            The exposure time in units of the EXPOSURE keyword in the SDFITS header
        """
        exp_ref_on = self._sdfits.index(bintable=self._bintable_index).iloc[self._refonrows]["EXPOSURE"].to_numpy()
        exp_ref_off = self._sdfits.index(bintable=self._bintable_index).iloc[self._refoffrows]["EXPOSURE"].to_numpy()
        exp_sig_on = self._sdfits.index(bintable=self._bintable_index).iloc[self._sigonrows]["EXPOSURE"].to_numpy()
        exp_sig_off = self._sdfits.index(bintable=self._bintable_index).iloc[self._sigoffrows]["EXPOSURE"].to_numpy()
        dur_ref_on = self._sdfits.index(bintable=self._bintable_index).iloc[self._refonrows]["DURATION"].to_numpy()
        dur_ref_off = self._sdfits.index(bintable=self._bintable_index).iloc[self._refoffrows]["DURATION"].to_numpy()
        dur_sig_on = self._sdfits.index(bintable=self._bintable_index).iloc[self._sigonrows]["DURATION"].to_numpy()
        dur_sig_off = self._sdfits.index(bintable=self._bintable_index).iloc[self._sigoffrows]["DURATION"].to_numpy()
        if not self._nocal:
            exp_ref = exp_ref_on + exp_ref_off
            exp_sig = exp_sig_on + exp_sig_off
            dur_ref = dur_ref_on + dur_ref_off
            dur_sig = dur_sig_on + dur_sig_off
        else:
            exp_ref = exp_ref_off
            exp_sig = exp_sig_off
            dur_ref = dur_ref_off
            dur_sig = dur_sig_off
        if self._smoothref > 1:
            nsmooth = self._smoothref
        else:
            nsmooth = 1.0
        self._exposure = exp_sig * exp_ref * nsmooth / (exp_sig + exp_ref * nsmooth)
        self._duration = dur_sig * dur_ref * nsmooth / (dur_sig + dur_ref * nsmooth)

    def _calc_delta_freq(self):
        """Get the array of channel frequency width

        df = [ 0.5*(df_ref_on + df_ref_off) + 0.5*(df_sig_on + df_sig_off) ] / 2


        where `df_ref_on` and `df_ref_off` are the channel widths of the reference spectra for cal='F' and cal='T',
        respectively and `df_sig_on` and `df_sig_off` are the channel widths of the signal spectra for cal='F' and cal='T',
        respectively.

        """
        df_ref_on = self._sdfits.index(bintable=self._bintable_index).iloc[self._refonrows]["CDELT1"].to_numpy()
        df_ref_off = self._sdfits.index(bintable=self._bintable_index).iloc[self._refoffrows]["CDELT1"].to_numpy()
        df_sig_on = self._sdfits.index(bintable=self._bintable_index).iloc[self._sigonrows]["CDELT1"].to_numpy()
        df_sig_off = self._sdfits.index(bintable=self._bintable_index).iloc[self._sigoffrows]["CDELT1"].to_numpy()
        if not self._nocal:
            df_ref = 0.5 * (df_ref_on + df_ref_off)
            df_sig = 0.5 * (df_sig_on + df_sig_off)
        else:
            df_ref = df_ref_off
            df_sig = df_sig_off
        self._delta_freq = 0.5 * (df_ref + df_sig)


class SubBeamNodScan(ScanBase):
    r"""
    Parameters
    ----------
    sigtp:  list of `~dysh.spectra.scan.TPScan`
        Signal total power scans
    reftp:  list of `~dysh.spectra.scan.TPScan`
        Reference total power scans
    fdnum: int
        The feed number
    ifnum : int
        The IF number
    plnum : int
        The polarization number
    calibrate: bool
        Whether or not to calibrate the data.
    smoothref: int
        the number of channels in the reference to boxcar smooth prior to calibration
    apply_flags : boolean, optional.  If True, apply flags before calibration.
    tscale : str, optional
        The brightness scale unit for the output scan, must be one of (case-insensitive)
                - 'Ta'  : Antenna Temperature
                - 'Ta*' : Antenna temperature corrected to above the atmosphere
                - 'Flux'  : flux density in Jansky
        If 'ta*' or 'flux' the zenith opacity must also be given. Default:'ta'
    ap_eff : float or None
        Aperture efficiency o be used when scaling data to brightness temperature of flux. The provided aperture
        efficiency must be a number between 0 and 1.  If None, `dysh` will calculate it as described in
        :meth:`~dysh.util.GBTGainCorrection.aperture_efficiency`. Only one of `ap_eff` or `surface_error`
        can be provided.
    surface_error: Quantity or None
        Surface rms error, in units of length (typically microns), to be used in the Ruze formula when calculating the
        aperture efficiency.  If None, `dysh` will use the known GBT surface error model.  Only one of `ap_eff` or `surface_error`
        can be provided.
    zenith_opacity: float, optional
        The zenith opacity to use in calculating the scale factors for the integrations.  Default:None
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
    vane : `~dysh.spectra.vane.VaneSpectrum` or None
        Vane calibration spectrum. This will be used to derive the system temperature.
    """

    def __init__(
        self,
        sigtp,
        reftp,
        fdnum,
        ifnum,
        plnum,
        calibrate=True,
        smoothref=1,
        apply_flags=False,
        tscale="ta",
        ap_eff=None,
        surface_error=None,
        zenith_opacity=None,
        tcal=None,
        observer_location=Observatory["GBT"],
        weights="tsys",
        nocal=False,
        vane: VaneSpectrum | None = None,
        **kwargs,
    ):
        ScanBase.__init__(
            self,
            sigtp[0]._sdfits,
            smoothref,
            apply_flags,
            observer_location,
            fdnum,
            ifnum,
            plnum,
            tcal=tcal,
            ap_eff=ap_eff,
            surface_error=surface_error,
            zenith_opacity=zenith_opacity,
        )
        if len(reftp) != len(sigtp):
            raise ValueError(
                f"Reference and signal total power arrays are different lengths: {len(reftp)} != {len(sigtp)}"
            )

        self._bintable_index = sigtp[0]._bintable_index
        self._scan = sigtp[0]._scan
        self._sigtp = sigtp
        self._reftp = reftp
        # If the user supplied channel= to subbeamnod(), then reftp will already have the correct channel range.
        self._nchan = len(reftp[0]._calibrated[0])
        self._nrows = np.sum([stp.nrows for stp in self._sigtp])
        self._nint = self._nrows
        self._vane = vane
        # Take the first reference scan for each sigtp as the row to use for creating metadata.
        meta_rows = []
        for r in self._sigtp:
            meta_rows.append(r._refoffrows[0])
        meta_rows = list(set(meta_rows))

        self._finish_initialization(calibrate, {"weights": weights}, meta_rows, tscale, zenith_opacity, tcal=tcal)

    def _calc_exposure(self):
        # This is done in calibrate() via assignment.
        pass

    def _calc_delta_freq(self):
        # This is done in calibrate() via assignment.
        pass

    def calibrate(self, **kwargs):  ##SUBBEAMNOD
        """Calibrate the SubBeamNodScan data"""
        nspect = len(self._reftp)
        self._tsys = np.empty(nspect, dtype=float)
        self._exposure = np.empty(nspect, dtype=float)
        self._duration = np.empty(nspect, dtype=float)
        self._delta_freq = np.empty(nspect, dtype=float)
        self._calibrated = np.ma.empty((nspect, self._nchan), dtype=float)

        if self._vane is not None:
            self._tcal[:] = self.get_vane_tcal()

        for i in range(nspect):
            sig = self._sigtp[i].timeaverage(weights=kwargs["weights"])
            ref = self._reftp[i].timeaverage(weights=kwargs["weights"])
            nsmooth = 1.0
            if self._smoothref > 1:
                ref = ref.smooth("box", self._smoothref, decimate=-1)
                nsmooth = self._smoothref
            # Set system temperature.
            self._tsys[i] = ref.meta["WTTSYS"]
            if self._vane is not None:
                self._tsys[i] = self._vane._get_tsys(ref.data, self._tcal[i])
            tsys = self._tsys[i]
            # Combine sig and ref.
            ta = ((sig - ref) / ref).flux.value * tsys
            self._exposure[i] = (
                sig.meta["EXPOSURE"]
                * ref.meta["EXPOSURE"]
                * nsmooth
                / (sig.meta["EXPOSURE"] + ref.meta["EXPOSURE"] * nsmooth)
            )
            self._duration[i] = (
                sig.meta["DURATION"]
                * ref.meta["DURATION"]
                * nsmooth
                / (sig.meta["DURATION"] + ref.meta["DURATION"] * nsmooth)
            )
            self._delta_freq[i] = sig.meta["CDELT1"]
            self._calibrated[i] = ta
