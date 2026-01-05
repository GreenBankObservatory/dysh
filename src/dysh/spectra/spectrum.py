"""
The Spectrum class to contain and manipulate spectra.
"""

import warnings
from copy import deepcopy

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import ITRS, SkyCoord, SpectralCoord, StokesCoord
from astropy.coordinates.spectral_coordinate import attach_zero_velocities
from astropy.io import registry
from astropy.io.fits import BinTableHDU, Column
from astropy.io.fits.verify import VerifyWarning
from astropy.modeling.fitting import LinearLSQFitter

# from astropy.nddata.ccddata import fits_ccddata_writer
from astropy.table import Table
from astropy.time import Time
from astropy.units.quantity import Quantity
from astropy.utils.masked import Masked
from astropy.wcs import WCS, FITSFixedWarning
from ndcube import NDCube
from scipy.stats import anderson
from specutils import Spectrum as Spectrum1D

from dysh.log import logger
from dysh.spectra import core

from ..coordinates import (  # is_topocentric,; topocentric_velocity_to_frame,
    KMS,
    Observatory,
    astropy_convenience_frame_names,
    astropy_frame_dict,
    change_veldef,
    frame_to_label,
    get_velocity_in_frame,
    make_target,
    replace_convention,
    sanitize_skycoord,
    veldef_to_convention,
)
from ..line import SpectralLineSearch
from ..line.search import all_cats
from ..log import HistoricalBase, log_call_to_history, log_call_to_result
from ..plot import specplot as sp
from ..util import (
    docstring_parameter,
    minimum_string_match,
)
from ..util.docstring_manip import copy_docstring
from . import (
    FWHM_TO_STDDEV,
    available_smooth_methods,
    baseline,
    curve_of_growth,
    decimate,
    exclude_to_spectral_region,
    get_spectral_equivalency,
    spectral_region_to_list_of_tuples,
    spectral_region_to_unit,
)

# from astropy.nddata import StdDevUncertainty

# Spectrum attributes to be ignored by Spectrum._copy_attributes
IGNORE_ON_COPY = [
    "_data",
    "_flux",
    "_meta",
    "_mask",
    "_weights",
    "_baseline_model",
    "_plotter",
    "_uncertainty",
    "_unit",
]


class Spectrum(Spectrum1D, HistoricalBase):
    """
    This class contains a spectrum and its attributes. It is built on
    `~specutils.Spectrum` with added attributes like baseline model.
    Note that `~specutils.Spectrum` can contain multiple spectra but
    we probably will not use that because the restriction that it can
    have only one spectral axis conflicts with slight Doppler shifts.
    See `~specutils.Spectrum` for the instantiation arguments.
    """

    @log_call_to_history
    def __init__(self, *args, **kwargs):
        HistoricalBase.__init__(self)
        self._target = kwargs.pop("target", None)
        if self._target is not None:
            self._target = sanitize_skycoord(self._target)
            self._velocity_frame = self._target.frame.name
        else:
            self._velocity_frame = None
        self._observer = kwargs.pop("observer", None)
        _ = kwargs.pop("psf", None)  # Hack to enable rdiv.
        Spectrum1D.__init__(self, *args, **kwargs)
        # Try making a target from meta. If it fails, don't worry about it.
        if False:
            if self._target is None:
                try:
                    self._target = make_target(self.meta)
                    self._velocity_frame = self._target.frame.name
                except Exception:
                    pass
        self._spectral_axis._target = self._target
        if self.velocity_convention is None and "VELDEF" in self.meta:
            self._spectral_axis._doppler_convention = veldef_to_convention(self.meta["VELDEF"])
        if "MEANTSYS" in self.meta and "EXPOSURE" in self.meta and "CDELT1" in self.meta:
            w = core.tsys_weight(self.meta["EXPOSURE"], self.meta["CDELT1"], self.meta["MEANTSYS"])
            self._weights = np.full(self.flux.shape, w)
        else:
            self._weights = np.ones(self.flux.shape)
        #  Get observer quantity from observer location.
        if "DATE-OBS" in self.meta:
            self._obstime = Time(self.meta["DATE-OBS"])
        elif "MJD-OBS" in self.meta:
            self._obstime = Time(self.meta["MJD-OBS"])
        else:
            self._obstime = None
        self._spectral_axis._observer = self.observer
        if self._spectral_axis._observer is not None:
            self._velocity_frame = self._spectral_axis._observer.name

        # if mask is not set via the flux input (NaNs in flux or flux.mask),
        # then set the mask to all False == good
        if self.mask is None:
            self._mask = np.full(np.shape(self.flux), False)
        self._baseline_model = None
        self._subtracted = False
        self._normalized = False
        self._exclude_regions = None
        self._include_regions = None  # do we really need this?
        self._plotter = None
        # `self._resolution` is the spectral resolution in channels.
        # This will be 1 for VEGAS, and >1 for the ACS.
        if "FREQRES" in self.meta and "CDELT1" in self.meta:
            self._resolution = self.meta["FREQRES"] / abs(self.meta["CDELT1"])
        else:
            self._resolution = 1

    def _spectrum_property(self, prop: str):
        """
        Utility method to return a header value as a property.

        Parameters
        ----------
        prop : str
            The property name, _case-sensitive_

        Returns
        -------
        value : Any or Non
            The property value or None if `prop` is not in the header.

        """
        return self.meta.get(prop, None)

    @property
    def nchan(self) -> int:
        """The number of channels in the Spectrum"""
        return len(self.frequency)

    @property
    def weights(self):
        """The channel weights of this spectrum"""
        return self._weights

    @property
    def exclude_regions(self):
        """The baseline exclusion region(s) of this spectrum"""
        return self._exclude_regions

    @property
    def baseline_model(self):
        """Returns the computed baseline model or None if it has not yet been computed."""
        return self._baseline_model

    @property
    def flux(self):
        """
        Converts the stored data and unit and mask into a `~astropy.units.quantity.Quantity` object.

        Returns
        -------
        `~astropy.units.quantity.Quantity`
            Spectral data as a quantity. Masked values are filled with NaN.
        """
        return Masked(self.data * self.unit, mask=self.mask).filled(np.nan)

    @log_call_to_history
    def baseline(self, degree, exclude=None, include=None, color="k", **kwargs):
        # fmt: off
        """
        Compute and optionally remove a baseline.
        The code uses `~astropy.modeling.fitting.Fitter` and `~astropy.modeling.polynomial` to compute the baseline.
        See the documentation for those modules for details.
        This method will set the `baseline_model` attribute to the fitted model function which can be evaluated over a domain.

        Note that `include` and `exclude` are mutually exclusive. If both are present, only `include` will be used.

        Parameters
        ----------
        degree : int
            The degree of the polynomial series, a.k.a. baseline order
        exclude : list of 2-tuples of int or `~astropy.units.quantity.Quantity`, or `~specutils.SpectralRegion`
            List of region(s) to exclude from the fit.  The tuple(s) represent a range in the form [lower,upper], inclusive.

            Examples:

            One channel-based region: [11,51]

            Two channel-based regions: [(11,51),(99,123)].

            One `~astropy.units.quantity.Quantity` region: [110.198*u.GHz,110.204*u.GHz].

            One compound `~specutils.SpectralRegion`: SpectralRegion([(110.198*u.GHz,110.204*u.GHz),(110.196*u.GHz,110.197*u.GHz)]).

            Default: no exclude region

        include : list of 2-tuples of int or `~astropy.units.quantity.Quantity`, or `~specutils.SpectralRegion`
            List of region(s) to include in the fit. The tuple(s) represent a range in the form [lower,upper], inclusive.
            See `exclude` for examples.
        color : str
            The color to plot the baseline model, if remove=False. Can be any type accepted by matplotlib.
        model : str
            One of 'polynomial', 'chebyshev', 'legendre', or 'hermite'
            Default: 'chebyshev'
        fitter : `~astropy.modeling.fitting.Fitter`
            The fitter to use. Default: `~astropy.modeling.fitting.LinearLSQFitter` (with `calc_uncertaintes=True`).
            Be careful when choosing a different fitter to be sure it is optimized for this problem.
        remove : bool
            If True, the baseline is removed from the spectrum.
            If False, the baseline will be computed and overlaid on the spectrum. Default: False
        """
        # fmt: on
        # @todo: Are exclusion regions OR'd with the existing mask? make that an option?
        kwargs_opts = {
            "remove": False,
            "model": "chebyshev",
            "fitter": LinearLSQFitter(calc_uncertainties=True),
        }
        kwargs_opts.update(kwargs)

        # `include` and `exclude` are mutually exclusive, but we allow `include`
        # if `include` is used, transform it to `exclude`.
        if include is not None:
            if exclude is not None:
                logger.warning(f"Warning: ignoring exclude={exclude}")
            exclude = core.include_to_exclude_spectral_region(include, self)
        self._baseline_model = baseline(self, degree, exclude, **kwargs)
        if kwargs_opts["remove"]:
            s = self.subtract(self._baseline_model(self.spectral_axis))
            self._data = s._data
            self._subtracted = True
        if self._plotter is not None:
            if kwargs_opts["remove"]:
                self._plotter._line.set_ydata(self._data)
                self._plotter.clear_overlays(blines=True)
                if not self._plotter._freezey:
                    self._plotter.freey()
            else:
                if self._plotter._xunit == "chan":
                    xval = np.arange(len(self.flux))
                else:
                    xval = self._plotter._sa
                bline_data = self._baseline_model(self.spectral_axis).to(self._plotter._yunit)
                self._plotter._axis.plot(xval, bline_data, c=color, gid="baseline")
                self._plotter.refresh()

    # baseline
    @log_call_to_history
    def undo_baseline(self):
        """
        Undo the most recently computed baseline. If the baseline
        has been subtracted, it will be added back to the data. The `baseline_model`
        attribute is set to None. Exclude regions are untouched.
        """
        if self._baseline_model is None:
            return
        if self._subtracted:
            if self._normalized:
                logger.warning("Cannot undo previously normalized baseline subtraction")
                return
            s = self.add(self._baseline_model(self.spectral_axis))
            self._data = s._data
            if self._plotter is not None:
                self._plotter._line.set_ydata(self._data)
                if not self._plotter._freezey:
                    self._plotter.freey()
        self._baseline_model = None

    @property
    def subtracted(self):
        """Has a baseline model been subtracted?"

        Returns
        -------
            True if a baseline model has been subtracted, False otherwise
        """
        return self._subtracted

    def _set_exclude_regions(self, exclude):
        """
        Set the mask for the regions to exclude.

        Parameters
        ----------
        exclude : list of 2-tuples of int or ~astropy.units.quantity.Quantity, or ~specutils.SpectralRegion
            List of region(s) to exclude from the fit.  The tuple(s) represent a range in the form [lower,upper], inclusive.
            In channel units.

            Examples: One channel-based region: [11,51],
                      Two channel-based regions: [(11,51),(99,123)].
                      One `~astropy.units.quantity.Quantity` region: [110.198*u.GHz,110.204*u.GHz].
                      One compound `~specutils.SpectralRegion`: SpectralRegion([(110.198*u.GHz,110.204*u.GHz),(110.196*u.GHz,110.197*u.GHz)]).

        """
        pass

    def list_to_spectral_region(self, inlist):
        # @todo utility code to convert a input list of channels or quantities to a spectral region with units of self.spectral_axis.unit.
        # This could go in core.py combine this with _set_exclude_regions
        pass

    def bshow(self):
        """Show the baseline model"""
        print(f"baseline model {self._baseline_model}")

    @copy_docstring(sp.SpectrumPlot.plot)
    def plot(self, **kwargs):
        """ """

        if self._plotter is None:
            self._plotter = sp.SpectrumPlot(self, **kwargs)
        self._plotter.plot(**kwargs)
        return self._plotter

    def get_selected_regions(self, unit=None):
        """Get selected regions from plot."""
        if self._plotter is None:
            raise TypeError("No plotter attached to spectrum. Use Spectrum.plot() first.")

        regions = self._plotter.get_selected_regions()

        # If there are no selected regions, tell the user and return.
        if len(regions) == 0:
            warnings.warn(
                "No selected regions.",
                UserWarning,
                stacklevel=2,
            )
            return

        if unit is not None:
            regions = exclude_to_spectral_region(regions, self)
            regions = spectral_region_to_unit(regions, self, unit=unit, append_doppler=True)
            regions = spectral_region_to_list_of_tuples(regions)

        return regions

    @property
    def obstime(self):
        return self._obstime

    @property
    def plotter(self):
        return self._plotter

    def stats(self, roll=0, qac=False):
        """
        Compute some statistics of this `Spectrum`.  The mean, rms, median,
        data minimum and data maximum are calculated.  Note this works
        with slicing, so, e.g.,  `myspectrum[45:153].stats()` will return
        the statistics of the slice.

        Parameters
        ----------
        roll : int
            Return statistics on a 'rolled' array differenced with the
            original array. If there is no correllaton between channels,
            a roll=1 would return an RMS sqrt(2) larger than that of the
            input array. Another advantage of rolled statistics it will
            remove most slow variations, thus RMS/sqrt(2) might be a better
            indication of the underlying RMS.
            Note that for roll > 1, the RMS is already corrected by sqrt(2),
            so they can be directly compared.
            Default: 0
        qac : bool
            If set, the returned simple string contains mean,rms,datamin,datamax
            for easier visual regression. Based on some legacy code.
            Default: False

        Returns
        -------
        stats : dict
            Dictionary consisting of (mean,median,rms,datamin,datamax)
        """

        # note the Spectrum class has special nanXXX functions for most, but not std()
        if roll == 0:
            mean = self.mean()
            median = self.median()
            rms = np.nanstd(self.flux)
            dmin = self.min()
            dmax = self.max()
            npt = len(self.flux)
        else:
            d = self[roll:] - self[:-roll]
            mean = d.mean()
            median = d.median()
            rms = np.nanstd(d.flux) / np.sqrt(2)
            dmin = d.min()
            dmax = d.max()
            npt = len(self.flux) - 2

        if qac:
            out = f"{mean.value} {rms.value} {dmin.value} {dmax.value}"
            return out

        # these two should be the same
        nan1 = np.isnan(self.data).sum()
        nan2 = self.mask.sum()
        if nan1 != nan2:
            logger.warning(f"Warning: {nan1} != {nan2}: inconsistency counters in mask usage")
        elif nan1 > 0:
            logger.info(f"Note: found {nan1} NaN (masked) values")

        out = {"mean": mean, "median": median, "rms": rms, "min": dmin, "max": dmax, "npt": npt, "nan": nan2}

        return out

    def radiometer(self, roll=0):
        """
        Check the radiometer equation, and return the dimensionless ratio of the
        measured vs. expected noise. Generally this number is 1.0 or higher, unless
        for example channels were hanning correlated, measured noise will be lower.

        User is responsible for selecting the channels, via e.g. indexing:

             r1 = sp0[1000:2000].radiometer()

        Parameters
        ----------
        roll : int
             Subtract the data of channel `i+roll` from channel `i`
             before computing the rms. This helps reduce artifacts due
             to bad baselines or channel-to-channel correlations. A
             value of 1 or 2 is recommended.
             See also the `~dysh.spectra.spectrum.Spectrum.roll` function
             where a series of `roll` values can be checked.
             The default is 0 (no roll).

        Returns
        -------
        ratio : real
            The ratio of measured to expected RMS.

        """
        dt = self.meta["EXPOSURE"]
        df = abs(self.meta["CDELT1"])
        tsys = self.meta["TSYS"]
        if roll == 0:
            rms0 = self.stats()["rms"].value
        else:
            rms0 = self.stats(roll=roll)["rms"].value
        rms1 = tsys / np.sqrt(df * dt)
        return rms0 / rms1

    def roll(self, rollmax=1):
        """
        Rolling data to check for channel correllations and channel-to-channel
        correllations. For all roll's from 1 to rollmax the RMS in the rolled
        data is compared to the raw RMS. For well behaved (and baseline subtracted)
        data the ratio of the raw RMS to the rolled RMS should approach 1.0.
        Note: rolling the data is shifting the data by "roll" channels.

        Parameters
        ----------
        rollmax : int
            Roll the data by values 1 through `rollmax`.
            The default is 1.

        Returns
        -------
        ratio : list
            The ratios of raw RMS by the rolled RMS. The list has a length of `rollmax`.
        """
        rms0 = self.stats()["rms"].value
        r = []
        for n in range(1, rollmax + 1):
            r.append(rms0 / self.stats(roll=n)["rms"].value)
        return r

    def snr(self, peak=True, flux=False, rms=None):
        """
        Signal-to-noise (S/N) ratio, measured either in channel or total flux mode.
        Make sure the spectrum has been baseline substracted, or the snr is
        meaningless.
        See also :meth:`~dysh.spectra.Spectrum.sratio`, the signal ratio.


        Parameters
        ----------
        peak : bool
            If True, the largest positive  deviation from the mean is compared to the rms.
            If False, the largest negative deviation from the mean is compared to the rms.
            For normal noise the returned snr value depends on the number of channels
            via the error function.

        flux : bool
            If True, the integrated flux over the spectrum is compared to the expected
            flux given pure noise.
            If False, channel based snr is computed, also controlled by the value of the
            peak in the spectrum.

            See also `Specutils Analysis <https://specutils.readthedocs.io/en/stable/analysis.html>`_.

        rms : None or `~astropy.units.quantity.Quantity`}
            If given, this is the RMS used in the S/N computations. By default it is
            determined from the `Spectrum.stats(roll=1)["rms"]` value of the `Spectrum`.

        Returns
        -------
        ratio : float
            The S/N, either flux or channel based
        """
        # @todo  could check if the data has a baseline solution
        s0 = self.stats()
        s1 = self.stats(roll=1)
        if rms is None:
            rms = s1["rms"]
        if flux:
            snr = np.nansum(self.flux) / (rms * np.sqrt(len(self.flux)))
        elif peak:
            snr = (s0["max"] - s0["mean"]) / rms
        else:
            snr = (s0["mean"] - s0["min"]) / rms
        return snr.value

    def sratio(self, mean=0.0):
        """
        Signal ratio:   (pSum+nSum)/(pSum-nSum)
        Here pSum and nSum are the sum of positive and negative values respectively
        in the spectrum. The `Spectrum` must have been baseline subtracted before using
        this function, otherwise the results will not make sense.

        Parameters
        ----------
        mean : float
            At your own risk, don't use this, should do a baseline subtraction before.
            If not, this could be used as a cheat.

        Returns
        -------
        ratio : real
            The signal ratio, between -1 and 1, 0 being pure noise.
        """
        sp = self.flux.value - mean
        psum = sp[sp > 0.0].sum()
        nsum = sp[sp < 0.0].sum()
        return (psum + nsum) / (psum - nsum)

    def normalness(self):
        """
        Compute the p-value if the noise in a spectrum is gaussian
        using the Anderson-Darling statistic
        The p-value gives the probability that the spectrum is gaussian.
        If p>0.05, the spectrum can be considered gaussian.
        See also "D'Agostino, R. B., & Stephens, M. A. (Eds.). (1986)
        Goodness-of-fit techniques" , table 4.9
        """
        anderson_test = anderson(self.data)
        #   see also  https://github.com/SJVeronese/nicci-package/
        z = anderson_test.statistic
        if z >= 0.6:
            p = np.exp(1.2937 - 5.709 * z + 0.0186 * z**2)
        elif z >= 0.34:
            p = np.exp(0.9177 - 4.279 * z - 1.38 * z**2)
        elif z >= 0.2:
            p = 1 - np.exp(-8.318 + 42.796 * z - 59.938 * z**2)
        else:
            p = 1 - np.exp(-13.436 + 101.14 * z - 223.73 * z**2)
        return p

    @log_call_to_history
    def decimate(self, n):
        """
        Decimate the `Spectrum` by n pixels.

        Parameters
        ----------
        n : int
            Decimation factor of the spectrum by returning every n-th channel.

        Returns
        -------
        s : `Spectrum`
            The decimated `Spectrum`.
        """
        new_data, new_meta = decimate(self.data * self.flux.unit, n, self.meta)
        s = Spectrum.make_spectrum(new_data, meta=new_meta, observer_location="from_meta")
        if self._baseline_model is not None:
            s._baseline_model = None

        return s

    @log_call_to_history
    def smooth(
        self,
        method="hanning",
        width=1,
        decimate=0,
        meta=None,
        mask=None,
        boundary="extend",
        nan_treatment="fill",
        fill_value=np.nan,
        preserve_nan=True,
    ):
        """
        Smooth or Convolve the `Spectrum`, optionally decimating it.

        Default smoothing is hanning.

        Parameters
        ----------
        method : string
            Smoothing method. Valid are: 'hanning', 'boxcar' and
            'gaussian'. Minimum match applies.
            The default is 'hanning'.
        width : int
            Effective width of the convolving kernel.
            For 'hanning' this should be 1, with a 0.25,0.5,0.25 kernel.
            For 'boxcar' an even value triggers an odd one with half the
            signal at the edges, and will thus not reproduce GBTIDL.
            For 'gaussian' this is the FWHM of the final frequency response.
            That is, the smoothed `Spectrum` will have an effective frequency resolution
            equal to CDELT1*`width`.
            The default is 1.
        decimate : int
            Decimation factor of the spectrum by returning every decimate channel.
            -1:   no decimation
            0:    use the width parameter
            >1:   user supplied decimation (use with caution)
        mask : None or `numpy.ndarray`
            A "mask" array.  Shape must match ``array``, and anything that is masked
            (i.e., not 0/`False`) will be set to NaN for the convolution.  If
            `None`, no masking will be performed unless ``array`` is a masked array.
            If ``mask`` is not `None` *and* ``array`` is a masked array, a pixel is
            masked if it is masked in either ``mask`` *or* ``array.mask``.
        boundary : str
            A flag indicating how to handle boundaries:
                * `None`
                    Set the ``result`` values to zero where the kernel
                    extends beyond the edge of the array.
                * 'fill'
                    Set values outside the array boundary to ``fill_value`` (default).
                * 'wrap'
                    Periodic boundary that wrap to the other side of ``array``.
                * 'extend'
                    Set values outside the array to the nearest ``array``
                    value.
        fill_value : float
            The value to use outside the array when using ``boundary='fill'``. Default value is ``NaN``.
        nan_treatment : {'interpolate', 'fill'}
            The method used to handle NaNs in the input ``array``:
                * ``'interpolate'``: ``NaN`` values are replaced with
                  interpolated values using the kernel as an interpolation
                  function. Note that if the kernel has a sum equal to
                  zero, NaN interpolation is not possible and will raise an
                  exception.
                * ``'fill'``: ``NaN`` values are replaced by ``fill_value``
                  prior to convolution.
        preserve_nan : bool
            After performing convolution, should pixels that were originally NaN
            again become NaN?

        Raises
        ------
        Exception
            If no valid smoothing method is given.
        ValueError
            If `width` is less than one.
            If `width` is less than the spectral resolution (in channels).
            If `decimate` is not an integer.

        Returns
        -------
        s : `Spectrum`
            The new, possibly decimated, convolved `Spectrum`.
        """

        valid_methods = available_smooth_methods()
        this_method = minimum_string_match(method, valid_methods)
        if width < 1:
            raise ValueError(f"`width` ({width}) must be >=1.")

        if this_method is None:
            raise Exception(f"smooth({method}): valid methods are {valid_methods}")
        md = np.ma.masked_array(self._data, self.mask)
        if decimate == 0:
            # Take the default decimation by `width`.
            decimate = int(abs(width))
            if not float(width).is_integer():
                logger.info(f"Adjusting decimation factor to be a natural number. Will decimate by {decimate}")
        if this_method == "gaussian":
            if width <= self._resolution:
                raise ValueError(
                    f"`width` ({width} channels) cannot be less than the current resolution ({self._resolution} channels)."
                )
            kwidth = np.sqrt(width**2 - self._resolution**2)  # Kernel effective width.
            stddev = kwidth / FWHM_TO_STDDEV

            new_data, new_meta = core.smooth(
                data=md,
                method=this_method,
                width=stddev,
                ndecimate=decimate,
                meta=self.meta,
                boundary=boundary,
                mask=mask,
                nan_treatment=nan_treatment,
                fill_value=fill_value,
                preserve_nan=preserve_nan,
            )
        else:
            new_data, new_meta = core.smooth(
                data=md,
                method=this_method,
                width=width,
                ndecimate=decimate,
                meta=self.meta,
                boundary=boundary,
                mask=mask,
                nan_treatment=nan_treatment,
                fill_value=fill_value,
                preserve_nan=preserve_nan,
            )

        s = Spectrum.make_spectrum(new_data * self.flux.unit, meta=new_meta, observer_location="from_meta")
        s._baseline_model = self._baseline_model

        # Update the spectral resolution in channel units.
        s._resolution = s.meta["FREQRES"] / abs(s.meta["CDELT1"])

        return s

    def shift(self, s, remove_wrap=True, fill_value=np.nan, method="fft"):
        """
        Shift the `Spectrum` by `s` channels in place.

        Parameters
        ----------
        s : float
            Number of channels to shift the `Spectrum` by.
        remove_wrap : bool
            If `False` keep channels that wrap around the edges.
            If `True` fill channels that wrap with `fill_value`.
        fill_value : float
            If `remove_wrap=True` fill channels that wrapped with this value.
        method : "fft"
            Method used to perform the fractional channel shift.
            "fft" uses a phase shift.
        """

        new_spec = self._copy()
        data = np.ma.masked_where(new_spec.mask, new_spec.data)
        new_data = core.data_shift(data, s, remove_wrap=remove_wrap, fill_value=fill_value, method=method)

        # Update data values.
        new_spec._data = new_data

        # Update metadata.
        new_spec.meta["CRPIX1"] += s

        # Update WCS.
        new_spec.wcs.wcs.crpix[0] += s

        # Update `SpectralAxis` values.
        # Radial velocity needs to be copied by hand.
        radial_velocity = deepcopy(new_spec._spectral_axis._radial_velocity)
        new_spectral_axis_values = new_spec.wcs.spectral.pixel_to_world(np.arange(new_spec.flux.shape[-1]))
        new_spec._spectral_axis = new_spec.spectral_axis.replicate(value=new_spectral_axis_values)
        new_spec._spectral_axis._radial_velocity = radial_velocity

        return new_spec

    def find_shift(self, other, units=None, frame=None):
        """
        Find the shift required to align this `Spectrum` with `other`.

        Parameters
        ----------
        other : `Spectrum`
            Target `Spectrum` to align to.
        units : {None, `astropy.units.quantity.Quantity`}
            Find the shift to align the two `Spectra` in these units.
            If `None`, the `Spectra` will be aligned using the units of
            `other`.
        frame : {None, str}
            Find the shift in this reference frame.
            If `None` will use the frame of `other`.

        Returns
        -------
        shift : float
            Number of channels that this `Spectrum` must be shifted to
            be aligned with `other`.
        """

        if not isinstance(other, Spectrum):
            raise ValueError("`other` must be a `Spectrum`.")

        if frame is not None and frame not in astropy_frame_dict.keys():
            raise ValueError(
                f"`frame` ({frame}) not recognized. Frame must be one of {', '.join(list(astropy_frame_dict.keys()))}"
            )
        else:
            frame = other._velocity_frame

        sa = self.spectral_axis.with_observer_stationary_relative_to(frame)
        tgt_sa = other.spectral_axis.with_observer_stationary_relative_to(frame)

        if units is None:
            units = tgt_sa.unit

        sa = sa.to(units)
        tgt_sa = tgt_sa.to(units)

        cdelt1 = sa[1] - sa[0]
        shift = ((sa[0] - tgt_sa[0]) / cdelt1).value

        return shift

    def align_to(self, other, units=None, frame=None, remove_wrap=True, fill_value=np.nan, method="fft"):
        """
        Align the `Spectrum` with respect to `other`.

        Parameters
        ----------
        other : `Spectrum`
            Target `Spectrum` to align to.
        units : {None, `astropy.units.quantity.Quantity`}
            Find the shift to align the two `Spectra` in these units.
            If `None`, the `Spectra` will be aligned using the units of
            `other`.
        frame : {None, str}
            Find the shift in this reference frame.
            If `None` will use the frame of `other`.
        remove_wrap : bool
            If `True` allow spectrum to wrap around the edges.
            If `False` fill channels that wrap with `fill_value`.
        fill_value : float
            If `wrap=False` fill channels that wrapped with this value.
        method : "fft"
            Method used to perform the fractional channel shift.
            "fft" uses a phase shift.
        """

        s = self.find_shift(other, units=units, frame=frame)
        return self.shift(s, remove_wrap=remove_wrap, fill_value=fill_value, method=method)

    @property
    def equivalencies(self):
        """Get the spectral axis equivalencies that can be used in converting the axis
        between km/s and frequency or wavelength"""
        equiv = deepcopy(u.spectral())
        sa = self.spectral_axis
        if sa.doppler_rest is not None:
            rfq = sa.doppler_rest
        elif "RESTFREQ" in self.meta:
            cunit1 = self.meta.get("CUNIT1", self.wcs.wcs.cunit[0])
            # @todo this could be done with a dict str->function
            rfq = self.meta["RESTFREQ"] * cunit1  # WCS wants no E
        else:
            rfq = None
        if rfq is not None:
            equiv.extend(get_spectral_equivalency(rfq, self.velocity_convention))
        return equiv

    @property
    def target(self):
        """
        The target object of this spectrum.

        Returns
        -------
        target : `~astropy.coordinates.sky_coordinate.SkyCoord`
            The sky coordinate object
        """
        return self._target

    @property
    def tscale(self):
        """
        The descriptive brightness unit of the data. One of
            - 'Raw' : raw value, e.g., count
            - 'Ta'  : Antenna Temperature in K
            - 'Ta*' : Antenna temperature corrected to above the atmosphere in K
            - 'Flux': flux density in Jansky

        Returns
        -------
        tscale : str or None
            Brightness unit string. If there is no TSCALE in the header, None is returned.

        """
        return self._spectrum_property("TSCALE")

    @property
    def tscale_fac(self):
        """
        The factor by which the data have been scale from antenna temperature to corrected antenna temperature
        or flux density. This is the average of the values by which the integrations in the spectrum have been scaled.

        Returns
        -------
        tscale_fac : float or None
            The scale factor. If there is no TSCALFAC in the header, None is returned.

        """
        return self._spectrum_property("TSCALFAC")

    @property
    def observer(self):
        """
        Returns
        -------
            observer : `~astropy.coordinates.BaseCoordinateFrame` or derivative
            The coordinate frame of the observer if present.
        """
        return self._observer

    @property
    def velocity_frame(self) -> str:
        """String representation of the velocity frame"""
        return self._velocity_frame

    @property
    def doppler_convention(self) -> str:
        """String representation of the velocity (Doppler) convention"""
        return self.velocity_convention

    @property
    def rest_value(self) -> Quantity:
        """Rest frequency used in velocity conversions.

        Returns
        -------
        ~astropy.units.quantity.Quantity.Quantity
            The rest frequency as a Quantity object
        """
        return self.spectral_axis.doppler_rest

    @rest_value.setter
    def rest_value(self, value: Quantity):
        """
        "Set the rest frequency property and update the `Spectrum` metadata.

        Parameters
        ----------
        value : ~astropy.units.quanityt.Quantity
            A frequency-like quantity.

        """
        self._spectral_axis._doppler_rest = value
        self.meta["RESTFREQ"] = value.to("Hz").value
        self.meta["RESTFRQ"] = value.to("Hz").value

    def axis_velocity(self, unit=KMS):
        """Get the spectral axis in velocity units.
        *Note*: This is not the same as `Spectrum.velocity`, which includes the source radial velocity.

        Parameters
        ----------
        unit : `~astropy.units.quantity.Quantity` or str that can be converted to Quantity
                The unit to which the axis is to be converted
        Returns
        -------
        velocity : `~astropy.units.quantity.Quantity.Quantity`
                The converted spectral axis velocity
        """
        return self._spectral_axis.to(unit)

    def velocity_axis_to(self, unit=KMS, toframe=None, doppler_convention=None):
        """
        Convert the spectral axis to `unit` in `toframe` using `doppler_convention`
        if converting from frequency/wavelength to velocity.

        Parameters
        ----------
        unit : `~astropy.units.quantity.Quantity` or str that can be converted to `~astropy.units.quantity.Quantity`
            The unit to which the spectral axis is to be converted.
        toframe : str
            The coordinate frame to convert to, e.g. 'hcrs', 'icrs'.
        doppler_convention : None or {'optical', 'radio', 'relativistic'}
            The Doppler convention to use when converting to velocity.
            One of 'optical', 'radio', or 'relativistic'.

        Returns
        -------
        velocity : `~astropy.units.quantity.Quantity`
            The converted spectral axis in units of `unit`.
        """
        if toframe is not None and toframe != self.velocity_frame:
            s = self.with_frame(toframe)
        else:
            s = self
        if doppler_convention is not None:
            return s._spectral_axis.to(unit=unit, doppler_convention=doppler_convention).to(unit)
        else:
            return s.axis_velocity(unit)

    def get_velocity_in_frame(self, toframe: str) -> Quantity:
        """Compute the radial velocity of the `Spectrum.target` in a new velocity frame.
        See :meth:`~dysh.coordinates.core.get_velocity_in_frame`.

        Parameters
        ----------
        toframe : str
            The coordinate frame to convert to, e.g. 'hcrs', 'icrs'.

        Returns
        -------
        radial_velocity : `~astropy.units.quantity.Quantity`
            The radial velocity of the source in `toframe`
        """
        if self._target is None:
            raise Exception("Can't calculate velocity because Spectrum.target is None")
        return get_velocity_in_frame(self._target, toframe, self._observer, self._obstime)

    @log_call_to_history
    def set_frame(self, toframe):
        """Set the sky coordinate and doppler tracking reference frame of this Spectrum. The header 'CTYPE1' will be changed accordingly.

        To make a copy of this Spectrum with new coordinate referece frmae instead, use `with_frame`.

        Parameters
        ----------
        toframe : str, `~astropy.coordinates.BaseCoordinateFrame`, or `~astropy.coordinates.sky_coordinate.SkyCoord`
            The coordinate reference frame identifying string, as used by astropy, e.g. 'hcrs', 'icrs',
            an actual coordinate system instance, or a sky coordinate instance.
        """

        tfl = toframe
        if isinstance(toframe, str):
            tfl = toframe.lower()
            tfl = astropy_convenience_frame_names.get(tfl, tfl)
            if "itrs" in tfl:
                if isinstance(self._observer, ITRS):
                    return  # nothing to be done, we already have the correct axis
                raise ValueError(
                    "For topographic or ITRS coordaintes, you must supply a full astropy Coordinate instance."
                )
            elif self._velocity_frame == tfl:
                return  # the frame is already the requested frame

        self._spectral_axis = self._spectral_axis.with_observer_stationary_relative_to(tfl)
        self._observer = self._spectral_axis.observer
        if isinstance(tfl, str):
            self._velocity_frame = tfl
        else:
            self._velocity_frame = tfl.name
        # While it is incorrect to change CTYPE1, it is reasonable to change VELDEF.
        # SDFITS defines CTYPE1 as always being the TOPO frequency.
        # See Issue #373 on GitHub.
        self.meta["VELDEF"] = change_veldef(self.meta["VELDEF"], self._velocity_frame)

    def with_frame(self, toframe):
        """Return a copy of this `Spectrum` with a new coordinate reference frame.

        Parameters
        ----------
        toframe : str, `~astropy.coordinates.BaseCoordinateFrame`, or `~astropy.coordinates.sky_coordinate.SkyCoord`
            The coordinate reference frame identifying string, as used by astropy, e.g. 'hcrs', 'icrs',
            an actual coordinate system instance, or a sky coordinate instance.

        Returns
        -------
        spectrum : `Spectrum`
            A new `Spectrum` object with the rest frame set to `toframe`.
        """

        s = self._copy()
        s.set_frame(toframe)
        return s

    @log_call_to_history
    def set_convention(self, doppler_convention):
        """Set the velocity convention of this `Spectrum`.  The spectral axis of this `Spectrum` will be replaced
        with a new spectral axis with the input velocity convention.  The header 'VELDEF' value will
        be changed accordingly.

        To make a copy of this `Spectrum` with a new velocity convention instead, use `Spectrum.with_velocity_convention`.

        Parameters
        ----------
        doppler_convention : str
            The velocity convention, one of 'radio', 'optical', 'relativistic'

        """
        # replicate() gives the same asnwer as
        # self._spectral_axis.to(unit=self._spectral_axis.unit, doppler_convention=doppler_convention)
        new_sp_axis = self.spectral_axis.replicate(doppler_convention=doppler_convention)
        self._spectral_axis = new_sp_axis
        self.meta["VELDEF"] = replace_convention(self.meta["VELDEF"], doppler_convention)

    def with_velocity_convention(self, doppler_convention):
        """Returns a copy of this `Spectrum` with the input velocity convention.  The header 'VELDEF' value will
        be changed accordingly.

        Parameters
        ----------
        doppler_convention : str
            The velocity convention, one of 'radio', 'optical', 'relativistic'

        Returns
        -------
        spectrum : `Spectrum`
            A new `Spectrum` object with `doppler_convention` as its Doppler convention.
        """
        s = self._copy(velocity_convention=doppler_convention)
        s.set_convention(doppler_convention)
        return s

    def savefig(self, file, **kwargs):
        """Save the plot"""
        raise Exception("savefig() has been moved to the SpecPlot class")

    def _write_table(self, fileobj, format, **kwargs):
        """
        Write this `Spectrum` as an ~astropy.table.Table. The output columns will be
            - frequency or velocity - in  the Spectrum's spectral axis units, or converted with `xaxis_unit` keyword
            - flux - in the Spectrum's flux units, or converted with `yaxis_unit` keyword
            - uncertainty - The channel-based uncertainty or zeroes if uncertainty was not defined
            - weights - the channel weights
            - mask - 0 is unmasked ('good'), 1 is masked ('bad')
            - baseline model value - the value of the baseline model at each x axis point, regardless of whether it
              has been subtracted from the Spectrum or not. This will be all zeroes if no baseline model
              has been computed.

        Parameters
        ----------
        fileobj : str, file-like or `pathlib.Path`
            File to write to.  If a file object, must be opened in a
            writeable mode.
        format : str
            The output format. Must be a format supported by ~astropy.table.Table.write
        kwargs : variable
            Additional keyword arguments supported by ~astropy.table.Table.write

        """
        # @todo support xaxis_unit, yaxis_unit, flux_unit
        # xaxis_unit : str or :class:`astropy.units.Unit`
        # yaxis_unit : str or :class:`astropy.units.Unit`
        # flux_unit : str or :class:`astropy.units.Unit`
        if self._baseline_model is None:
            bl = np.zeros_like(self.flux)
            bldesc = "Fitted baseline value at given channel (was not defined)"
        else:
            bl = self._baseline_model(self._spectral_axis)
            bldesc = "Fitted baseline value at given channel"
        mask = self.mask * 1
        if self.uncertainty is None:
            unc = np.zeros(self.flux.shape)  # , mask=False)
            udesc = "Flux uncertainty (was not defined)"
        else:
            unc = self.uncertainty.quantity
            udesc = "Flux uncertainty"
        outarray = [self.spectral_axis, self.flux, unc, self.weights, mask, bl]
        description = [
            "Spectral axis",
            "Flux",
            udesc,
            "Channel weights",
            "Mask 0=unmasked, 1=masked",
            bldesc,
        ]
        # remove FITS reserve keywords
        meta = deepcopy(self.meta)
        meta.pop("NAXIS1", None)
        meta.pop("TDIM7", None)
        meta.pop("TUNIT7", None)
        meta["HISTORY"] = self.history  # use property not _history to ensure ascii
        meta["COMMENT"] = self.comments
        if format == "mrt":
            ulab = "e_flux"  # MRT convention that error on X is labeled e_X
        else:
            ulab = "uncertainty"
        outnames = ["spectral_axis", "flux", ulab, "weight", "mask", "baseline"]
        if "ipac" in format:
            # IPAC format wants a crazy dictionary style.
            d = {}
            for k, v in meta.items():
                d[k] = {"value": v}
            t = Table(outarray, names=outnames, meta={"keywords": d}, descriptions=description)
        else:
            t = Table(outarray, names=outnames, meta=meta, descriptions=description)
        # for now ignore complaints about keywords until we clean them up.
        # There are some that are more than 8 chars that should be fixed in GBTFITSLOAD
        warnings.simplefilter("ignore", VerifyWarning)
        t.write(fileobj, format=format, **kwargs)

    def _copy(self, **kwargs):
        """
        Perform deep copy operations on each attribute of the ``Spectrum``
        object.
        This overrides the ``specutils.Spectrum`` method so that
        target and observer attributes get copied.
        """
        alt_kwargs = dict(
            flux=deepcopy(self.flux),
            spectral_axis=deepcopy(self.spectral_axis),
            uncertainty=deepcopy(self.uncertainty),
            wcs=deepcopy(self.wcs),
            mask=deepcopy(self.mask),
            meta=deepcopy(self.meta),
            unit=deepcopy(self.unit),
            velocity_convention=deepcopy(self.velocity_convention),
            rest_value=deepcopy(self.rest_value),
            target=deepcopy(self.target),
            observer=deepcopy(self.observer),
        )

        alt_kwargs.update(kwargs)

        s = self.__class__(**alt_kwargs)
        # s.topocentric = is_topocentric(meta["CTYPE1"])  # GBT-specific to use CTYPE1 instead of VELDEF
        return s

    @classmethod
    def fake_spectrum(cls, nchan=1024, seed=None, normal=True, use_wcs=True, **kwargs):
        """
        Create a fake spectrum with gaussian noise, useful for simple testing.
        A default header is created, which may be modified with kwargs.

        Parameters
        ----------
        nchan : int, optional
            Number of channels. The default is 1024.
        seed : {None, int, array_like[ints], `numpy.random.SeedSequence`, `numpy.random.BitGenerator`, `numpy.random.Generator`}, optional
            A seed to initialize the `BitGenerator`. If None, then fresh, unpredictable entropy will be pulled from the OS.
            If an int or array_like[ints] is passed, then all values must be non-negative and will be passed to
            `SeedSequence` to derive the initial `BitGenerator` state. One may also pass in a `SeedSequence` instance.
            Additionally, when passed a `BitGenerator`, it will be wrapped by `Generator`. If passed a `Generator`, it will
            be returned unaltered.
            The default is `None`.
        normal : bool, optional
            If set, the noise is distributed normal with a mean 0.1 and dispersion 0.1.
            If False, the noise is uniformly distributed between 0 and 1.
            The default is True.
        use_wcs : bool, optional
            If set, create a WCS object from the metadata.
            This is computationally expensive, so setting it to False can speed things up when no WCS is needed.
            The default is True.
        **kwargs: dict or key=value
            Metadata to put in the header.  If the key exists already in
            the default header, it will be replaced. Otherwise the key and value will be
            added to the header. Keys are case insensitive.

        Returns
        -------
        spectrum : `Spectrum`
            The spectrum object.
        """

        rng = np.random.default_rng(seed)
        if normal:
            data = rng.normal(0.1, 0.1, nchan) * u.K
        else:
            data = rng.random(nchan) * u.K
        meta = {
            "OBJECT": "NGC2415",
            "BANDWID": 23437500.0,
            "DATE-OBS": "2021-02-10T07:38:37.50",
            "DURATION": 0.9982445,
            "EXPOSURE": 44.949832229522286,  # fixed by radiometer equation
            "TSYS": 17.930595470605255,
            "CTYPE1": "FREQ-OBS",
            "CRVAL1": 1402544936.7749996,
            "CRPIX1": float(nchan) / 2.0,
            "CDELT1": -715.2557373046875,
            "CTYPE2": "RA",
            "CRVAL2": 114.23878994411744,
            "CTYPE3": "DEC",
            "CRVAL3": 35.24315395841497,
            "CRVAL4": -6,
            "OBSERVER": "A. Dysh User",
            "OBSID": "unknown",
            "SCAN": 152,
            "OBSMODE": "OnOff:PSWITCHON:TPWCAL",
            "FRONTEND": "Rcvr1_2",
            "TCAL": 1.4551637172698975,
            "VELDEF": "OPTI-HEL",
            "VFRAME": 15264.39118499772,
            "RVSYS": 3775382.910954342,
            "OBSFREQ": 1402544936.7749996,
            "LST": 42101.90296246562,
            "AZIMUTH": 285.9514963267411,
            "ELEVATIO": 42.100623613548194,
            "TAMBIENT": 270.4,
            "PRESSURE": 696.2290227048372,
            "HUMIDITY": 0.949,
            "RESTFREQ": 1420405751.786,
            "FREQRES": 715.2557373046875,
            "EQUINOX": 2000.0,
            "RADESYS": "FK5",
            "TRGTLONG": 114.2375,
            "TRGTLAT": 35.24194444444444,
            "SAMPLER": "A2_0",
            "FEED": 1,
            "SRFEED": 0,
            "FEEDXOFF": 0.0,
            "FEEDEOFF": 0.0,
            "SUBREF_STATE": 1,
            "SIDEBAND": "L",
            "PROCSEQN": 1,
            "PROCSIZE": 2,
            "PROCSCAN": "ON",
            "PROCTYPE": "SIMPLE",
            "LASTON": 152,
            "LASTOFF": 0,
            "TIMESTAMP": "2021_02_10_07:38:37",
            "QD_BAD": -1,
            "QD_METHOD": "",
            "VELOCITY": 3784000.0,
            "DOPFREQ": 1420405751.786,
            "ADCSAMPF": 3000000000.0,
            "VSPDELT": 65536.0,
            "VSPRVAL": 19.203125,
            "VSPRPIX": 16385.0,
            "SIG": "T",
            "CAL": "F",
            "CALTYPE": "LOW",
            "CALPOSITION": "Unknown",
            "IFNUM": 0,
            "PLNUM": 0,
            "FDNUM": 0,
            "HDU": 1,
            "BINTABLE": 0,
            "ROW": 2,
            "DATE": "2022-02-14T17:25:01",
            "ORIGIN": "NRAO Green Bank",
            "TELESCOP": "NRAO_GBT",
            "INSTRUME": "VEGAS",
            "SDFITVER": "sdfits ver1.22",
            "FITSVER": "1.9",
            "CTYPE4": "STOKES",
            "PROJID": "TGBT21A_501_11",
            "BACKEND": "VEGAS",
            "SITELONG": -79.83983,
            "SITELAT": 38.43312,
            "SITEELEV": 824.595,
            "EXTNAME": "SINGLE DISH",
            "FITSINDEX": 0,
            "PROC": "OnOff",
            "OBSTYPE": "PSWITCHON",
            "SUBOBSMODE": "TPWCAL",
            "INTNUM": 0,
            "CUNIT1": "Hz",
            "CUNIT2": "deg",
            "CUNIT3": "deg",
            "RESTFRQ": 1420405751.786,
            "MEANTSYS": 17.16746070048293,
            "WTTSYS": 17.16574907094451,
            "TSCALE": "Ta*",
            "TSCALFAC": 0.54321,
            "TUNIT7": "K",
            "BUNIT": "K",
            "AP_EFF": 0.71,
            "SURF_ERR": 230.0,
            "SE_UNIT": "micron",
            "TAU_Z": 0.08,
        }
        for k, v in kwargs.items():
            meta[k.upper()] = v
        # @todo   fix for radiometer equation"EXPOSURE" "TSYS": "CDELT1"
        return Spectrum.make_spectrum(data, meta, observer_location=Observatory["GBT"], use_wcs=use_wcs)

    # @todo allow observer or observer_location.  And/or sort this out in the constructor.
    @classmethod
    def make_spectrum(cls, data, meta, use_wcs=True, observer_location=None, observer=None):
        # , shift_topo=False):
        """Factory method to create a `Spectrum` object from a data and header.  The the data are masked,
        the `Spectrum` mask will be set to the data mask.

        Parameters
        ----------
        data :  `~numpy.ndarray`
            The data array. See `~specutils.Spectrum`
        meta : dict
            The metadata, typically derived from an SDFITS header.
            Required items in `meta` are 'CTYPE[123]','CRVAL[123]', 'CUNIT[123]', 'VELOCITY', 'EQUINOX', 'RADESYS'
        use_wcs : bool
            If True, create a WCS object from `meta`
            Default: True
        observer_location : `~astropy.coordinates.EarthLocation` or str
            Location of the observatory. See `~dysh.coordinates.Observatory`.
            This will be transformed to `~astropy.coordinates.ITRS` using the time of observation DATE-OBS or MJD-OBS in `meta`.
            If this parameter is given the special str value 'from_meta', then an observer_location
            will be created from SITELONG, SITELAT, and SITEELEV in the meta dictionary.
        observer : `~astropy.coordinates.BaseCoordinateFrame`
            Coordinate frame for the observer.
            Will be ignored if DATE-OBS or MJD-OBS are present in `meta` and `observer_location` is not `None`.

        Returns
        -------
        spectrum : `Spectrum`
            The spectrum object
        """
        # @todo generic check_required method since I now have this code in two places (coordinates/core.py).
        # @todo requirement should be either DATE-OBS or MJD-OBS, but make_target() needs to be updated
        # in that case as well.
        _required = set(
            [
                "CRVAL1",
                "CRVAL2",
                "CRVAL3",
                "CTYPE1",
                "CTYPE2",
                "CTYPE3",
                "CUNIT1",
                "CUNIT2",
                "CUNIT3",
                "VELOCITY",
                "EQUINOX",
                "RADESYS",
                "DATE-OBS",
                "VELDEF",
                "RESTFRQ",
            ]
        )

        # Otherwise we change the input meta, which could lead to confusion.
        _meta = deepcopy(meta)

        # RADECSYS is also a valid column name. See issue #287
        # https://github.com/GreenBankObservatory/dysh/issues/287
        if "RADECSYS" in _meta.keys() and "RADESYS" not in _meta.keys():
            _meta["RADESYS"] = deepcopy(_meta["RADECSYS"])
            del _meta["RADECSYS"]

        missing = set(_required).difference(set(_meta.keys()))
        if len(missing) > 0:
            raise ValueError(f"Header (meta) is missing one or more required keywords: {missing}")

        # @todo WCS is expensive.
        # Possibly figure how to calculate spectral_axis instead.
        # @todo allow fix=False in WCS constructor?
        if use_wcs:
            with warnings.catch_warnings():
                # Skip warnings about DATE-OBS being converted to MJD-OBS.
                warnings.filterwarnings("ignore", category=FITSFixedWarning)
                # Skip warnings FITS keywords longer than 8 chars or containing
                # illegal characters (like _).
                warnings.filterwarnings("ignore", category=VerifyWarning)
                wcs_meta_keys = [
                    "CRPIX1",
                    "CTYPE1",
                    "CDELT1",
                    "CRVAL1",
                    "CUNIT1",
                    "CRVAL2",
                    "CTYPE2",
                    "CUNIT2",
                    "CRVAL3",
                    "CTYPE3",
                    "CUNIT3",
                    "CTYPE4",
                    "CRVAL4",
                    "DATE-OBS",
                ]
                try:
                    wcs_meta = {k: _meta[k] for k in wcs_meta_keys}
                except KeyError as exc:
                    raise KeyError(f"Missing item for {exc} in meta.") from exc
                wcs = WCS(header=wcs_meta)
                # It would probably be safer to add NAXISi to meta.
                if wcs.naxis > 3:
                    wcs.array_shape = (0, 0, 0, len(data))
                # For some reason these aren't identified while creating the WCS object.
                if "SITELONG" in _meta.keys():
                    wcs.wcs.obsgeo[:3] = _meta["SITELONG"], _meta["SITELAT"], _meta["SITEELEV"]
                # Reset warnings.
        else:
            wcs = None
        target = make_target(_meta)
        vc = veldef_to_convention(_meta["VELDEF"])

        # Define an observer as needed.
        if observer is not None:
            obsitrs = observer
        elif observer_location is not None and (_meta.get("DATE-OBS") or _meta.get("MJD-OBS")) is not None:
            obstime = Time(_meta.get("DATE-OBS") or _meta.get("MJD-OBS"))
            if observer_location == "from_meta":
                try:
                    observer_location = Observatory.get_earth_location(
                        _meta["SITELONG"], _meta["SITELAT"], _meta["SITEELEV"]
                    )
                except KeyError as ke:
                    raise Exception(f"Not enough info to create observer_location: {ke}")  # noqa: B904
            obsitrs = SpectralCoord._validate_coordinate(
                attach_zero_velocities(observer_location.get_itrs(obstime=obstime))
            )
        else:
            logger.warning(
                "'meta' does not contain DATE-OBS or MJD-OBS. Spectrum won't be convertible to certain coordinate"
                " frames"
            )
            obsitrs = None

        if hasattr(data, "mask"):
            s = cls(
                flux=data,
                wcs=wcs,
                meta=_meta,
                velocity_convention=vc,
                radial_velocity=target.radial_velocity,
                rest_value=_meta["RESTFRQ"] * u.Hz,
                observer=obsitrs,
                target=target,
                mask=data.mask,
            )
        else:
            s = cls(
                flux=data,
                wcs=wcs,
                meta=_meta,
                velocity_convention=vc,
                radial_velocity=target.radial_velocity,
                rest_value=_meta["RESTFRQ"] * u.Hz,
                observer=obsitrs,
                target=target,
            )
        # For some reason, Spectrum.spectral_axis created with WCS do not inherit
        # the radial velocity. In fact, they get no radial_velocity attribute at all!
        # This method creates a new spectral_axis with the given radial velocity.
        if observer_location is None and observer is None:
            s.set_radial_velocity_to(target.radial_velocity)  # open
        return s

    def _arithmetic_apply(self, other, op, handle_meta, **kwargs):
        if isinstance(other, NDCube):
            result = op(other, **{"handle_meta": handle_meta}, **kwargs)
        else:
            result = op(other, **{"handle_meta": handle_meta, "meta_other_meta": False}, **kwargs)
        self._copy_attributes(result)
        return result

    def _copy_attributes(self, other):
        """
        Copy `Spectrum` attributes after
        an arithmetic operation.
        Only copy attributes that are not modified by the arithmetic.
        I.e., do not copy the "_flux" attribute.
        """
        for k, v in vars(self).items():
            if k not in IGNORE_ON_COPY:
                vars(other)[k] = deepcopy(v)

    def __add__(self, other):
        op = self.add
        handle_meta = self._add_meta
        if not isinstance(other, (NDCube, u.Quantity)):
            try:
                other = u.Quantity(other, unit=self.unit)
            except TypeError:
                return NotImplemented
        result = self._arithmetic_apply(other, op, handle_meta)
        return result

    __radd__ = __add__

    def __sub__(self, other):
        op = self.subtract
        handle_meta = self._add_meta
        if not isinstance(other, (NDCube, u.Quantity)):
            try:
                other = u.Quantity(other, unit=self.unit)
            except TypeError:
                return NotImplemented
        result = self._arithmetic_apply(other, op, handle_meta)
        return result

    def __rsub__(self, other):
        return -1 * (self - other)

    def __mul__(self, other):
        op = self.multiply
        handle_meta = self._mul_meta
        if not isinstance(other, NDCube):
            other = u.Quantity(other)
        result = self._arithmetic_apply(other, op, handle_meta)
        return result

    __rmul__ = __mul__

    def __div__(self, other):
        op = self.divide
        handle_meta = self._div_meta
        if not isinstance(other, NDCube):
            other = u.Quantity(other)
        result = self._arithmetic_apply(other, op, handle_meta)
        return result

    def __truediv__(self, other):
        op = self.divide
        handle_meta = self._div_meta
        if not isinstance(other, NDCube):
            other = u.Quantity(other)
        result = self._arithmetic_apply(other, op, handle_meta)
        return result

    def __rtruediv__(self, other):
        op = self.divide
        handle_meta = self._rdiv_meta
        if not isinstance(other, NDCube):
            other = u.Quantity(other)
        result = self._arithmetic_apply(
            other, op, handle_meta, operand2=self, compare_wcs="no, but I don't like your defaults", wcs_use_self=False
        )
        # Set compare_wcs to something, so self._arithmetic_wcs is used.
        return result

    def _add_meta(self, operand, operand2, **kwargs):
        kwargs.setdefault("other_meta", True)
        meta = deepcopy(operand)
        if kwargs["other_meta"]:
            meta["EXPOSURE"] = operand["EXPOSURE"] + operand2["EXPOSURE"]
            meta["DURATION"] = operand["DURATION"] + operand2["DURATION"]

        return meta

    def _mul_meta(self, operand, operand2=None, **kwargs):
        # TBD
        return deepcopy(operand)

    def _div_meta(self, operand, operand2=None, **kwargs):
        # TBD
        return deepcopy(operand)

    def _rdiv_meta(self, operand, operand2=None, **kwargs):
        return deepcopy(operand2)

    def _arithmetic_wcs(self, operation, operand, compare_wcs, **kwargs):
        # Had to overrride
        # astropy/nddata/mixins/ndarithmetic.NDArithmeticMixin._arithmetic_wcs
        # since it does not provide enough flexibility to take the wcs from
        # operand2. Notice that this function is called with self as either
        # self or operand, and operand as either operand or operand2 depending
        # on the type of operand2 in NDArithmeticMixin._prepare_then_do_arithmetic.
        kwarg_opts = {"use_self": True}
        kwarg_opts.update(kwargs)
        if kwarg_opts["use_self"]:
            return deepcopy(self.wcs)
        else:
            return deepcopy(operand.wcs)

    def __getitem__(self, item):
        def q2idx(q, wcs, spectral_axis, coo, sto):
            """Quantity to index."""
            if "velocity" in u.get_physical_type(q):
                # Convert from velocity to channel.
                idx = vel2idx(q, wcs, spectral_axis, coo, sto)
            elif "length" in u.get_physical_type(q):
                # Convert wavelength to channel.
                idx = wav2idx(q, wcs, spectral_axis, coo, sto)
            elif "frequency" in u.get_physical_type(q):
                # Convert from frequency to channel.
                idxs = wcs.world_to_pixel(coo, q, sto)
                idx = int(np.round(idxs[0]))
            return idx

        def vel2idx(vel, wcs, spectral_axis, coo, sto):
            eq = get_spectral_equivalency(spectral_axis.doppler_rest, spectral_axis.doppler_convention)
            # Make `vel` a `SpectralCoord`.
            vel_sp = spectral_axis.to(unit=vel.unit).replicate(value=vel.value, unit=vel.unit)
            with u.set_enabled_equivalencies(eq):
                idxs = wcs.world_to_pixel(coo, vel_sp, sto)
            return int(np.round(idxs[0]))

        def wav2idx(wav, wcs, spectral_axis, coo, sto):
            wav_sp = spectral_axis.to(unit=wav.unit).replicate(value=wav.value, unit=wav.unit)
            with u.set_enabled_equivalencies(u.spectral()):
                idxs = wcs.world_to_pixel(coo, wav_sp, sto)
            return int(np.round(idxs[0]))

        # Only slicing along the spectral coordinate.
        # Assumes that the WCS is in frequency.
        if isinstance(item, slice):
            if item.step is not None and item.step != 1:
                raise NotImplementedError(
                    "Cannot slice with a step other than 1. Try smoothing and decimating if you want to downsample."
                )

            spectral_axis = self.spectral_axis
            # Use the WCS to convert from world to pixel values.
            wcs = self.wcs
            # We need a sky location to convert incorporating velocity shifts.
            coo = SkyCoord(
                wcs.wcs.crval[wcs.wcs.lng] * wcs.wcs.cunit[wcs.wcs.lng],
                wcs.wcs.crval[wcs.wcs.lat] * wcs.wcs.cunit[wcs.wcs.lat],
                frame="fk5",
            )
            # Same for the Stokes axis.
            sto = StokesCoord(0)
            start = item.start
            stop = item.stop

            # Start.
            if isinstance(start, u.Quantity):
                start_idx = q2idx(start, wcs, spectral_axis, coo, sto)
            else:
                start_idx = start

            # Stop.
            if isinstance(stop, u.Quantity):
                stop_idx = q2idx(stop, wcs, spectral_axis, coo, sto)
            else:
                stop_idx = stop

        else:
            raise NotImplementedError("Use a slice for slicing.")

        # WCS slicing does not do well if the start is negative.
        if start is not None and not isinstance(start, u.Quantity) and start < 0:
            start_idx = len(self.data) + start

        # Sort the start and stop indices if they are not None nor
        # negative integers.
        if (
            start is not None
            and stop is not None
            and not (
                (not isinstance(start, u.Quantity) and start < 0) or (not isinstance(stop, u.Quantity) and stop < 0)
            )
        ):
            start_idx, stop_idx = np.sort([start_idx, stop_idx])

        # Slicing uses NumPY ordering by default.
        sliced_wcs = wcs[0:1, 0:1, 0:1, start_idx:stop_idx]

        new_flux = self.flux[start_idx:stop_idx]

        # Update meta.
        meta = self.meta.copy()
        head = sliced_wcs.to_header()
        for k in ["CRPIX1", "CRVAL1"]:
            meta[k] = head[k]
        meta["BANDWID"] = abs(meta["CDELT1"]) * len(new_flux)  # Hz

        # New Spectrum.
        new_spectrum = self.make_spectrum(
            Masked(new_flux, self.mask[start_idx:stop_idx]),
            meta=meta,
            observer_location=Observatory[meta["TELESCOP"]],
        )
        new_spectrum._weights = self._weights[start_idx:stop_idx]
        return new_spectrum

    @log_call_to_result
    def average(self, spectra, weights: str | np.ndarray | None = "tsys", align=False):
        r"""
        Average this `Spectrum` with `spectra`.
        The resulting `average` will have an exposure equal to the sum of the exposures,
        and coordinates and system temperature equal to the weighted average of the coordinates and system temperatures.

        Parameters
        ----------
        spectra : list of `Spectrum`
            Spectra to be averaged. They must have the same number of channels.
            No checks are done to ensure they are aligned.
        weights: None, str or ~numpy.ndarray
            If None, the channel weights will be equal and set to unity.

            If 'tsys' the channel weights will be calculated as:

             :math:`w = t_{exp} \times \delta\nu/T_{sys}^2`

            If 'spectral', the weights in each Spectrum will be used.

            If an array, it must have shape `(len(spectra)+1,)` or `(len(spectra)+1,nchan)` where `nchan` is the
            number of channels in the spectra.
            The first element ofthe weights array will be applied to the current spectrum.
        align : bool
            If `True` align the `spectra` to itself.
            This uses `Spectrum.align_to`.

        Returns
        -------
        average : `Spectrum`
            Averaged spectra.
        """

        if type(spectra) is not list:
            spectra = [spectra]

        spectra += [self]

        return average_spectra(spectra, weights=weights, align=align, history=self.history)

    def cog(
        self,
        vc=None,
        width_frac=None,
        bchan=None,
        echan=None,
        flat_tol=0.1,
        fw=1,
        xunit="km/s",
        vframe=None,
        doppler_convention=None,
    ) -> dict:
        """
        Curve of growth (CoG) analysis based on Yu et al. (2020) [1]_.

        Parameters
        ----------
        vc : float
            Central velocity of the line.
            If not provided, it will be estimated from the moment 1 of the `Spectrum`.
        width_frac : list
            List of fractions of the total flux at which to compute the line width.
            If 0.25 and 0.85 are not included, they will be added to estimate the concentration
            as defined in [1]_.
        bchan : int
            Beginning channel where there is signal.
            If not provided it will be estimated using `fw` times the width of the line at the largest `width_frac`.
        echan : int
            End channel where there is signal.
            If not provided it will be estimated using `fw` times the width of the line at the largest `width_frac`.
        flat_tol : float
            Tolerance used to define the flat portion of the curve of growth.
            The curve of growth will be considered flat when it's slope is within `flat_tol` times the standard deviation of the slope from zero.
        fw : float
            When estimating the line-free range, use `fw` times the largest width.
        xunit : str or `~astropy.units.quantity.Quantity`
            Units for the x axis when computing the CoG.
        vframe : None or str
            Velocity frame to use.
            The results will be provided in this velocity frame.
            If None, the velocity frame of the `Spectrum` will be used.
            The velocity frame of the `Spectrum` won't be changed.
        doppler_convention : None or {'radio', 'optical', 'relativistic'}
            The Doppler velocity covention to use when converting to velocity units.
            If a string, one of 'optical', 'radio', or 'relativistic'.
            If None, it will use the `Spectrum.doppler_convention`.

        Returns
        -------
        results : dict
            Dictionary with the flux (:math:`F`), width (:math:`V`), flux asymmetry (:math:`A_F`), slope asymmetry (:math:`A_C`), concentration (:math:`C_V`),
            rms, central velocity ("vel"), redshifted flux (:math:`F_r`), blueshifted flux (:math:`F_b`),`bchan` and `echan`, and errors on the
            flux ("flux_std"), width ("width_std"), central velocity ("vel_std"), redshifted flux ("flux_r_std"), and blueshifted flux ("flux_b_std").
            The rms is the standard deviation in the range outside of (bchan,echan).

        .. [1] `N. Yu, L. Ho & J. Wang, "On the Determination of Rotation Velocity and Dynamical Mass of Galaxies Based on Integrated H I Spectra"
           <https://ui.adsabs.harvard.edu/abs/2020ApJ...898..102Y/abstract>`_.
        """
        if width_frac is None:
            width_frac = [0.25, 0.65, 0.75, 0.85, 0.95]
        y = self.flux
        x = self.velocity_axis_to(unit=xunit, toframe=vframe, doppler_convention=doppler_convention)
        vframe = x.observer.name
        logger.info(f"Velocity frame: {frame_to_label[vframe]}")
        doppler_convention = x.doppler_convention
        logger.info(f"Doppler convention: {doppler_convention}")
        rdict = curve_of_growth(x, y, vc=vc, width_frac=width_frac, bchan=bchan, echan=echan, flat_tol=flat_tol, fw=fw)
        rdict.update({"vframe": vframe, "doppler_convention": doppler_convention})
        return rdict

    def _min_max_freq(self):
        """Return the sorted min and max frequency (in Hz) of the spectrum, regardless of the units of its axis"""
        start_freq = self.spectral_axis.quantity[0].to("Hz", equivalencies=u.spectral())
        end_freq = self.spectral_axis.quantity[-1].to("Hz", equivalencies=u.spectral())
        return Quantity(np.sort([start_freq.value, end_freq.value]), unit=start_freq.unit)

    @docstring_parameter(str(all_cats()))
    def query_lines(
        self, chemical_name: str | None = None, intensity_lower_limit: float | None = None, cat: str = "gbtlines"
    ) -> Table:
        """
        Query locally or remotely for lines and return a table object. The query returns lines
        with rest frequencies in the range of this Spectrum's spectral_axis.

        **Note:** If the search parameters result in no matches, a zero-length Table will be returned.

        Parameters
        ----------
        chemical_name : str, optional
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
        intensity_lower_limit : float, optional
                Lower limit on the intensity in the logarithmic CDMS/JPL scale.  This corresponds to the 'intintensity' column in the returned table.
        cat : str or Path
            The catalog to use.  One of: {0}  (minimum string match) or a valid Path to a local astropy-compatible table.  The local table
            must have all the columns listed in the `columns` parameter.

                - `'gbtlines'` is a local catalog of spectral lines between 300 MHz and 120 GHz with CDMS/JP log(intensity) > -9.

                - `'gbtrecomb'` is a local catalog of H, He, and C recombination lnes between 300 MHz and 120 GHz.

        Returns
        -------
        ~astropy.table.Table
            An astropy table containing the results of the search

        """
        minf, maxf = self._min_max_freq()
        return SpectralLineSearch.query_lines(
            min_frequency=minf,
            max_frequency=maxf,
            intensity_lower_limit=intensity_lower_limit,
            cat=cat,
            intensity_type="CDMS/JPL (log)",
        )

    @docstring_parameter(str(all_cats()))
    def recomb(self, line, cat: str = "gbtrecomb") -> Table:
        """
        Search for recombination lines of H, He, and C in the frequency range of this Spectrum.

        Parameters
        ----------
        line : str
           A string describing the line or series to search for, e.g. "Hydrogen", "Halpha", "Hebeta", "C", "carbon".

        cat : str or Path
            The catalog to use.  One of: {0}  (minimum string match) or a valid Path to a local astropy-compatible table.  The local table
            must have all the columns listed in the `columns` parameter.

                - `'gbtlines'` is a local catalog of spectral lines between 300 MHz and 120 GHz with CDMS/JP log(intensity) > -9.

                - `'gbtrecomb'` is a local catalog of H, He, and C recombination lnes between 300 MHz and 120 GHz.

        Returns
        -------
        ~astropy.table.Table
            An astropy table containing the results of the search

        """
        minf, maxf = self._min_max_freq()
        return SpectralLineSearch.recomb(
            min_frequency=minf,
            max_frequency=maxf,
            line=line,
            cat=cat,
        )

    @docstring_parameter(str(all_cats()))
    def recomball(self, cat: str = "gbtrecomb") -> Table:
        """
        Fetch all recombination lines of H, He, C in the frequency range of this Spectrum from the catalog.

        Parameters
        ----------
        cat : str or Path
            The catalog to use.  One of: {0}  (minimum string match) or a valid Path to a local astropy-compatible table.  The local table
            must have all the columns listed in the `columns` parameter.

                - `'gbtlines'` is a local catalog of spectral lines between 300 MHz and 120 GHz with CDMS/JP log(intensity) > -9.

                - `'gbtrecomb'` is a local catalog of H, He, and C recombination lnes between 300 MHz and 120 GHz.

        Returns
        -------
        ~astropy.table.Table
            An astropy table containing the results of the search

        """
        minf, maxf = self._min_max_freq()
        return SpectralLineSearch.recomball(min_frequency=minf, max_frequency=maxf, cat=cat)

    def meta_as_table(self):
        """
        Return `Spectrum.meta` as an `~astropy.table.Table`.
        """

        meta_array = {k: [v] for k, v in self.meta.items()}
        # TDIM7 gets lost somewhere. Add it back.
        if "TDIM7" not in meta_array.keys():
            meta_array["TDIM7"] = [f"({len(self.data), 1, 1, 1})"]
        d = {}
        d["HISTORY"] = self.history
        d["COMMENT"] = self.comments
        return Table(meta_array, meta=d)

    def _make_bintable(self, flags: bool) -> BinTableHDU:
        """
        Create a :class:`~astropy.io.fits.BinaryTableHDU` from the data of this `Spectrum`.
        """

        cd = BinTableHDU(data=self.meta_as_table(), name="SINGLE DISH").columns
        data_format = f"{np.shape(self.data)[0]}E"
        cd.add_col(Column(name="DATA", format=data_format, array=[self.data]))
        # Re-arrange so DATA is column 7.
        cd1 = cd[:6] + cd[-1] + cd[6:-1]
        if flags:
            flags = self.mask.astype(np.uint8)
            flag_format = f"{np.shape(flags)[0]}B"
            cd1.add_col(Column(name="FLAGS", format=flag_format, array=[flags]))
        return BinTableHDU.from_columns(cd1, name="SINGLE DISH")

    def _write_sdfits(
        self, fileobj, flags: bool = True, output_verify="exception", overwrite: bool = False, checksum: bool = False
    ) -> None:
        """
        Write this `Spectrum` as an SDFITS.

        Parameters
        ----------
        fileobj : str, file-like or `pathlib.Path`
            File to write to.  If a file object, must be opened in a
            writeable mode.
        flags : bool, optional
            Write the mask as flags.
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

        self._make_bintable(flags=flags).writeto(
            name=fileobj, output_verify=output_verify, overwrite=overwrite, checksum=checksum
        )


# @todo figure how how to document write()
####################################################################
# There is probably a less brute-force way to do this but I haven't
# been able to figure it out.  astropy.io.registry tools are not
# well explained.  register_writer 'consumes' the format keyword, so
# it cannot be passed along via a single overaarching write method,
# e.g., spectrum_writer()
####################################################################
def ascii_spectrum_writer_basic(spectrum, fileobj, **kwargs):
    spectrum._write_table(fileobj, format="ascii.basic", **kwargs)


def ascii_spectrum_writer_commented_header(spectrum, fileobj, **kwargs):
    spectrum._write_table(fileobj, format="ascii.commented_header", **kwargs)


def ascii_spectrum_writer_fixed_width(spectrum, fileobj, **kwargs):
    spectrum._write_table(fileobj, format="ascii.fixed_width", **kwargs)


def ascii_spectrum_writer_ipac(spectrum, fileobj, **kwargs):
    spectrum._write_table(fileobj, format="ascii.ipac", **kwargs)


def spectrum_writer_votable(spectrum, fileobj, **kwargs):
    spectrum._write_table(fileobj, format="votable", **kwargs)


def spectrum_writer_ecsv(spectrum, fileobj, **kwargs):
    spectrum._write_table(fileobj, format="ascii.ecsv", **kwargs)


def spectrum_writer_mrt(spectrum, fileobj, **kwargs):
    spectrum._write_table(fileobj, format="mrt", **kwargs)


def spectrum_writer_fits(spectrum, fileobj, **kwargs):
    spectrum._write_table(fileobj, format="fits", **kwargs)


def spectrum_writer_sdfits(spectrum, fileobj, **kwargs):
    spectrum._write_sdfits(fileobj, **kwargs)


def _read_table(fileobj, format, **kwargs):
    # GBTIDL format example:
    #    Scan:    152          NGC2415 2021-02-10 +07 38 37.5
    #                             Ta
    #     GHz-HEL                 YY
    #   1.4143356980054205     -0.1042543
    #   1.4143349827132639      0.0525000
    #
    # Note we don't get the full coordinates, just RA!
    if format == "gbtidl":
        # Read it into a pandas dataframe, which will have one column
        # with a name like 'Scan:    152          NGC2415 2021-02-10 +07 38 37.5'
        # We can split out the data array
        # We do not use e.g., Table.read(fileobs,data_start=3) because we need to
        # get the header info and don't want to open the file twice.
        # Note that functions like Table.read, pandas.read*, np.load_txt also natively
        # recognize compressed files, which is why we aren't using raw Python I/O
        # which does not.
        # t = Table.read(fileobj, format="ascii.fixed_width")  # ,data_start=3,header_start=3)
        df = pd.read_table(fileobj)
        df2 = df[df.columns[0]][2:].str.split(expand=True).astype(float).add_prefix("col")
        spectral_axis = df2["col0"]
        flux = df2["col1"]
        # Parse the first line of the header and put into meta
        _tmp, scan, target, date, ra = df.columns[0].split(maxsplit=4)
        meta = {}
        meta["SCAN"] = int(scan)
        meta["OBJECT"] = target
        meta["DATE-OBS"] = Time(date).to_string()
        c = SkyCoord(ra + " 0 0 0.0", unit=(u.hour, u.deg))
        meta["RA"] = c.ra.degree
        # Parse the 2nd link of the header to get the veldef and flux unit
        # Sometimes velocity_convention isn't there!
        h2 = df[df.columns[0]][0].split()
        if len(h2) == 2:
            velocity_convention = h2[0][0:4]
            flux_unit = h2[1]
        else:
            velocity_convention = "None"
            flux_unit = h2[0]
        if flux_unit == "Ta":
            fu = u.K
        elif flux_unit == "Counts":
            fu = u.ct
        else:
            fu = u.dimensionless_unscaled
        # parse the 3rd line column names, We only care about the ferquency axis
        h3 = df[df.columns[0]][1].split()
        units, vd = h3[0].split("-")
        spectral_axis = spectral_axis.values * u.Unit(units)  # noqa: PD011
        meta["VELDEF"] = velocity_convention + "-" + vd
        meta["POL"] = h3[1]

        s = Spectrum(flux=flux.values * fu, spectral_axis=spectral_axis, meta=meta)  # noqa: PD011
        return s

    t = Table.read(fileobj, format=format, **kwargs)
    f = t["flux"].value * t["flux"].unit
    return Spectrum.make_spectrum(f, meta=t.meta, observer_location="from_meta")


def spectrum_reader_ecsv(fileobj, **kwargs):
    return _read_table(fileobj, format="ascii.ecsv", **kwargs)


def spectrum_reader_fits(fileobj, **kwargs):
    return _read_table(fileobj, format="fits", **kwargs)


def spectrum_reader_gbtidl(fileobj, **kwargs):
    return _read_table(fileobj, format="gbtidl", **kwargs)


with registry.delay_doc_updates(Spectrum):
    # WRITERS
    registry.register_writer("ascii.basic", Spectrum, ascii_spectrum_writer_basic)
    registry.register_writer("basic", Spectrum, ascii_spectrum_writer_basic)
    registry.register_writer("ascii.commented_header", Spectrum, ascii_spectrum_writer_commented_header)
    registry.register_writer("commented_header", Spectrum, ascii_spectrum_writer_commented_header)
    registry.register_writer("ascii.fixed_width", Spectrum, ascii_spectrum_writer_fixed_width)
    registry.register_writer("fixed_width", Spectrum, ascii_spectrum_writer_fixed_width)
    registry.register_writer("ascii.ipac", Spectrum, ascii_spectrum_writer_ipac)
    registry.register_writer("ipac", Spectrum, ascii_spectrum_writer_ipac)
    registry.register_writer("votable", Spectrum, spectrum_writer_votable)
    registry.register_writer("ecsv", Spectrum, spectrum_writer_ecsv)
    # UnifiedOutputRegistry.write uses all caps if format not specified
    # and it believes the desired output is ecsv
    registry.register_writer("ECSV", Spectrum, spectrum_writer_ecsv)
    registry.register_writer("mrt", Spectrum, spectrum_writer_mrt)
    registry.register_writer("fits", Spectrum, spectrum_writer_fits)
    registry.register_writer("sdfits", Spectrum, spectrum_writer_sdfits)
    # UnifiedOutputRegistry.write uses retrurns tabular-fits if format
    # not specified and it believes the desired output is fits.
    registry.register_writer("tabular-fits", Spectrum, spectrum_writer_fits)

    # READERS
    registry.register_reader("fits", Spectrum, spectrum_reader_fits)
    registry.register_reader("gbtidl", Spectrum, spectrum_reader_gbtidl)
    registry.register_reader("ascii.ecsv", Spectrum, spectrum_reader_ecsv)
    registry.register_reader("ecsv", Spectrum, spectrum_reader_ecsv)

    # We aren't going to support these since they don't have easily digestible metadata
    # if they have metadata at all.
    # registry.register_writer("ascii.basic", Spectrum, ascii_spectrum_reader_basic)
    # registry.register_writer("basic", Spectrum, ascii_spectrum_reader_basic)
    # registry.register_writer("ascii.commented_header", Spectrum, ascii_spectrum_reader_commented_header)
    # registry.register_writer("commented_header", Spectrum, ascii_spectrum_reader_commented_header)
    # registry.register_writer("ascii.fixed_width", Spectrum, ascii_spectrum_reader_fixed_width)
    # registry.register_writer("fixed_width", Spectrum, ascii_spectrum_reader_fixed_width)
    # registry.register_writer("ascii.ipac", Spectrum, ascii_spectrum_reader_ipac)
    # registry.register_writer("ipac", Spectrum, ascii_spectrum_reader_ipac)
    # registry.register_writer("votable", Spectrum, spectrum_reader_votable)
    # registry.register_writer("mrt", Spectrum, spectrum_reader_mrt)


def average_spectra(spectra, weights="tsys", align=False, history=None):
    r"""
    Average `spectra`. The resulting `average` will have an exposure equal to the sum of the exposures,
    and coordinates and system temperature equal to the weighted average of the coordinates and system temperatures.

    Parameters
    ----------
    spectra : list of `Spectrum`
        Spectra to be averaged. They must have the same number of channels.
        No checks are done to ensure they are aligned.
    weights: None, str or ~numpy.ndarray
        If None, the channel weights will be equal and set to unity.

        If 'tsys' the channel weights will be calculated as:

         :math:`w = t_{exp} \times \delta\nu/T_{sys}^2`

        If 'spectral', the `weights` array in each Spectrum will be used.

        If an array, it must have shape `(len(spectra),)` or `(len(spectra),nchan)` where `nchan` is the number of channels in the spectra.
    align : bool
        If `True` align the `spectra` to the first element.
        This uses `Spectrum.align_to`.
    history : `dysh.log.HistoricalBase`
        History to append to the averaged spectra.

    Returns
    -------
    average : `Spectrum`
        Averaged spectra.
    """

    nspec = len(spectra)
    nchan = len(spectra[0].data)
    shape = (nspec, nchan)
    _data = np.empty(shape, dtype=float)
    _mask = np.zeros(shape, dtype=bool)
    data_array = np.ma.MaskedArray(_data, mask=_mask, dtype=float, fill_value=np.nan)
    wts = np.empty(shape, dtype=float)
    exposures = np.empty(nspec, dtype=float)
    durations = np.empty(nspec, dtype=float)
    tsyss = np.empty(nspec, dtype=float)
    xcoos = np.empty(nspec, dtype=float)
    ycoos = np.empty(nspec, dtype=float)
    ap_eff = np.empty(nspec, dtype=float)
    surface_error = np.empty(nspec, dtype=float)
    zenith_opacity = np.empty(nspec, dtype=float)
    observer = spectra[0].observer
    units = spectra[0].flux.unit
    seunit = spectra[0].meta.get("SE_UNIT", "")
    pols = []
    for i, s in enumerate(spectra):
        if not isinstance(s, Spectrum):
            raise ValueError(f"Element {i} of `spectra` is not a `Spectrum`. {type(s)}")
        if units != s.flux.unit:
            raise ValueError(
                f"Element {i} of `spectra` has units {s.flux.unit}, but the first element has units {units}."
            )
        if align:
            if i > 0:
                s = s.align_to(spectra[0])
        data_array[i] = s.data
        data_array[i].mask = s.mask

        # Assign a single weight to all channels in the Spectrum s
        if isinstance(weights, np.ndarray):
            wts[i] = weights[i]
        elif weights == "tsys":
            wts[i] = core.tsys_weight(s.meta["EXPOSURE"], s.meta["CDELT1"], s.meta["TSYS"])
        elif weights == "spectral":
            wts[i] = s.weights
        else:
            wts[i] = 1.0

        exposures[i] = s.meta["EXPOSURE"]
        durations[i] = s.meta["DURATION"]
        tsyss[i] = s.meta["TSYS"]
        xcoos[i] = s.meta["CRVAL2"]
        ycoos[i] = s.meta["CRVAL3"]
        ap_eff[i] = s.meta["AP_EFF"]
        surface_error[i] = s.meta["SURF_ERR"]
        # if data are in Ta units, then there  normally wouldn't be a zenith opacity provided
        zenith_opacity[i] = s.meta.get("TAU_Z", -1)
        pols.append(s.meta["CRVAL4"])
    _mask = np.isnan(data_array.data) | data_array.mask
    wts = np.ma.masked_invalid(wts)
    data_array = np.ma.MaskedArray(data_array, mask=_mask, fill_value=np.nan)
    data, sum_of_weights = np.ma.average(data_array, axis=0, weights=wts, returned=True)
    tsys = np.ma.average(tsyss, axis=0, weights=wts[:, 0])
    xcoo = np.ma.average(xcoos, axis=0, weights=wts[:, 0])
    ycoo = np.ma.average(ycoos, axis=0, weights=wts[:, 0])
    ap = np.ma.average(ap_eff, axis=0, weights=wts[:, 0])
    se = np.ma.average(surface_error, axis=0, weights=wts[:, 0])
    zenith_opacity = np.ma.masked_where(zenith_opacity < 0, zenith_opacity)
    ze = np.ma.average(zenith_opacity, axis=0, weights=wts[:, 0])
    exposure = exposures[~wts.mask[:, 0]].sum(axis=0)
    duration = durations[~wts.mask[:, 0]].sum(axis=0)

    new_meta = deepcopy(spectra[0].meta)
    new_meta["TSYS"] = tsys
    new_meta["EXPOSURE"] = exposure
    new_meta["DURATION"] = duration
    new_meta["CRVAL2"] = xcoo
    new_meta["CRVAL3"] = ycoo
    new_meta["AP_EFF"] = ap
    new_meta["SURF_ERR"] = se
    new_meta["SE_UNIT"] = seunit
    if not hasattr(ze, "mask"):
        new_meta["TAU_Z"] = ze

    upols = set(pols)  # unique crval4's being averaged (polarizations)
    numpols = len(upols)
    if numpols == 1:  # only one pol being averaged, no change needed
        new_meta["CRVAL4"] = next(iter(upols))
    elif numpols == 2:  # two pols being averaged, check that it is XX and YY or LL and RR, invalid otherwise
        if (upols == {-5, -6}) or (upols == {-1, -2}):
            new_meta["CRVAL4"] = 1
        else:
            new_meta["CRVAL4"] = 0
    elif numpols >= 3:  # 3 or more pols, invalid
        new_meta["CRVAL4"] = 0

    averaged = Spectrum.make_spectrum(Masked(data * units, data.mask), meta=new_meta, observer=observer)
    averaged._weights = sum_of_weights
    if history is not None:
        # Keep previous history first.
        averaged._history = history + averaged._history

    return averaged
