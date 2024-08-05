"""
The Spectrum class to contain and manipulate spectra.
"""

import warnings
from copy import deepcopy

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord, SpectralCoord, StokesCoord
from astropy.coordinates.spectral_coordinate import NoVelocityWarning
from astropy.io import registry
from astropy.io.fits.verify import VerifyWarning
from astropy.modeling.fitting import LinearLSQFitter
from astropy.nddata.ccddata import fits_ccddata_writer
from astropy.table import Table
from astropy.time import Time
from astropy.wcs import WCS, FITSFixedWarning
from ndcube import NDCube
from specutils import Spectrum1D

from dysh.spectra import core

from ..coordinates import (  # is_topocentric,; topocentric_velocity_to_frame,
    KMS,
    Observatory,
    astropy_frame_dict,
    change_ctype,
    get_velocity_in_frame,
    make_target,
    replace_convention,
    sanitize_skycoord,
    veldef_to_convention,
)

# from ..plot.dev.iPlotter import SpectrumPlot as sp
from ..util import minimum_string_match
from . import baseline, get_spectral_equivalency

# from astropy.nddata import StdDevUncertainty


class Spectrum(Spectrum1D):
    """
    This class contains a spectrum and its attributes. It is built on
    `~specutils.Spectrum1D` with added attributes like baseline model.
    Note that `~specutils.Spectrum1D` can contain multiple spectra but
    we probably will not use that because the restriction that it can
    have only one spectral axis conflicts with slight Doppler shifts.
    See `~specutils.Spectrum1D` for the instantiation arguments.
    """

    def __init__(self, *args, **kwargs):
        # print(f"ARGS={args}")
        self._target = kwargs.pop("target", None)
        if self._target is not None:
            # print(f"self._target is {self._target}")
            self._target = sanitize_skycoord(self._target)
            self._velocity_frame = self._target.frame.name
        else:
            self._velocity_frame = None
        # @todo - have _observer_location attribute instead?
        # and observer property returns getITRS(observer_location,obstime)
        self._observer_location = kwargs.pop("observer_location", None)
        self._observer = kwargs.pop("observer", None)
        if self._observer is not None and self._observer_location is not None:
            raise Exception("You can only specify one of observer_location or observer")
        Spectrum1D.__init__(self, *args, **kwargs)
        # super(Spectrum1D, self).__init__(*args, **kwargs)
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
        # if "uncertainty" not in kwargs:
        # self._uncertainty = StdDevUncertainty(np.ones_like(self.data), unit=self.flux.unit)
        # Weights are Non here because can be defined
        # separately from uncertainty by the user. Not really recommended
        # but they can do so. See self.weights()
        # self._weights = None
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
        # @todo fix resolution, including what units to use
        self._resolution = 1  # placeholder

    @property
    def weights(self):
        """The channel weights of this spectrum"""
        # Cannot take power of NDUncertainty, so convert
        # to a Quanity first.
        # if self._weights is None:
        #    return self._uncertainty.quantity**-2
        return self._weights

    @property
    def exclude_regions(self):
        """The baseline exclusion region(s) of this spectrum"""
        return self._exclude_regions

    def _toggle_sections(self, nchan, s):
        """helper routine to toggle between an include= and exclude=
        only works in channel (0..nchan-1) units
        sections s need to be a list of (start_chan,end_chan) tuples,
        for example [(100,200),(500,600)] would be an include=
        An exclude= needs to start with 0
        channels need to be ordered low to high, but there is no check
        for this yet!
        """
        ns = len(s)
        s1 = []
        e = 0  #  set this to 1 if you want to be exact complementary
        if s[0][0] == 0:
            # print("toggle_sections: edged")
            for i in range(ns - 1):
                s1.append((s[i][1] + e, s[i + 1][0] - e))
        else:
            # print("toggle_sections: internal")
            s1.append((0, s[0][0]))
            for i in range(ns - 1):
                s1.append((s[i][1], s[i + 1][0]))
            s1.append((s[ns - 1][1], nchan - 1))
        return s1

    ##@todo
    # def exclude_region(self,region):
    # where region is SpectralRegion, channels, velocity, etc.  See core.py baseline method.
    #
    # def region_to_mask():
    #  set spectrum mask to True inside exclude_regions. normally we don't do this for baselining

    @property
    def baseline_model(self):
        """Returns the computed baseline model or None if it has not yet been computed."""
        return self._baseline_model

    def baseline(self, degree, exclude=None, include=None, **kwargs):
        # fmt: off
        """
        Compute and optionally remove a baseline.  The model for the
        baseline can be either a
        `1D polynomial model <https://docs.astropy.org/en/latest/api/astropy.modeling.polynomial.Polynomial1D.html>`_ or a
        `1D Chebyshev polynomial of the first kind <https://docs.astropy.org/en/latest/api/astropy.modeling.polynomial.Chebyshev1D.html>`_.
        The code uses `astropy.modeling`
        and `astropy.fitter` to compute the baseline.  See the documentation for those modules.
        This method will set the `baseline_model` attribute to the fitted model function which can be evaluated over a domain.

        Note that include= and exclude= are mutually exclusive.

        Parameters
        ----------
            degree : int
                The degree of the polynomial series, a.k.a. baseline order
            exclude : list of 2-tuples of int or ~astropy.units.Quantity, or ~specutils.SpectralRegion
                List of region(s) to exclude from the fit.  The tuple(s) represent a range in the form [lower,upper], inclusive.
                In channel units.

                Examples:

                One channel-based region: [11,51]

                Two channel-based regions: [(11,51),(99,123)].

                One ~astropy.units.Quantity region: [110.198*u.GHz,110.204*u.GHz].

                One compound `~specutils.SpectralRegion`: SpectralRegion([(110.198*u.GHz,110.204*u.GHz),(110.196*u.GHz,110.197*u.GHz)]).

                Default: no exclude region

            include: list of 2-tuples of int (currently units not supported yet, pending issue 251/260)

            model : str
                One of 'polynomial' 'chebyshev', 'legendre', or 'hermite'
                Default: 'chebyshev'
            fitter  :  `~astropy.fitting._FitterMeta`
                The fitter to use. Default: `~astropy.fitter.LinearLSQFitter` (with `calc_uncertaintes=True`).
                Be care when choosing a different fitter to be sure it is optimized for this problem.
            remove : bool
                If True, the baseline is removed from the spectrum. Default: False
            normalize : bool
                If True, the frequency axis is internally rescaled from 0..1
                to avoid roundoff problems (and make the coefficients slightly more
                understandable). This is usually needed for a polynomial, though overkill
                for the others who do their own normalization.
                CAVEAT:   with normalize=True, you cannot undo a baseline fit.
                Default: False

        """
        # fmt: on
        # @todo: Are exclusion regions OR'd with the existing mask? make that an option?
        kwargs_opts = {
            "remove": False,
            "normalize": False,
            "model": "chebyshev",
            "fitter": LinearLSQFitter(calc_uncertainties=True),
        }
        kwargs_opts.update(kwargs)

        if kwargs_opts["normalize"]:
            print("Warning: baseline fit done in [0,1) space, even though it might say Hz (issue ###)")
            spectral_axis = deepcopy(self._spectral_axis)  # save the old axis
            self._normalized = True  # remember it's now normalized
            nchan = len(spectral_axis)
            for i in range(nchan):
                self._spectral_axis[i] = (i * 1.0 / nchan) * u.Hz  # would like to use "u.chan" units - not working yet
            # some @todo here about single setter, units u.chan etc.

        # include= and exclude= are mutually exclusive, but we allow include=
        # if include is used, transform it to exclude=
        if include != None:
            if exclude != None:
                print(f"Warning: ignoring exclude={exclude}")
            nchan = len(self._spectral_axis)
            exclude = self._toggle_sections(nchan, include)

        self._baseline_model = baseline(self, degree, exclude, **kwargs)

        if kwargs_opts["remove"]:
            s = self.subtract(self._baseline_model(self.spectral_axis))
            self._data = s._data
            self._subtracted = True

        if kwargs_opts["normalize"]:
            self._spectral_axis = spectral_axis
            del spectral_axis

    # baseline

    def undo_baseline(self):
        """
        Undo the most recently computed baseline. If the baseline
        has been subtracted, it will be added back. The `baseline_model`
        attribute is set to None. Exclude regions are untouched.
        """
        if self._baseline_model is None:
            return
        if self._subtracted:
            if self._normalized:
                warnings.warn("Cannot undo previously normalized baseline subtraction")
                return
            s = self.add(self._baseline_model(self.spectral_axis))
            self._data = s._data
            self._baseline_model = None

    def _set_exclude_regions(self, exclude):
        """
        Set the mask for the regions to exclude.

        Parameters
        ----------
            exclude : list of 2-tuples of int or ~astropy.units.Quantity, or ~specutils.SpectralRegion
                List of region(s) to exclude from the fit.  The tuple(s) represent a range in the form [lower,upper], inclusive.
                In channel units.

                Examples: One channel-based region: [11,51],
                          Two channel-based regions: [(11,51),(99,123)].
                          One ~astropy.units.Quantity region: [110.198*u.GHz,110.204*u.GHz].
                          One compound ~specutils.SpectralRegion: SpectralRegion([(110.198*u.GHz,110.204*u.GHz),(110.196*u.GHz,110.197*u.GHz)]).

        """
        pass

    def list_to_spectral_region(self, inlist):
        # @todo utility code to convert a input list of channels or quantities to a spectral region with units of self.spectral_axis.unit.
        # This could go in core.py combine this with _set_exclude_regions
        pass

    def bshow(self):
        """Show the baseline model"""
        print(f"baseline model {self._baseline_model}")

    def plot(self, **kwargs):
        self._plotter = sp(self, **kwargs)

    @property
    def obstime(self):
        return self._obstime

    @property
    def plotter(self):
        return self._plotter

    def stats(self, roll=0):
        """
        Compute some statistics of this `Spectrum`.  The mean, rms,
        data minimum and data maximum are calculated.  Note this works
        with slicing, so, e.g.,  `myspectrum[45:153].stats()` will return
        the statistics of the slice.

        Parameters
        ----------
            roll : int
                Return statistics on a 'rolled' array differenced with the
                origibnal array. If no correllaton between subsequent point,
                a roll=1 would return an RMS  sqrt(2) larger than that of the
                input array. Another advantage of rolled statistics it will
                remove most slow variations, thus RMS/sqrt(2) might be a better
                indication of the underlying RMS.
        Returns
        -------
        stats : tuple
            Tuple consisting of (mean,rms,datamin,datamax)


        """
        # @todo: maybe make this a dict return value a dict
        if roll == 0:
            mean = self.mean()
            rms = self.data.std()
            dmin = self.min()
            dmax = self.max()
        else:
            d = self._data[roll:] - self._data[:-roll]
            mean = d.mean()
            rms = d.std()
            dmin = d.min()
            dmax = d.max()
        return (mean, rms, dmin, dmax)

    def _decimate(self, n, offset=0):
        """Decimate a spectrum by n pixels, starting at pixel offset
        @todo deprecate?   decimation is in smooth, do we need a separate one?
        """
        nchan = len(self._data)
        print("_decimate deprecated", nchan, offset, n)
        idx = np.arange(offset, nchan, n)
        new_data = self._data[idx] * u.K  # why units again?
        s = Spectrum.make_spectrum(new_data, meta=self.meta)
        s._spectral_axis = self._spectral_axis[idx]
        # @todo  fix WCS
        return s

    def smooth(self, method="hanning", width=1, decimate=0, kernel=None):
        """
        Smooth or Convolve a spectrum, optionally decimating it.

        Default smoothing is hanning.

        Parameters
        ----------
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
            Decimation factor of the spectrum by returning every width channel.
            -1:   no decimation
            0:    use the width parameter
            >1:   user supplied decimation (use with caution)
            The default is 0, meaning decimation is by 'width'
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
        s : Spectrum
            The new, possibly decimated, convolved spectrum.
            The meta data are currently passed on,and hence will
            contain some original WCS parameters.
        """
        nchan = len(self._data)
        decimate = int(decimate)

        # @todo  see also core.smooth() for valid_methods
        valid_methods = ["hanning", "boxcar", "gaussian"]
        this_method = minimum_string_match(method, valid_methods)
        if this_method == None:
            raise Exception(f"smooth({method}): valid methods are {valid_methods}")

        if this_method == "gaussian":
            stddev = np.sqrt(width**2 - self._resolution**2) / 2.35482
            s1 = core.smooth(self._data, this_method, stddev)
        else:
            s1 = core.smooth(self._data, this_method, width)

        new_data = s1 * self.flux.unit
        if decimate >= 0:
            if decimate == 0:
                # take the default decimation by 'width'
                idx = np.arange(0, nchan, width)
                new_resolution = 1  # new resolution in the new pixel width
                cell_shift = 0.5 * (width - 1) * (self._spectral_axis.value[1] - self._spectral_axis.value[0])
            else:
                # user selected decimation (but should be <= width)
                if decimate > width:
                    raise Exception(f"Cannot decimate with more than width {width}")
                idx = np.arange(0, nchan, decimate)
                new_resolution = np.sqrt(width**2 - decimate**2)  # @todo dangerous?
                cell_shift = 0.5 * (decimate - 1) * (self._spectral_axis.value[1] - self._spectral_axis.value[0])
            # print("    cell_shift:", cell_shift)
            new_data = new_data[idx]
            new_meta = deepcopy(self.meta)
            new_meta["CDELT1"] = width * self.meta["CDELT1"]  # @todo etc ???
            s = Spectrum.make_spectrum(new_data, meta=new_meta)
            s._spectral_axis = self._spectral_axis[idx]
            for i in range(len(s._spectral_axis)):  # grmpf, no proper setter
                s._spectral_axis.value[i] += cell_shift
            if self._baseline_model is not None:
                print("Warning: removing baseline_model")
                s._baseline_model = None  # was already None
            s._resolution = new_resolution
            # @todo  fix WCS
        else:
            s = Spectrum.make_spectrum(new_data, meta=self.meta)
            s._baseline_model = self._baseline_model  # it never got copied
            s._resolution = width
            # @todo   resolution is not well defined multiple methods are used in succession
        return s

    @property
    def equivalencies(self):
        """Get the spectral axis equivalencies that can be used in converting the axis
        between km/s and frequency or wavelength"""
        equiv = u.spectral()
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
        return self._target

    @property
    def observer(self):
        """
        Returns
        -------
            observer : `~astropy.coordinates.BaseCoordinateFrame` or derivative
            The coordinate frame of the observer if present.
        """
        if self._observer is None and self._observer_location is not None:
            return SpectralCoord._validate_coordinate(self._observer_location.get_itrs(obstime=self._obstime))
        else:
            return self._observer

    @property
    def velocity_frame(self):
        """String representation of the velocity frame"""
        return self._velocity_frame

    @property
    def doppler_convention(self):
        """String representation of the velocity (Doppler) convention"""
        return self.velocity_convention

    def axis_velocity(self, unit=KMS):
        """Get the spectral axis in velocity units.
        *Note*: This is not the same as `Spectrum.velocity`, which includes the source radial velocity.

        Parameters
        ----------
        unit : `~astropy.units.Quantity` or str that can be converted to Quantity
                The unit to which the axis is to be converted
        Returns
        -------
        velocity : `~astropy.units.Quantity`
                The converted spectral axis velocity
        """
        return self._spectral_axis.to(unit)

    def velocity_axis_to(self, unit=KMS, toframe=None, doppler_convention=None):
        """
        Parameters
        ----------
        unit : `~astropy.units.Quantity` or str that can be converted to Quantity
            The unit to which the axis is to be converted

        toframe : str
            The coordinate frame to convert to, e.g. 'hcrs', 'icrs'

        doppler_convention : str
            The Doppler velocity covention to use, one of 'optical', 'radio', or 'rest'

        Returns
        -------
        test_spectrum.pyvelocity : `~astropy.units.Quantity`
            The converted spectral axis velocity
        """
        if toframe is not None and toframe != self.velocity_frame:
            self.set_frame(toframe)
        if doppler_convention is not None:
            return self._spectral_axis.to(unit=unit, doppler_convention=doppler_convention).to(unit)
        else:
            return self.axis_velocity(unit)

    def get_velocity_shift_to(self, toframe):
        if self._target is None:
            raise Exception("Can't calculate velocity because Spectrum.target is None")
        return get_velocity_in_frame(self._target, toframe, self._observer, self._obstime)

    def set_frame(self, toframe):
        # @todo VELDEF should be changed as well?
        """Set the sky coordinate and doppler tracking reference frame of this Spectrum. The header 'CTYPE1' will be changed accordingly.

        To make a copy of this Spectrum with new coordinate referece frmae instead, use `with_frame`.

        Parameters
        ----------
        toframe - str
            The coordinate reference frame identifying string, as used by astropy, e.g. 'hcrs', 'icrs', etc.
        """
        if "topo" in toframe:
            actualframe = self.observer
        else:
            actualframe = astropy_frame_dict.get(toframe, toframe)
        # print(f"actual frame is {actualframe} {type(actualframe)}")
        self._spectral_axis = self._spectral_axis.with_observer_stationary_relative_to(actualframe)
        self._meta["CTYPE1"] = change_ctype(self._meta["CTYPE1"], toframe)
        if isinstance(actualframe, str):
            self._velocity_frame = actualframe
        else:
            self._velocity_frame = actualframe.name

    def with_frame(self, toframe):
        """Return a copy of this Spectrum with a new coordinate reference frame.

        Parameters
        ----------
        toframe - str
            The coordinate reference frame identifying string, as used by astropy, e.g. 'hcrs', 'icrs', etc.

        Returns
        -------
        spectrum : `~dysh.spectra.Spectrum`
            A new Spectrum object
        """

        s = self._copy()
        s.set_frame(toframe)
        return s

    def set_convention(self, doppler_convention):
        """Set the velocity convention of this Spectrum.  The spectral axis of this Spectrum will be replaced
        with a new spectral axis with the input velocity convention.  The header 'VELDEF' value will
        be changed accordingly.

        To make a copy of this Spectrum with a new velocity convention instead, use `with_velocity_convention`.

        Parameters
        ----------
        doppler_convention - str
            The velocity convention, one of 'radio', 'optical', 'relativistic'

        """
        # replicate() gives the same asnwer as
        # self._spectral_axis.to(unit=self._spectral_axis.unit, doppler_convention=doppler_convention)
        new_sp_axis = self.spectral_axis.replicate(doppler_convention=doppler_convention)
        self._spectral_axis = new_sp_axis
        self.meta["VELDEF"] = replace_convention(self.meta["VELDEF"], doppler_convention)

    def with_velocity_convention(self, doppler_convention):
        """Returns a copy of this Spectrum with the input velocity convention.  The header 'VELDEF' value will
        be changed accordingly.

        Parameters
        ----------
        doppler_convention - str
            The velocity convention, one of 'radio', 'optical', 'relativistic'

        Returns
        -------
        spectrum : `~dysh.spectra.Spectrum`
            A new Spectrum object
        """
        if False:
            # this doesn't work.
            # for some reason, the axis velocity
            # still contains the difference between TOPO and observed frame
            s = self.__class__(
                flux=self.flux,
                wcs=self.wcs,
                meta=self.meta,
                velocity_convention=doppler_convention,
                target=self._target,
                observer=self._spectral_axis.observer,
            )
            s.meta["VELDEF"] = replace_convention(self.meta["VELDEF"], doppler_convention)
        s = self._copy(velocity_convention=doppler_convention)
        s.set_convention(doppler_convention)
        return s

    def savefig(self, file, **kwargs):
        """Save the plot"""
        if self._plotter is None:
            raise Exception("You have to invoke plot() first")
        self._plotter.figure.savefig(file, **kwargs)

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
        description = ["Spectral axis", "Flux", udesc, "Channel weights", "Mask 0=unmasked, 1=masked", bldesc]
        # remove FITS reserve keywords
        meta = deepcopy(self.meta)
        meta.pop("NAXIS1")
        meta.pop("TDIM7")
        meta.pop("TUNIT7")
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
        # There are some that are more than 8 chars that should be fixed in GBTFISLOAD
        warnings.simplefilter("ignore", VerifyWarning)
        t.write(fileobj, format=format, **kwargs)

    def _copy(self, **kwargs):
        """
        Perform deep copy operations on each attribute of the ``Spectrum``
        object.
        This overrides the ``specutils.Spectrum1D`` method so that
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
    def fake_spectrum(cls, nchan=1024, **kwargs):
        """
        Create a fake spectrum, useful for simple testing. A default header is
        created, which may be modified with kwargs.

        Parameters
        ----------
        nchan : int, optional
            Number of channels. The default is 1024.

        **kwargs: dict or key=value
            Metadata to put in the header.  If the key exists already in
            the default header, it will be replaced. Otherwise the key and value will be
            added to the header. Keys are case insensitive.

        Returns
        -------
        spectrum : `~dysh.spectra.Spectrum`
            The spectrum object
        """
        data = np.random.rand(nchan) * u.K
        meta = {
            "OBJECT": "NGC2415",
            "BANDWID": 23437500.0,
            "DATE-OBS": "2021-02-10T07:38:37.50",
            "DURATION": 0.9982445,
            "EXPOSURE": 732.1785161896237,
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
            "OBSERVER": "Michael Fanelli",
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
            "RESTFREQ": 1420405751.7,
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
            "DOPFREQ": 1420405751.7,
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
            "RESTFRQ": 1420405751.7,
            "MEANTSYS": 17.16746070048293,
            "WTTSYS": 17.16574907094451,
        }
        for k, v in kwargs.items():
            meta[k.upper()] = v
        return Spectrum.make_spectrum(data, meta, observer_location=Observatory["GBT"])

    # @todo allow observer or observer_location.  And/or sort this out in the constructor.
    @classmethod
    def make_spectrum(cls, data, meta, use_wcs=True, observer_location=None):
        # , shift_topo=False):
        """Factory method to create a Spectrum object from a data and header.

        Parameters
        ----------
        data :  `~numpy.ndarray`
            The data array. See `~specutils.Spectrum1D`
        meta : dict
            The metadata, typically derived from an SDFITS header.
            Required items in `meta` are 'CTYPE[123]','CRVAL[123]', 'CUNIT[123]', 'VELOCITY', 'EQUINOX', 'RADESYS'
        use_wcs : bool
            If True, create a WCS object from `meta`
        observer_location : `~astropy.coordinates.EarthLocation` or str
            Location of the observatory. See `~dysh.coordinates.Observatory`.
            This will be transformed to `~astropy.coordinates.ITRS` using the time of observation DATE-OBS or MJD-OBS in `meta`.
            If this parameter is given the special str value 'from_meta', then an observer_location
            will be created from SITELONG, SITELAT, and SITEELEV in the meta dictionary.

        Returns
        -------
        spectrum : `~dysh.spectra.Spectrum`
            The spectrum object
        """
        # @todo add resolution being the channel separation, unless we use the FREQRES column
        warnings.simplefilter("ignore", NoVelocityWarning)
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

        # for k in _required:
        #    print(f"{k} {k in header}")
        if not _required <= meta.keys():
            raise ValueError(f"Header (meta) is missing one or more required keywords: {_required}")

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
                wcs = WCS(header=meta)
                # It would probably be safer to add NAXISi to meta.
                wcs.array_shape = (0, 0, 0, len(data))
                # For some reason these aren't identified while creating the WCS object.
                wcs.wcs.obsgeo[:3] = meta["SITELONG"], meta["SITELAT"], meta["SITEELEV"]
                # Reset warnings.
        else:
            wcs = None
        # is_topo = is_topocentric(meta["CTYPE1"])  # GBT-specific to use CTYPE1 instead of VELDEF
        target = make_target(meta)
        vc = veldef_to_convention(meta["VELDEF"])
        # vf = astropy_frame_dict[decode_veldef(meta["VELDEF"])[1]]
        # could be clever:
        # obstime = Time(meta.get("DATE-OBS",meta.get("MJD-OBS",None)))
        # if obstime is not None: blah blah
        if "DATE-OBS" in meta:
            obstime = Time(meta["DATE-OBS"])
        elif "MJD-OBS" in meta:
            obstime = Time(meta["MJD-OBS"])
        if "DATE-OBS" in meta or "MJD-OBS" in meta:
            if observer_location == "from_meta":
                try:
                    observer_location = Observatory.get_earth_location(
                        meta["SITELONG"], meta["SITELAT"], meta["SITEELEV"]
                    )
                except KeyError as ke:
                    raise Exception(f"Not enough info to create observer_location: {ke}")
            if observer_location is None:
                obsitrs = None
            else:
                obsitrs = SpectralCoord._validate_coordinate(observer_location.get_itrs(obstime=obstime))
        else:
            warnings.warn(
                "'meta' does not contain DATE-OBS or MJD-OBS. Spectrum won't be convertible to certain coordinate"
                " frames"
            )
            obsitrs = None
        s = cls(
            flux=data,
            wcs=wcs,
            meta=meta,
            velocity_convention=vc,
            radial_velocity=target.radial_velocity,
            rest_value=meta["RESTFRQ"] * u.Hz,
            observer=obsitrs,
            # observer_location=observer_location,
            target=target,
        )
        # For some reason, Spectrum1D.spectral_axis created with WCS do not inherit
        # the radial velocity. In fact, they get no radial_velocity attribute at all!
        # This method creates a new spectral_axis with the given radial velocity.
        if observer_location is None:
            s.set_radial_velocity_to(target.radial_velocity)  # open
        return s

    def _arithmetic_apply(self, other, op, handle_meta, **kwargs):
        if isinstance(other, NDCube):
            result = op(other, **{"handle_meta": handle_meta})
        else:
            result = op(other, **{"handle_meta": handle_meta, "meta_other_meta": False})
        self._shallow_copy_attributes(result)
        return result

    def _shallow_copy_attributes(self, other):
        other._target = self._target
        other._observer = self._observer
        other._velocity_frame = self._velocity_frame
        other._obstime = self._obstime
        other._baseline_model = self._baseline_model
        other._exclude_regions = self._exclude_regions
        other._mask = self._mask
        other._subtracted = self._subtracted
        other.spectral_axis._doppler_convention = self.doppler_convention

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

    def _add_meta(self, operand, operand2, **kwargs):
        # print(kwargs)
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
            with u.set_enabled_equivalencies(u.spectral()):
                wav_sp = spectral_axis.to(unit=wav.unit).replicate(value=wav.value, unit=wav.unit)
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

        # Slicing uses NumPY ordering by default.
        sliced_wcs = wcs[0:1, 0:1, 0:1, start_idx:stop_idx]

        # Update meta.
        meta = self.meta.copy()
        head = sliced_wcs.to_header()
        for k in ["CRPIX1", "CRVAL1"]:
            meta[k] = head[k]

        # New Spectrum.
        return self.make_spectrum(
            self.flux[start_idx:stop_idx], meta=meta, observer_location=Observatory[meta["TELESCOP"]]
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
        tmp, scan, target, date, ra = df.columns[0].split(maxsplit=4)
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
        spectral_axis = spectral_axis.values * u.Unit(units)
        meta["VELDEF"] = velocity_convention + "-" + vd
        meta["POL"] = h3[1]
        s = Spectrum(flux=flux.values * fu, spectral_axis=spectral_axis, meta=meta)
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
    registry.register_writer("mrt", Spectrum, spectrum_writer_mrt)
    registry.register_writer("fits", Spectrum, spectrum_writer_fits)
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
