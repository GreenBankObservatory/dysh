"""
The Spectrum class to contain and manipulate spectra.
"""

import warnings
from copy import deepcopy

import astropy.units as u
import numpy as np
from astropy.coordinates import SpectralCoord
from astropy.coordinates.spectral_coordinate import NoVelocityWarning
from astropy.io import registry
from astropy.io.fits.verify import VerifyWarning
from astropy.modeling.fitting import LinearLSQFitter
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
    decode_veldef,
    frame_to_label,
    get_velocity_in_frame,
    make_target,
    replace_convention,
    sanitize_skycoord,
    veldef_to_convention,
)
from ..plot import specplot as sp
from ..util import minimum_string_match
from . import baseline, get_spectral_equivalency


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
            self._target.sanitized = True
            self._velocity_frame = self._target.frame.name
        else:
            self._velocity_frame = None
        # @todo - have _observer_location attribute instead?
        # and observer property returns getITRS(observer_location,obstime)
        self._observer = kwargs.pop("observer", None)
        Spectrum1D.__init__(self, *args, **kwargs)
        self._spectral_axis._target = self._target
        self._spectral_axis._observer = self._observer
        if self._observer is not None:
            self._velocity_frame = self._observer.name
        if "DATE-OBS" in self.meta:
            self._obstime = Time(self.meta["DATE-OBS"])
        else:
            self._obstime = None
        # if "CTYPE1" in self.meta:
        #    # may not need these attributes anymore.
        #    if self._target is not None:
        #        self._target.topocentric = is_topocentric(self.meta["CTYPE1"])
        #    self.topocentric = is_topocentric(self.meta["CTYPE1"])

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
    def exclude_regions(self):
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
        if self._plotter is None:
            self._plotter = sp.SpectrumPlot(self, **kwargs)
        self._plotter.plot(**kwargs)

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
        @todo deprecate?   decimation is in smooth
        """
        nchan = len(self._data)
        print("PJT decimate", nchan, offset, n)
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



         smooth a spectrum
            method:    hanning, boxcar, gaussian, fft (not implemented)
            width:     in pixels
            decimate:  -1  none
                        0  use the width parameter
                       >0  use the decimate factor explicitly
            kernel:    give your own array to convolve with (not implemented)
        """
        nchan = len(self._data)
        decimate = int(decimate)
        # print("PJT smooth",method,nchan,width,decimate)
        # print("    old resolution: ",self._resolution)

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

        new_data = s1 * u.K  #     self._data.flux._unit   # didn't work
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
            print("    cell_shift:", cell_shift)
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
            # print("new resolution: ",s._resolution)
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
        velocity : `~astropy.units.Quantity`
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
        Write this `Spectrum` as an ~astropy.table.Table.

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

        flux = self.flux
        axis = self.spectral_axis
        mask = self.mask
        t = Table([axis, flux, mask], names=["spectral_axis", "flux", "mask"], meta=self.meta)
        if self.uncertainty is not None:
            t.add_column(self.uncertainty._array, name="uncertainty")
        # f=kwargs.pop("format")
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

        observer_location : `~astropy.coordinates.EarthLocation`
            Location of the observatory. See `~dysh.coordinates.Observatory`.
            This will be transformed to `~astropy.coordinates.ITRS` using the time of observation DATE-OBS or MJD-OBS in `meta`.

        Returns
        -------
        spectrum : `~dysh.spectra.Spectrum`
            The spectrum object
        """
        # @todo add resolution being the channel separation, unless we use the FREQRES column
        warnings.simplefilter("ignore", NoVelocityWarning)
        # @todo generic check_required method since I now have this code in two places (coordinates/core.py).
        # should we also require DATE-OBS or MJD-OBS?
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
            # skip warnings about DATE-OBS being converted to MJD-OBS
            warnings.filterwarnings("ignore", category=FITSFixedWarning)
            # skip warnings FITS keywords longer than 8 chars or containing
            # illegal characters (like _)
            warnings.filterwarnings("ignore", category=VerifyWarning)
            wcs = WCS(header=meta)
            # reset warnings?
        else:
            wcs = None
        # is_topo = is_topocentric(meta["CTYPE1"])  # GBT-specific to use CTYPE1 instead of VELDEF
        target = make_target(meta)
        vc = veldef_to_convention(meta["VELDEF"])
        vf = astropy_frame_dict[decode_veldef(meta["VELDEF"])[1]]
        # could be clever:
        # obstime = Time(meta.get("DATE-OBS",meta.get("MJD-OBS",None)))
        # if obstime is not None: blah blah
        if "DATE-OBS" in meta:
            obstime = Time(meta["DATE-OBS"])
        elif "MJD-OBS" in meta:
            obstime = Time(meta["MJD-OBS"])
        if "DATE-OBS" in meta or "MJD-OBS" in meta:
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
            data,
            wcs=wcs,
            meta=meta,
            velocity_convention=vc,
            radial_velocity=target.radial_velocity,
            rest_value=meta["RESTFRQ"] * u.Hz,
            observer=obsitrs,
            target=target,
        )
        # For some reason, Spectrum1D.spectral_axis created with WCS do not inherit
        # the radial velocity. In fact, they get no radial_velocity attribute at all!
        # This method creates a new spectral_axis with the given radial velocity.
        if observer_location is None:
            s.set_radial_velocity_to(target.radial_velocity)
        # I THINK THIS IS NO LONGER NEEDED
        # if shift_topo:  # and is_topo
        #    vshift = topocentric_velocity_to_frame(target, vf, observer=Observatory["GBT"], obstime=obstime)
        #
        return s

    def _arithmetic_apply(self, other, op, handle_meta, **kwargs):
        if isinstance(other, NDCube):
            result = op(other, **{"handle_meta": handle_meta})
        elif isinstance(other, u.Quantity):
            result = op(other, **{"handle_meta": handle_meta, "meta_other_meta": False})
        elif not isinstance(other, u.Quantity):
            try:
                other = u.Quantity(other, unit=self.unit)
                result = op(other, **{"handle_meta": handle_meta, "meta_other_meta": False})
            except TypeError:
                return NotImplemented
        # result._target = self._target
        # result._observer = self._observer
        # result._velocity_frame = self._velocity_frame
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
        other.spectral_axis.doppler_convention = self.doppler_convention

    def __add__(self, other):
        op = self.add
        handle_meta = self._add_meta
        result = self._arithmetic_apply(other, op, handle_meta)
        return result

    def __sub__(self, other):
        op = self.subtract
        handle_meta = self._add_meta
        result = self._arithmetic_apply(other, op, handle_meta)
        return result

    def __mul__(self, other):
        op = self.multiply
        handle_meta = self._mul_meta
        result = self._arithmetic_apply(other, op, handle_meta)
        return result

    # @todo replace with __truediv__. See issue #241
    def __div__(self, other):
        op = self.divide
        handle_meta = self._div_meta
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


with registry.delay_doc_updates(Spectrum):
    registry.register_writer("ascii.basic", Spectrum, ascii_spectrum_writer_basic)
    registry.register_writer("basic", Spectrum, ascii_spectrum_writer_basic)
    registry.register_writer("ascii.commented_header", Spectrum, ascii_spectrum_writer_commented_header)
    registry.register_writer("commented_header", Spectrum, ascii_spectrum_writer_commented_header)
    registry.register_writer("ascii.fixed_width", Spectrum, ascii_spectrum_writer_fixed_width)
    registry.register_writer("fixed_width", Spectrum, ascii_spectrum_writer_fixed_width)
    registry.register_writer("ascii.ipac", Spectrum, ascii_spectrum_writer_ipac)
    registry.register_writer("ipac", Spectrum, ascii_spectrum_writer_ipac)
    registry.register_writer("votable", Spectrum, spectrum_writer_votable)
