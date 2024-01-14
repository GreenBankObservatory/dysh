import warnings
from copy import deepcopy

import astropy.units as u
import numpy as np
from astropy.coordinates import SpectralCoord
from astropy.io import registry
from astropy.io.fits.verify import VerifyWarning
from astropy.modeling.fitting import LinearLSQFitter
from astropy.table import Table
from astropy.time import Time
from astropy.wcs import WCS, FITSFixedWarning
from specutils import Spectrum1D

from ..coordinates import (  # is_topocentric,; topocentric_velocity_to_frame,
    KMS,
    Observatory,
    astropy_frame_dict,
    change_ctype,
    decode_veldef,
    get_velocity_in_frame,
    make_target,
    replace_convention,
    sanitize_skycoord,
    veldef_to_convention,
)
from ..plot import specplot as sp
from . import baseline, get_spectral_equivalency


class Spectrum(Spectrum1D):
    """
    This class contains a spectrum and its attributes. It is built on
    `~specutils.Spectrum1D` with added attributes like baseline model.
    Note that `~specutils.Spectrum1D` can contain multiple spectra but
    we probably will not use that because the restriction that it can
    have only one spectral axis conflicts with slight Doppler shifts.
    See `~specutils.Spectrum1D` for the instantiation arguments.

    *Note:* `velocity_convention` should be one of {'radio', 'optical', 'relativistic'}; the  old `~specutils.Spectrum1D` documentation is wrong (there should not be a 'doppler\\_' prefix).
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
        # @TODO - have _observer_location attribute instead
        # and observer property returns getITRS(observer_location,obstime)
        self._observer = kwargs.pop("observer", None)
        Spectrum1D.__init__(self, *args, **kwargs)
        self._spectral_axis._target = self._target
        self._spectral_axis._observer = self._observer
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
        self._exclude_regions = None
        self._plotter = None

    @property
    def exclude_regions(self):
        return self._exclude_regions

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

    def baseline(self, degree, exclude=None, **kwargs):
        """
        Compute and optionally remove a baseline.  The model for the
        baseline can be either a
        `1D polynomial model <https://docs.astropy.org/en/latest/api/astropy.modeling.polynomial.Polynomial1D.html>`_ or a
        `1D Chebyshev polynomial of the first kind <https://docs.astropy.org/en/latest/api/astropy.modeling.polynomial.Chebyshev1D.html>`_.  The code uses `astropy.modeling`
        and `astropy.fitter` to compute the baseline.  See the documentation for those modules.  This method will set the `baseline_model` attribute to the fitted model function which can be evaluated over a domain.

        Parameters
        ----------
            degree : int
                The degree of the polynomial series, a.k.a. baseline order
            exclude : list of 2-tuples of int or ~astropy.units.Quantity, or ~specutils.SpectralRegion
                List of region(s) to exclude from the fit.  The tuple(s) represent a range in the form [lower,upper], inclusive.
                In channel units.

                Examples: One channel-based region: [11,51], Two channel-based regions: [(11,51),(99,123)]. One ~astropy.units.Quantity region: [110.198*u.GHz,110.204*u.GHz]. One compound ~specutils.SpectralRegion: SpectralRegion([(110.198*u.GHz,110.204*u.GHz),(110.196*u.GHz,110.197*u.GHz)]).

                Default: no exclude region

                TODO: Are these OR'd with the existing mask? make that an option
                TODO: Allow these to be Quantities (spectral axis units or equivalent). See list_to_spectral_region()
            model : str
                One of 'polynomial' or 'chebyshev', Default: 'polynomial'
            fitter  :  `~astropy.fitting._FitterMeta`
                The fitter to use. Default: `~astropy.fitter.LinearLSQFitter` (with `calc_uncertaintes=True`).  Be care when choosing a different fitter to be sure it is optimized for this problem.
            remove : bool
                If True, the baseline is removed from the spectrum. Default: False

        """
        kwargs_opts = {
            "remove": False,
            "model": "polynomial",
            "fitter": LinearLSQFitter(calc_uncertainties=True),
        }
        kwargs_opts.update(kwargs)

        self._baseline_model = baseline(self, degree, exclude, **kwargs)
        if kwargs_opts["remove"]:
            s = self.subtract(self._baseline_model(self.spectral_axis))
            self._data = s._data
            self._subtracted = True

    def _undo_baseline(self):
        """
        Undo the most recently computed baseline. If the baseline
        has been subtracted, it will be added back. The `baseline_model`
        attribute is set to None. Exclude regions are untouched.
        """
        if self._baseline_model is None:
            return
        if self._subtracted:
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

                        Examples: One channel-based region: [11,51], Two channel-based regions: [(11,51),(99,123)]. One ~astropy.units.Quantity region: [110.198*u.GHz,110.204*u.GHz]. One compound ~specutils.SpectralRegion: SpectralRegion([(110.198*u.GHz,110.204*u.GHz),(110.196*u.GHz,110.197*u.GHz)]).

        """
        pass

    def list_to_spectral_region(self, inlist):
        # todo utility code to convert a input list of channels or quantities to a spectral region with units of self.spectral_axis.unit. This could go in core.py
        # combine this with _set_exclude_regions
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

    def stats(self):
        """
        Compute some statistics of this `Spectrum`.  The mean, rms,
        data minimum and data maximum are calculated.  Note this works
        with slicing, so, e.g.,  `myspectrum[45:153].stats()` will return
        the statistics of the slice.

        Returns
        -------
        stats : tuple
            Tuple consisting of (mean,rms,datamin,datamax)
            TODO: maybe make this a dict

        """
        mean = self.mean()
        rms = self.data.std()
        dmin = self.min()
        dmax = self.max()
        return (mean, rms, dmin, dmax)

    def smooth(self):
        # todo use specutils.manipulation.smoothing
        # option to smooth baseline too?
        pass

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
            rfq = self.meta["RESTFREQ"] * cunit1
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
        return self._observer

    @property
    def velocity_frame(self):
        return self._velocity_frame

    # This is already in specutils.OneDSpectrumMixin
    # @property
    # def velocity_convention(self):
    #    return self._spectral_axis.doppler_convention

    @property
    def doppler_convention(self):
        return self.velocity_convention

    def axis_velocity(self, unit=KMS):
        """Get the spectral axis in velocity units.
        *Note*: This is not the same as `Spectrum.velocity`, which includes the source radial velocity.
        """
        return self._spectral_axis.to(unit)

    # not needed
    # def velocity_axis_to(self, unit=KMS, toframe=None, doppler_convention=None):
    #    if toframe is not None:
    #        self.set_frame(toframe)
    #    if doppler_convention is not None:
    #        return self._spectral_axis.to(unit=unit, doppler_convention=doppler_convention).to(unit)
    #    else:
    #        return self.velocity.to(unit)

    def get_velocity_shift_to(self, toframe):
        if self._target is None:
            raise Exception("Can't calculate velocity because Spectrum.target is None")
        return get_velocity_in_frame(self._target, toframe, self._observer, self._obstime)

    def set_frame(self, toframe):
        # @TODO VELDEF should be changed as well?
        """Set the sky coordinate and doppler tracking reference frame of this Spectrum. The header 'CTYPE1' will be changed accordingly.

        To make a copy of this Spectrum with new coordinate referece frmae instead, use `with_frame`.

        Parameters
        ----------
        toframe - str
            The coordinate reference frame identifying string, as used by astropy, e.g. 'hcrs', 'icrs', etc.
        """

        if "topo" in toframe:
            actualframe = self.observer  # ???
        else:
            actualframe = astropy_frame_dict.get(toframe, toframe)
        # print(f"actual frame is {actualframe}")
        self._spectral_axis = self._spectral_axis.with_observer_stationary_relative_to(actualframe)
        self._meta["CTYPE1"] = change_ctype(self._meta["CTYPE1"], toframe)
        # @todo shouldn't have two variables for the same thing.
        # Still needed?
        # self.topocentric = self._target.topocentric = "topo" in toframe

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

        s = deepcopy(self)
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
        s = deepcopy(self)
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

    @classmethod
    def make_spectrum(cls, data, meta, use_wcs=True, observer_location=None):
        # , shift_topo=False):
        """Factory method to create a Spectrum object from a data and header.

        Parameters
        ----------
        data :  `~numpy.ndarray`
            The data array. See `~specutils.Spectrum1D`
        meta : dict
            The metadata, typically derived from an SDFITS header.  Required items in `meta` are 'CTYPE[123]','CRVAL[123]', 'CUNIT[123]', 'VELOCITY', 'EQUINOX', 'RADESYS'
        use_wcs : bool
            If True, create a WCS object from `meta`

        observer_location : `~astropy.coordinates.EarthLocation`
            Location of the observatory. See `~dysh.coordinates.Observatory`. This will be transformed to `~astropy.coordinates.ITRS` using the time of observation DATE-OBS or MJD-OBS in `meta`.

        Returns
        -------
        spectrum : `~dysh.spectra.Spectrum`
            The spectrum object
        """
        # @TODO generic check_required method since I now have this code in two places (coordinates/core.py).
        # should we also require DATE-OBS or MJD-OBS?
        _required = set([
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
        ])

        # for k in _required:
        #    print(f"{k} {k in header}")
        if not _required <= meta.keys():
            raise ValueError(f"Header (meta) is missing one or more required keywords: {_required}")

        # @TODO WCS is expensive.
        # Possibly figure how to calculate spectral_axis instead.
        # @TODO allow fix=False in WCS constructor?
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


# @TODO figure how how to document write()
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
