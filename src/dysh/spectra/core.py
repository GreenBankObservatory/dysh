"""
Core functions for spectral data.
"""

from copy import deepcopy
from functools import reduce

import astropy.units as u
import numpy as np
from astropy.convolution import (
    Box1DKernel,
    Gaussian1DKernel,
    Trapezoid1DKernel,
    convolve,
)
from astropy.modeling.fitting import LinearLSQFitter
from astropy.modeling.polynomial import Chebyshev1D, Hermite1D, Legendre1D, Polynomial1D
from scipy import ndimage
from specutils import SpectralRegion
from specutils.fitting import fit_continuum
from specutils.fitting.fitmodels import _strip_units_from_model
from specutils.utils import QuantityModel

from ..log import logger
from ..util import grouper, merge_ranges, minimum_string_match, powerof2

# note that these methods always return odd number in the kernel
_available_smooth_methods = {
    "boxcar": Box1DKernel,  # (e.g. width=2 gives hanning)
    "hanning": Trapezoid1DKernel,  # only for width=1
    "gaussian": Gaussian1DKernel,
}


def available_smooth_methods():
    """The list of smooth methods that dysh understands. These can be passed to various
    smooth routines via their `method` keywords.

    Returns
    -------
    methods: list
        The method names allowable to `smooth`

    """
    return list(_available_smooth_methods.keys())


FWHM_TO_STDDEV = np.sqrt(8 * np.log(2.0))


# @todo: allow data to be SpectrumList or array of Spectrum
def average(data, axis=0, weights=None):
    """Average a group of spectra or scans.

    Parameters
    ----------
    data : `~numpy.ndarray`
        The spectral data, typically with shape (nspect,nchan).
    axis : int
        The axis over which to average the data.  Default axis=0 will return the average spectrum
        if shape is (nspect,nchan)
    weights : `~numpy.ndarray`
        The weights to use in averaging.  These might typically be system temperature based.
        The weights array must be the length of the axis over which the average is taken.
        Default: None will use equal weights

    Returns
    -------
    average : `~numpy.ndarray`
        The average along the input axis
    """
    # Spectra that are blanked will have all channels set to NaN
    # indices = ~np.isnan(data)
    # Find indices with ANY spectra (rows) that have a NaN or Inf
    # goodindices = ~np.ma.fix_invalid(data).mask.any(axis=1)
    # c = data[goodindices]
    # if len(c) == 0:
    #    return np.nan
    # Find indices that have any spectra with all channels = NaN
    # badindices = np.where(np.isnan(data).all(axis=1))
    goodindices = find_non_blanks(data)
    if weights is not None:
        return np.average(data[goodindices], axis, weights[goodindices])
    else:
        return np.average(data[goodindices], axis, weights)


def integration_isnan(data):
    """Helper function to calculate a boolean array that indicates whether
    a collection of integrations is blanked.

    Parameters
    ----------
    data : `~numpy.ndarray`
        The spectral data, typically with shape (nspect,nchan).

    Returns
    -------
    blanks : `~numpy.ndarray`
        Array with length nspect with value True where an integration is blanked.
    """

    return np.isnan(data).all(axis=1)


def find_non_blanks(data):
    """
    Finds the indices of integrations that are not blanked.

    Parameters
    ----------
    data : `~numpy.ndarray`
        The spectral data, typically with shape (nspect,nchan).

    Returns
    -------
    blanks : `~numpy.ndarray`
        Array with indices of non-blanked integrations.
    """

    return np.where(~integration_isnan(data))


def find_blanks(data):
    """
    Finds the indices of blanked integrations

    Parameters
    ----------
    data : `~numpy.ndarray`
        The spectral data, typically with shape (nspect,nchan).

    Returns
    -------
    blanks : `~numpy.ndarray`
        Array with indices of blanked integrations.
    """

    return np.where(integration_isnan(data))


def find_nonblank_ints(cycle1, cycle2, cycle3=None, cycle4=None):
    """
    Find the indices of integrations that are not blanked.

    Parameters
    ----------
    cycle1 : `~numpy.ndarray`
        Data for cycle 1. For example, signal with the noise diode off.
    cycle2 : `~numpy.ndarray`
        Data for cycle 2. For example, reference with the noise diode off.
    cycle3 : `~numpy.ndarray`
        Data for cycle 3. For example, signal with the noise diode on.
        Default is `None`.
    cycle4 : `~numpy.ndarray`
        Data for cycle 4. For example, reference with the noise diode on.
        Default is `None`.

    Returns
    -------
    goodrows : `~numpy.array`
        Indices of the non-blanked rows.
    """

    nb1 = find_non_blanks(cycle1)
    nb2 = find_non_blanks(cycle2)
    if cycle3 is not None:
        nb3 = find_non_blanks(cycle3)
    else:
        nb3 = nb1
    if cycle4 is not None:
        nb4 = find_non_blanks(cycle4)
    else:
        nb4 = nb2
    goodrows = reduce(np.intersect1d, (nb1, nb2, nb3, nb4))

    if len(goodrows) != len(cycle1):
        nblanks = len(cycle1) - len(goodrows)
        logger.info(f"Ignoring {nblanks} blanked integration(s).")

    return goodrows


def sort_spectral_region_subregions(spectral_region):
    """ """

    for i, s in enumerate(spectral_region._subregions):
        spectral_region._subregions[i] = sorted(s)


def invert_spectral_region(sr, refspec):
    """
    Invert an spectral region. The spectral region is sorted and the ranges has been merged previously.

    Parameters
    ----------
    sr : `~specutils.SpectralRegion`
        Spectral region to invert.
    refspec : `~spectra.spectrum.Spectrum`
        Spectrum to define the boundaries for the inversion.

    Returns
    -------
    `~specutils.SpectralRegion`
        Inverted spectral region.
    """
    bl = refspec.spectral_axis.to(sr.bounds[0].unit, equivalencies=refspec.equivalencies).quantity.min()
    bu = refspec.spectral_axis.to(sr.bounds[0].unit, equivalencies=refspec.equivalencies).quantity.max()
    sr._reorder()
    ll = [bl, *list(sum(sr._subregions, ())), bu]
    return exclude_to_spectral_region(list(grouper(ll, 2)), refspec)


def include_to_exclude_spectral_region(include, refspec):
    """
    Convert an inclusion region to an exclude region.

    Parameters
    ----------
    include : list of 2-tuples of int or `~astropy.units.Quantity`, or `~specutils.SpectralRegion`
        List of region(s) to exclude from the fit.  The tuple(s) represent a range in the form [lower,upper], inclusive.
        Examples:

        One channel-based region:

        >>> [11,51]

        Two channel-based regions:

        >>> [(11,51),(99,123)]

        One `~astropy.units.Quantity` region:

        >>> [110.198*u.GHz,110.204*u.GHz].

        One compound `~specutils.SpectralRegion`:

        >>> SpectralRegion([(110.198*u.GHz,110.204*u.GHz),(110.196*u.GHz,110.197*u.GHz)]).

    refspec: `~spectra.spectrum.Spectrum`
        The reference spectrum whose spectral axis will be used
        when converting between `include` and axis units (e.g. channels to GHz).

    Returns
    -------
    exclude_region : `~specutils.SpectralRegion`
        `include` as a region to be excluded.
    """
    spectral_region = exclude_to_spectral_region(include, refspec)
    sort_spectral_region_subregions(spectral_region)
    # Merge include ranges.
    spectral_region = exclude_to_spectral_region(list(merge_ranges(spectral_region._subregions)), refspec)
    # Invert.
    return invert_spectral_region(spectral_region, refspec)


def exclude_to_spectral_region(exclude, refspec):
    """Convert `exclude` to a `~specutils.SpectralRegion`.

    Parameters
    ----------
    exclude : list of 2-tuples of int or `~astropy.units.Quantity`, or `~specutils.SpectralRegion`
        List of region(s) to exclude. The tuple(s) represent a range in the form [lower,upper], inclusive.
        Examples:

        One channel-based region:

        >>> [11,51]

        Two channel-based regions:

        >>> [(11,51),(99,123)]

        One `~astropy.units.Quantity` region:

        >>> [110.198*u.GHz,110.204*u.GHz].

        One compound `~specutils.SpectralRegion`:

        >>> SpectralRegion([(110.198*u.GHz,110.204*u.GHz),(110.196*u.GHz,110.197*u.GHz)]).

    refspec: `~spectra.spectrum.Spectrum`
        The reference spectrum whose spectral axis will be used
        when converting between `exclude` and axis units (e.g. channels to GHz).

    Returns
    -------
    sr : `~specutils.SpectralRegion`
        A `~specutils.SpectralRegion` corresponding to `exclude`.
    """

    sa = refspec.spectral_axis

    # A single SpectralRegion was given.
    if isinstance(exclude, SpectralRegion):
        sr = exclude
    # `list` of `int` or `Quantity` or `SpectralRegion` was given.
    else:
        # If user provided a single list or list of lists,
        # we have to turn it into a list of tuples.
        # If SpectralRegion took a list argument, we wouldn't have to do this.
        # List of lists.
        if isinstance(exclude[0], list):
            exclude = list(map(tuple, exclude))
        # List of scalars or quantities.
        elif not isinstance(exclude[0], tuple):
            it = iter(exclude)
            exclude = list(zip(it, it, strict=False))
        try:
            sr = SpectralRegion(exclude)
            # The above will error if the elements are not quantities.
            # In that case use the spectral axis to define the exclusion regions.
        except ValueError:
            # Make sure all the channels are within bounds.
            exclude = np.array(exclude, dtype=int)
            exclude[exclude >= len(sa)] = len(sa) - 1
            sr = SpectralRegion(sa.quantity[exclude])

    return sr


def spectral_region_to_unit(spectral_region, refspec, unit=None, append_doppler=True):
    """
    Change the unit of `spectral_region` to `unit` using the equivalencies of `refspec`.
    If no `unit` is provided, it will change to the units of `refspec._spectral_axis`.

    Parameters
    ----------
    spectral_region : `~specutils.SpectralRegion`
        `~specutils.SpectralRegion` whose units will be converted.
    refspec : `~spectra.spectrum.Spectrum`
        The reference spectrum whose spectral axis will be used
        when converting to `unit` (e.g. channels to GHz).
    unit : str or `~astropy.units.Quantity`
        The target units for `spectral_region`.
    append_doppler : bool
        Add the `doppler_convention` and `doppler_rest` attributes to the columns of the `~astropy.table.QTable`
        used to convert the `spectral_region` units.

    Returns
    -------
    spectral_region : `~specutils.SpectralRegion`
        SpectralRegion with units of `unit`.
    """

    qt = spectral_region.as_table()

    if unit is None:
        unit = refspec._spectral_axis.unit

    if append_doppler:
        append_doppler_to_spectral_region_qtable(qt, refspec)

    lb = qt["lower_bound"].to(unit, equivalencies=refspec.equivalencies)
    ub = qt["upper_bound"].to(unit, equivalencies=refspec.equivalencies)

    return SpectralRegion(list(zip(lb, ub, strict=False)))


def append_doppler_to_spectral_region_qtable(qtable, refspec):
    """
    Set the `doppler_convention` and `doppler_rest` attributes to the columns of `qtable`.

    Parameters
    ----------
    qtable : `~astropy.table.QTable`
        Table with `~specutils.SpectralAxis` as columns.
    refspec : `~spectra.spectrum.Spectrum`
        The reference spectrum whose spectral axis will be used to set the attributes.
    """
    transfer = {
        "doppler_convention": "doppler_convention",
        "doppler_rest": "rest_value",
    }
    for col in qtable.columns:
        for k, v in transfer.items():
            a = refspec.__getattribute__(v)
            qtable[col].__setattr__(k, a)


def spectral_region_to_list(spectral_region):
    """
    Turn `spectral_region` into a list of `~specutils.SpectralRegion`.
    Each subregion in `spectral_region` will be a list element.

    Parameters
    ----------
    spectral_region : `~specutils.SpectralRegion`
        `~specutils.SpectralRegion` to convert into a list of `~specutils.SpectralRegion`.

    Returns
    -------
    region_list : list of `~specutils.SpectralRegion`
        Subregions of `spectral_region` in a list of `~specutils.SpectralRegion`.
    """

    region_list = []

    # The continuum fitting routines use a list of `SpectralRegion` as input.
    for r in spectral_region.subregions:
        if r[0] != r[1]:
            region_list.append(SpectralRegion([r]))

    return region_list


def spectral_region_to_list_of_tuples(spectral_region):
    """
    Convert `spectral_region` into a list of tuples
    compatible with `~dysh.spectra.Spectrum.baseline`
    `include` or `exclude` arguments.

    Parameters
    ----------
    spectral_region : `~specutils.SpectralRegion`
        Region to convert to a list of tuples.
    """
    o = []
    for sr in spectral_region.subregions:
        o.append((sr[0], sr[1]))
    return o


def region_to_axis_indices(region, refspec):
    """
    Parameters
    ----------
    region : `~specutils.SpectralRegion`
    refspec: `~spectra.spectrum.Spectrum`
        The reference spectrum whose spectral axis will be used
        when converting between exclude and axis units (e.g., channels to GHz).

    Returns
    -------
    indices : 2-tuple of int
        The array indices in `refspec` corresponding to `region.bounds`
    """
    # Spectral region to indices in an input spectral axis.
    # @todo needs to work for multiple spectral regions? or just loop outside this call
    sa = refspec.spectral_axis
    if region.lower.unit != sa.unit:
        # @todo if they are conformable, then allow it and convert
        raise Exception(f"Axis units of region [{region.lower.unit}] and refspec [{sa.unit}] not identical")
    b = [x.value for x in region.bounds]
    indices = np.abs(np.subtract.outer(sa.value, b)).argmin(0)
    return indices


def exclude_to_region_list(exclude, spectrum, clip_exclude=True):
    """
    Convert an exclusion region, `exclude`, to a list of `~specutils.SpectralRegion`.
    This is used for baseline fitting.

    Parameters
    ----------
    exclude : list of 2-tuples of int or `~astropy.units.Quantity`, or `~specutils.SpectralRegion`
        List of region(s) to exclude. The tuple(s) represent a range in the form [lower,upper], inclusive.
        Examples:

        One channel-based region:

        >>> [11,51]

        Two channel-based regions:

        >>> [(11,51),(99,123)]

        One `~astropy.units.Quantity` region:

        >>> [110.198*u.GHz,110.204*u.GHz].

        One compound `~specutils.SpectralRegion`:

        >>> SpectralRegion([(110.198*u.GHz,110.204*u.GHz),(110.196*u.GHz,110.197*u.GHz)]).

    spectrum : `~spectra.spectrum.Spectrum`
        The reference spectrum whose spectral axis will be used
        when converting between `exclude` and axis units (e.g. channels to GHz).
    clip_exclude : bool
        Whether to clip the edges of the exclude regions when they are outside
        `spectrum.spectral_axis`.

    Returns
    -------
    region_list : list of `~specutils.SpectralRegion`
        Regions defined in `exclude` as a list of `~specutils.SpectralRegion`.
    """

    spectral_region = exclude_to_spectral_region(exclude, spectrum)
    spectral_region = spectral_region_to_unit(spectral_region, spectrum)
    sort_spectral_region_subregions(spectral_region)
    if clip_exclude:
        clip_spectral_region_subregions(spectral_region, spectrum)
    region_list = spectral_region_to_list(spectral_region)

    return region_list


def clip_spectral_region_subregions(spectral_region, spectrum):
    """
    Clip the values of the `spectral_region.subregions` if they extend
    outside the `spectrum.spectral_axis`.
    """
    sa_max = spectrum.spectral_axis.quantity.max()
    sa_min = spectrum.spectral_axis.quantity.min()
    for i, s in enumerate(spectral_region.subregions):
        if s[0] < sa_min:
            logger.info(f"{s[0]} is below the minimum spectral axis {sa_min}. Replacing.")
            spectral_region._subregions[i] = (sa_min, s[1])
        if s[1] > sa_max:
            logger.info(f"{s[1]} is above the maximum spectral axis {sa_max}. Replacing.")
            spectral_region._subregions[i] = (s[0], sa_max)


def baseline(spectrum, order, exclude=None, exclude_region_upper_bounds=True, **kwargs):
    """Fit a baseline to `spectrum`.
    The code uses `~astropy.modeling.fitting.Fitter` and `~astropy.modeling.polynomial` to compute the baseline.
    See the documentation for those modules for details.

    Parameters
    ----------
    spectrum : `~dysh.spectra.spectrum.Spectrum`
        The input spectrum.
    order : int
        The order of the polynomial series, a.k.a. baseline order.
    exclude : list of 2-tuples of int or `~astropy.units.Quantity`, or `~specutils.SpectralRegion`
        List of region(s) to exclude from the fit. The tuple(s) represent a range in the form [lower,upper], inclusive.
        Examples:

        One channel-based region:

        >>> [11,51]

        Two channel-based regions:

        >>> [(11,51),(99,123)]

        One `~astropy.units.Quantity` region:

        >>> [110.198*u.GHz,110.204*u.GHz].

        One compound `~specutils.SpectralRegion`:

        >>> SpectralRegion([(110.198*u.GHz,110.204*u.GHz),(110.196*u.GHz,110.197*u.GHz)]).

        Default: no exclude region

    model : str
        One of 'polynomial' or 'chebyshev', Default: 'polynomial'
    fitter : `~astropy.modeling.fitting.Fitter`
        The fitter to use. Default: `~astropy.modeling.fitter.LinearLSQFitter` (with `calc_uncertaintes=True`).
        Be careful when choosing a different fitter to be sure it is optimized for this problem.
    exclude_region_upper_bounds : bool
        Makes the upper bound of any excision region(s) inclusive.
        Allows excising channel 0 for lower-sideband data, and the last channel for upper-sideband data.
    clip_exclude : bool
        Whether to clip the exclude or include regions when they extend outside the `spectrum.spectral_axis`.

    Returns
    -------
    model : `~specutils.utils.QuantityModel`
        Best fit model.
        See `~specutils.fitting.fit_continuum`.
    """
    kwargs_opts = {
        "model": "chebyshev",
        "fitter": LinearLSQFitter(calc_uncertainties=True),
        "clip_exclude": True,
        "exclude_action": "replace",
    }
    kwargs_opts.update(kwargs)

    available_models = {
        "chebyshev": Chebyshev1D,
        "hermite": Hermite1D,
        "legendre": Legendre1D,
        "polynomial": Polynomial1D,
    }
    model = minimum_string_match(kwargs_opts["model"], list(available_models.keys()))
    if model is None:
        raise ValueError(f"Unrecognized input model {kwargs['model']}. Must be one of {list(available_models.keys())}")
    sa_min = spectrum.spectral_axis.min().value
    sa_max = spectrum.spectral_axis.max().value
    selected_model = available_models[model](degree=order, domain=(sa_max, sa_min))

    _valid_exclude_actions = ["replace", "append", None]
    if kwargs_opts["exclude_action"] not in _valid_exclude_actions:
        raise ValueError(
            f"Unrecognized exclude region action {kwargs['exclude_region']}. Must be one of {_valid_exclude_actions}"
        )
    fitter = kwargs_opts["fitter"]
    p = spectrum
    if np.isnan(p.data).all():
        # @todo handle masks
        return None  # or raise exception
    if exclude is not None:
        regionlist = exclude_to_region_list(exclude, spectrum, clip_exclude=kwargs_opts["clip_exclude"])
        if kwargs_opts["exclude_action"] == "replace":
            p._exclude_regions = regionlist
        elif kwargs_opts["exclude_action"] == "append":
            p._exclude_regions.extend(regionlist)
            regionlist = p._exclude_regions
    else:
        # use the spectrum's preset exclude regions if they
        # exist (they will be a list of SpectralRegions or None)
        regionlist = p._exclude_regions

    logger.info(f"EXCLUDING {regionlist}")

    fitted_model = fit_continuum(
        spectrum=p,
        model=selected_model,
        fitter=fitter,
        exclude_regions=regionlist,
        exclude_region_upper_bounds=exclude_region_upper_bounds,
    )

    # astropy's Polynomial1D model does not work well when using
    # a domain other than (-1,1), so we need to remove the units
    # of the model. We use `specutils` methods to do so. This has the
    # downside that the uncertainties of the model are removed.
    if not isinstance(fitted_model, QuantityModel):
        strip_model, _, _ = _strip_units_from_model(fitted_model, p)
        strip_model.domain = fitted_model.domain
        strip_model.window = fitted_model.window
        fitted_model = QuantityModel(strip_model, p.spectral_axis.unit, p.flux.unit)

    return fitted_model


def mean_tsys(calon, caloff, tcal, mode=0, fedge=0.1, nedge=None):
    """
    Get the system temperature from the neighboring calon and caloff, which reflect the state of the noise diode.
    We define an extra way to set the edge size, nedge, if you prefer to use
    number of edge channels instead of the inverse fraction.  This implementation recreates GBTIDL's `dcmeantsys`.

    Parameters
    ----------
        calon : `~numpy.ndarray`-like
            ON calibration

        caloff  :  `~numpy.ndarray`-like
            OFF calibration

        tcal  :  `~numpy.ndarray`-like
            calibration temperature

        mode : int
            mode=0  Do the mean before the division
            mode=1  Do the mean after the division

        fedge : float
            Fraction of edge channels to exclude at each end, a number between 0 and 1. Default: 0.1, meaning the central 80% bandwidth is used

        nedge : int
            Number of edge channels to exclude. Default: None, meaning use `fedge`

    Returns
    -------
        meanTsys : `~numpy.ndarray`-like
            The mean system temperature
    """
    # @todo Pedro thinks about a version that takes a spectrum with multiple SpectralRegions to exclude.
    nchan = len(calon)
    if nedge is None:
        nedge = int(nchan * fedge)
    # Python uses exclusive array ranges while GBTIDL uses inclusive ones.
    # Therefore we have to add a channel to the upper edge of the range
    # below in order to reproduce exactly what GBTIDL gets for Tsys.
    # See github issue #28.
    # Define the channel range once.
    chrng = slice(nedge, -(nedge - 1), 1)

    # Make them doubles. Probably not worth it.
    caloff = caloff.astype("d")
    calon = calon.astype("d")

    if mode == 0:  # mode = 0 matches GBTIDL output for Tsys values
        meanoff = np.nanmean(caloff[chrng])
        meandiff = np.nanmean(calon[chrng] - caloff[chrng])
        meanTsys = meanoff / meandiff * tcal + tcal / 2.0
    else:
        meanTsys = np.nanmean(caloff[chrng] / (calon[chrng] - caloff[chrng]))
        meanTsys = meanTsys * tcal + tcal / 2.0

    # meandiff can sometimes be negative, which makes Tsys negative!
    # GBTIDL also takes abs(Tsys) because it does sqrt(Tsys^2)
    return np.abs(meanTsys)


def sq_weighted_avg(a, axis=0, weights=None):
    # @todo make a generic moment or use scipy.stats.moment
    r"""Compute the mean square weighted average of an array (2nd moment).

    :math:`v = \sqrt{\frac{\sum_i{w_i~a_i^{2}}}{\sum_i{w_i}}}`

    Parameters
    ----------
    a : `~numpy.ndarray`
        The data to average.
    axis : int
        The axis over which to average the data.  Default: 0
    weights : `~numpy.ndarray` or None
        The weights to use in averaging.  The weights array must be the
        length of the axis over which the average is taken.  Default:
        `None` will use equal weights.

    Returns
    -------
    average : `~numpy.ndarray`
        The square root of the squared average of a along the input axis.
    """
    # if weights is None:
    #    w = np.ones_like(a)
    # else:
    #    w = weights
    a2 = a * a
    v = np.sqrt(np.average(a2, axis=axis, weights=weights))
    return v


def tsys_weight(exposure, delta_freq, tsys):
    r"""Compute the system temperature based weight(s) using the expression:
        :math:`w = t_{exp} \times \delta_\nu / T_{sys}^2,`

    If `exposure`, `delta_freq`, or `tsys` parameters are given as `~astropy.units.Quantity`,
    they will be converted to seconds, Hz, K, respectively for the calculation.
    (Not that this really matters since weights are relative to each other)


     >>> import astropy.units as u
     >>> [3*u.s, 4*u.s]   ### WRONG

    Rather, if using `~astropy.units.Quantity`, they have to be `~astropy.units.Quantity` objects, e.g.,

     >>> [3, 4]*u.s ### RIGHT!

    Parameters
    ----------
     exposure : `~numpy.ndarray`, float, or `~astropy.units.Quantity`
         The exposure time, typically given in seconds.
     delta_freq : `~numpy.ndarray`, float, or `~astropy.units.Quantity`
         The channel width in frequency units.
     tsys : `~numpy.ndarray`, float, or `~astropy.units.Quantity`
         The system temperature, typically in K.

    Returns
    -------
         weight : `~numpy.ndarray`
             The weights array
    """

    # Quantitys work with abs and power!
    # Using `numpy.power` like this results in increased
    # precision over the calculation used by GBTIDL:
    # weight = abs(delta_freq) * exposure / tsys**2.
    weight = abs(delta_freq) * exposure * np.power(tsys, -2.0)
    if isinstance(weight, u.Quantity):
        return weight.value
    else:
        return weight


def mean_data(data, fedge=0.1, nedge=None, median=False):
    """
    Special mean of data to exclude the edges like mean_tsys(), with
    an option to use the median instead of the mean.

    Parameters
    ----------
    data : `~numpy.ndarray`
        The spectral data.
    fedge : float, optional
        Fraction of edge channels to exclude at each end, a number between 0 and 1.
        If `nedge` is used, this parameter is not used.
        Default: 0.1, meaning the central 80% bandwidth is used.
    nedge : int, optional
        Number of edge channels to exclude. nedge cannot be 0.
        Default: None, meaning use `fedge`.
    median : boolean, optional
        Use the median instead of the mean.
        The default is False.

    Returns
    -------
    meandata : float

    """

    nchan = len(data)
    if nedge is None:
        nedge = int(nchan * fedge)
    chrng = slice(nedge, -(nedge - 1), 1)
    if median:
        meandata = np.nanmedian(data[chrng])
    else:
        meandata = np.nanmean(data[chrng])
    return meandata


def get_spectral_equivalency(restfreq, velocity_convention):
    # Yeesh, the doppler_convention parameter for SpectralAxis.to does not match the doppler_convention list for Spectrum1D!
    # This is actually bug in Spectrum1D documentation https://github.com/astropy/specutils/issues/1067
    if "radio" in velocity_convention:
        return u.doppler_radio(restfreq)
    elif "optical" in velocity_convention:
        return u.doppler_optical(restfreq)
    elif "relativistic" in velocity_convention:
        return u.doppler_relativistic(restfreq)
    elif "redshift" in velocity_convention:
        return u.doppler_redshift()
    else:
        raise ValueError(f"Unrecognized velocity convention {velocity_convention}")


def fft_pad(y):
    """
    Pad signal `y` to the next power of 2 using its edge values.

    Parameters
    ----------
    y : 1D `~numpy.ndarray`
        Signal to be padded.

    Returns
    -------
    padded : 1D `~numpy.ndarray`
        Padded signal.
    nskip : int
        Number of samples added to each end of the padded signal.
    """

    nch = len(y)
    pow2 = powerof2(nch)
    pow2 += 1
    newsize = 2**pow2

    padded = np.ma.empty(newsize, dtype=y.dtype)
    npad = newsize - nch
    nskip = npad // 2
    padded[nskip : nskip + nch] = y
    padded.mask[nskip : nskip + nch] = y.mask
    padded[0:nskip] = padded[nskip]
    padded[nskip + nch :] = padded[nskip + nch]

    return padded, nskip


def fft_shift(
    y,
    shift,
    pad=True,
    window=True,
    nan_treatment="fill",
    fill_value=0,
    keep_nan=True,
):
    """
    Shift a signal `y` by `shift` amount using a phase shift.
    This requires taking the inverse FFT of the signal, shifting its phase,
    and then taking the FFT to shift the signal.

    Parameters
    ----------
    y : 1D `~numpy.ndarray`
        Signal to be shifted. Only 1D supported.
    shift : float
        Amount to shift the signal by in channels.
    pad : bool
        Pad the signal to prevent aliasing.
    window : bool
        Apply a Welch window to prevent ringing.
    nan_treatment : 'fill'
        'fill' replaces NaN values with `fill_value` before the FFT.
    fill_value : float
        Value used to replace NaN values. Used if `nan_treatment=='fill'`.
    keep_nan : bool
        If `True` the output will keep the NaN values. If `False` NaN values will be filled.

    Returns
    -------
    new_y : 1D `~numpy.ndarray`
        Phase shifted version of the original signal `y`.
    """

    nch = len(y)
    y = deepcopy(y)

    # Pad if requested.
    # Padding reduces aliasing, but it requires more resources.
    if pad:
        padded, nskip = fft_pad(y)
    else:
        padded, nskip = deepcopy(y), 0

    new_len = len(padded)

    phase_shift = 2.0 * np.pi * shift / new_len

    yf = padded.copy()
    nan_mask = np.isnan(padded.data) | padded.mask

    if nan_treatment == "fill":
        yf[nan_mask] = fill_value
    else:
        return NotImplemented

    # Take the inverse FFT.
    ifft = np.fft.ifft(yf.filled(fill_value))

    # Split into amplitude and phase.
    amp = abs(ifft)
    phs = np.angle(ifft)

    # Array of indices for computing the phase shift at each channel.
    harr = np.arange(0, new_len // 2)
    # Index the phase shift samples this way to
    # get a Hermite-symmetric complex exponential.
    farr = np.hstack((harr, harr - new_len // 2))

    # Shift and apply window if needed.
    phs += farr * phase_shift
    if window:
        amp *= 1.0 - (farr / float(new_len / 2.0)) ** 2

    # Combine again and take the FFT.
    shifted_ifft = amp * np.cos(phs) + 1j * amp * np.sin(phs)
    shifted_y = np.fft.fft(shifted_ifft)

    # Remove padding and take the real part.
    new_y = np.ma.masked_invalid(shifted_y.real[nskip : nskip + nch])
    new_nan_mask = nan_mask[nskip : nskip + nch]

    # Shift NaN elements and mask.
    if keep_nan:
        if abs(shift) < 1:
            new_nan_mask = mask_fshift(new_nan_mask, shift)
            new_y[new_nan_mask] = np.nan
            new_y.mask[new_nan_mask] = True
        else:
            new_nan_mask = np.roll(new_nan_mask, int(shift))
            new_y[new_nan_mask] = np.nan
            new_y.mask[new_nan_mask] = True

    return new_y


def decimate(data, n, meta=None):
    """
    Decimate a data array by `n` pixels.

    Parameters
    ----------
    data : `~numpy.ndarray` or `~astropy.quantity.Quantity`
        The data to decimate.
    n : int
        Decimation factor of the spectrum by returning every n-th channel.
    meta : dict
        Metadata dictionary with CDELT1, CRVAL1, NAXIS1, and CRPIX1 which will be recalculated.

    Returns
    -------
    tuple : `~numpy.ndarray` or `~astropy.quantity.Quantity` and dict
        A tuple of the decimated `data` and updated metadata (or None if no `meta` given).
    """

    if not float(n).is_integer():
        raise ValueError(f"`n` ({n}) must be an integer.")

    nchan = len(data)
    idx = np.arange(0, nchan, n)
    new_data = data[idx]  # this will also decimate the mask if data has a mask attribute.
    if meta is not None:
        new_meta = deepcopy(meta)
        new_cdelt1 = meta["CDELT1"] * n
        cell_shift = 0.5 * (n - 1) * meta["CDELT1"]

        new_meta["CDELT1"] = new_cdelt1
        new_meta["CRPIX1"] = 1.0 + (meta["CRPIX1"] - 1) / n + 0.5 * (n - 1) / n
        new_meta["CRVAL1"] += cell_shift
        new_meta["NAXIS1"] = len(new_data)
    else:
        new_meta = None

    return (new_data, new_meta)


# @todo it would be nice if this could take a 2-D array of N spectra. astropy.convolve can handle it.
def smooth(
    data,
    method="hanning",
    width=1,
    ndecimate=0,
    meta=None,
    kernel=None,
    mask=None,
    boundary="extend",
    nan_treatment="fill",
    fill_value=np.nan,
    preserve_nan=True,
):
    """
    Smooth or Convolve spectrum, optionally decimating it.
    A number of methods from astropy.convolution can be selected
    with the method= keyword.

    Default smoothing is hanning.

    Parameters
    ----------
    data : `~numpy.ndarray`
        Input data array to smooth. Note smoothing array does not need a
        WCS since it is channel based.
    method : string, optional
        Smoothing method. Valid are: 'hanning', 'boxcar' and
        'gaussian'. Minimum match applies.
        The default is 'hanning'.
    width : int or float, optional
        Effective width of the convolving kernel.  Should ideally be an
        odd number.
        For 'hanning' this should be 1, with a 0.25,0.5,0.25 kernel.
        For 'boxcar' an even value triggers an odd one with half the
        signal at the edges, and will thus not reproduce GBTIDL.
        For 'gaussian' this is the FWHM of the final beam in channels (float). We normally
        assume the input beam has FWHM=1, pending resolution on cases
        where CDELT1 is not the same as FREQRES.
        The default is 1.
    ndecimate : int, optional
        Decimation factor of the spectrum by returning every `ndecimate`-th channel.
    meta: dict
         metadata dictionary with CDELT1, CRVAL1, CRPIX1, NAXIS1, and FREQRES which will be recalculated if necessary
    kernel : `~numpy.ndarray`, optional
        A numpy array which is the kernel by which the signal is convolved.
        Use with caution, as it is assumed the kernel is normalized to
        one, and is symmetric. Since width is ill-defined here, the user
        should supply an appropriate number manually.
        NOTE: not implemented yet.
        The default is None.
    mask : None or `~numpy.ndarray`, optional
        A "mask" array.  Shape must match ``array``, and anything that is masked
        (i.e., not 0/`False`) will be set to NaN for the convolution.  If
        `None`, no masking will be performed unless ``array`` is a masked array.
        If ``mask`` is not `None` *and* ``array`` is a masked array, a pixel is
        masked if it is masked in either ``mask`` *or* ``array.mask``.
    boundary : str, optional
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
    fill_value : float, optional
        The value to use outside the array when using ``boundary='fill'``. Default value is ``NaN``.
    nan_treatment : {'interpolate', 'fill'}, optional
        The method used to handle NaNs in the input ``array``:
            * ``'interpolate'``: ``NaN`` values are replaced with
              interpolated values using the kernel as an interpolation
              function. Note that if the kernel has a sum equal to
              zero, NaN interpolation is not possible and will raise an
              exception.
            * ``'fill'``: ``NaN`` values are replaced by ``fill_value``
              prior to convolution.
    preserve_nan : bool, optional
        After performing convolution, should pixels that were originally NaN
        again become NaN?

    Raises
    ------
    Exception
        If no valid smoothing method is given.

    Returns
    -------
    s : `~numpy.ndarray`
        The new convolved spectrum.

    """
    if kernel is not None:
        raise NotImplementedError("Custom kernels are not yet implemented.")

    asm = available_smooth_methods()
    method = minimum_string_match(method, asm)
    if method is None:
        raise ValueError(f"Unrecognized input method {method}. Must be one of {asm}")

    if not float(ndecimate).is_integer():
        raise ValueError("`decimate ({ndecimate})` must be an integer.")

    kernel = _available_smooth_methods[method](width)
    # Notes:
    # 1. the boundary='extend' matches  GBTIDL's  /edge_truncate CONVOL() method
    # 2. no need to pass along a mask to convolve if the data have a mask already. astropy will obey the data mask
    # 3. However, astropy convolve will not return a masked array even if the input data have a mask.
    # 4. We create an input mask if the data do not have one and we ensure input NaNs that are smoothed to output Nans get masked.
    # 5. We then mask any NaN on output by modifying the input mask
    if hasattr(data, "mask"):
        mask = data.mask
    else:
        mask = np.full(data.shape, False)

    new_data = convolve(
        data,
        kernel,
        boundary=boundary,
        mask=mask,
        nan_treatment=nan_treatment,
        fill_value=fill_value,
        preserve_nan=preserve_nan,
    )
    mask[np.where(np.isnan(new_data))] = True
    new_data = np.ma.masked_array(new_data, mask)
    new_meta = deepcopy(meta)
    if new_meta is not None:
        if method == "gaussian":
            width = width * FWHM_TO_STDDEV
        new_meta["FREQRES"] = np.sqrt((width * new_meta["CDELT1"]) ** 2 + new_meta["FREQRES"] ** 2)
    if ndecimate > 0:
        new_data, new_meta = decimate(new_data, n=ndecimate, meta=new_meta)

    return new_data, new_meta


def data_ishift(y, ishift, axis=-1, remove_wrap=True, fill_value=np.nan):
    """
    Shift `y` by `ishift` channels, where `ishift` is a natural number.

    Parameters
    ----------
    y : array
        Data to be shifted.
    ishift : int
        Amount to shift data by.
    axis : int
        Axis along which to apply the shift.
    remove_wrap : bool
        Replace channels that wrap around with `fill_value`.
    fill_value : float
        Value used to replace the data in channels that wrap around after the shift.

    Returns
    -------
    new_y : array
        Shifted `y`.
    """

    new_y = np.roll(y, ishift, axis=axis)

    if remove_wrap:
        if ishift < 0:
            new_y[ishift:] = fill_value
        else:
            new_y[:ishift] = fill_value

    return new_y


def data_fshift(y, fshift, method="fft", pad=False, window=True):
    """
    Shift `y` by `fshift` channels, where `abs(fshift)<1`.

    Parameters
    ----------
    y : array
        Data to be shifted.
    fshift : float
        Amount to shift the data by.
        abs(fshift) must be less than 1.
    method : "fft" | "interpolate"
        Method to use for shifting.
        "fft" uses a phase shift.
        "interpolate" uses `scipy.ndimage.shift`.
    pad : bool
        Pad the data during the phase shift.
        Only used if `method="fft"`.
    window : bool
        Apply a Welch window during phase shift.
        Only used if `method="fft"`.
    """

    if abs(fshift) >= 1:
        raise ValueError("abs(fshift) must be less than one: {fshift}")

    if method == "fft":
        new_y = fft_shift(y, fshift, pad=pad, window=window)
    elif method == "interpolate":
        new_y = ndimage.shift(y.filled(0.0), [fshift])
        mask = mask_fshift(y.mask, fshift)
        new_y = np.ma.masked_where(mask, new_y)

    return new_y


def data_shift(y, s, axis=-1, remove_wrap=True, fill_value=np.nan, method="fft", pad=False, window=True):
    """
    Shift `y` by `s` channels.

    Parameters
    ----------
    y : array
        Data to be shifted.
    s : float
        Amount to shift the data by.
    axis : int
        Axis along which to apply the shift.
    remove_wrap : bool
        Replace channels that wrap around with `fill_value`.
    fill_value : float
        Value used to replace the data in channels that wrap around after the shift.
    method : "fft" | "interpolate"
        Method to use for shifting.
        "fft" uses a phase shift.
        "interpolate" uses `scipy.ndimage.shift`.
    pad : bool
        Pad the data during the phase shift.
        Only used if `method="fft"`.
    window : bool
        Apply a Welch window during phase shift.
        Only used if `method="fft"`.
    """

    ishift = int(np.round(s))  # Integer shift.
    fshift = s - ishift  # Fractional shift.

    logger.debug(f"Shift: s={s}  ishift={ishift} fshift={fshift}")

    if ishift != 0:
        # Apply integer shift.
        y_new = data_ishift(y, ishift, axis=axis, remove_wrap=remove_wrap, fill_value=fill_value)
    else:
        y_new = deepcopy(y)

    if fshift != 0:
        # Apply fractional shift.
        y_new = data_fshift(y_new, fshift, method=method, pad=pad, window=window)

    return y_new


def mask_fshift(mask, fshift):
    """
    Shift `mask` by `fshift` channels.
    This should only be used when `abs(fshift)<1`.
    It expands the mask using binary dilation
    to account for the spread of the masked values when
    they are shifted by a fractional number of channels.

    Parameters
    ----------
    mask : array_like
        Array with masked values. Either ones and zeros or True and False.
    fshift : float
        Amount to shift by.

    Returns
    -------
    new_mask : array_like
        `mask` shifted by `fshift` channels.

    Raises
    ------
    ValueError
        If `abs(fshift)` is greather than 1.
    """
    if abs(fshift) >= 1:
        raise ValueError(f"abs(fshift) greater than 1 ({fshift=})")
    if fshift < 0:
        structure = np.array([1, 1, 0], dtype=bool)
    elif fshift > 0:
        structure = np.array([0, 1, 1], dtype=bool)
    else:
        structure = np.array([0, 1, 0], dtype=bool)
    new_mask = ndimage.binary_dilation(mask, structure=structure)
    return new_mask


def cog_slope(c, flat_tol=0.1):
    """
    Find slope for curve of growth analysis.

    Parameters
    ----------
    c : array
        Curve of growth values.
    flat_tol : float
        Tolerance to define the flat portion of the curve of growth.
        The flat portion is that which is zero within `flat_tol` times the rms
        of the slope.

    Returns
    -------
    slope : array
        Slope of `c`.
    slope_rms : float
        Standard deviation of the slope.
    flat_idx0 : int
        Index where the slope is consistent with zero.
    """
    slope = np.diff(c)
    slope_rms = slope.std()
    flat = abs(slope) < flat_tol * slope_rms
    flat_idx = np.where(flat)[0]
    try:
        flat_idx0 = flat_idx[0]
    except IndexError:
        flat_idx0 = -1

    return slope, slope_rms, flat_idx0


def cog_flux(c, flat_tol=0.1):
    """
    Find the flux from the curve of growth.

    Parameters
    ----------
    c : array
        Curve of growth values.
    flat_tol : float
        Tolerance to define the flat portion of the curve of growth.
        The flat portion is that which is zero within `flat_tol` times the rms
        of the slope.

    Returns
    -------
    flux : float
        The median of the curve of growth over its flat portion.
    flux_std : float
        The standard deviation of the curve of growth over the flat portion.
    slope : array
        The median value of the slope for the curve of growth before it becomes flat.
    """
    slope, _slope_rms, flat_idx0 = cog_slope(c, flat_tol)
    flux = np.nanmedian(c.filled(np.nan)[flat_idx0:])
    flux_std = np.nanstd(c[flat_idx0:])
    slope = np.nanmedian(slope.filled(np.nan)[:flat_idx0])
    return flux, flux_std, slope


def curve_of_growth(x, y, vc=None, width_frac=None, bchan=None, echan=None, flat_tol=0.1, fw=1) -> dict:
    """
    Curve of growth analysis based on Yu et al. (2020) [1]_.

    Parameters
    ----------
    x : array
        Velocity values.
    y : array
        Flux values.
    vc : float
        Central velocity of the line in the same units as `x`.
        If not provided, it will be estimated from the moment 1 of the `x` and `y` values.
    width_frac : list
        List of fractions of the total flux at which to compute the line width.
        If 0.25 and 0.85 are not included, they will be added to estimate the concentration
        as defined in [1]_.
    bchan : None or int
        Beginning channel where there is signal.
        If not provided it will be estimated using `fw` times the width of the line at the largest `width_frac`.
    echan : None or int
        End channel where there is signal.
        If not provided it will be estimated using `fw` times the width of the line at the largest `width_frac`.
    flat_tol : float
        Tolerance used to define the flat portion of the curve of growth.
        The curve of growth will be considered flat when it's slope is within `flat_tol` times the standard deviation of the slope from zero.
    fw : float
        When estimating the line-free range, use `fw` times the largest width.

    Returns
    -------
    results : dict
        Dictionary with the flux (:math:`F`), width (:math:`V`), flux asymmetry (:math:`A_F`), slope asymmetry (:math:`A_C`), concentration (:math:`C_V`),
        rms, central velocity ("vel"), redshifted flux (:math:`F_r`), blueshifted flux (:math:`F_b`),`bchan` and `echan`, and errors on the flux ("flux_std"), width ("width_std"), central velocity ("vel_std"), redshifted flux ("flux_r_std"), and blueshifted flux ("flux_b_std").
        The rms is the standard deviation in the range outside of (bchan,echan).

    .. [1] `N. Yu, L. Ho & J. Wang, "On the Determination of Rotation Velocity and Dynamical Mass of Galaxies Based on Integrated H I Spectra"
       <https://ui.adsabs.harvard.edu/abs/2020ApJ...898..102Y/abstract>`_.
    """

    # Save units for latter.
    # These will be stripped,
    # as quantities and masked arrays are not compatible.
    y_unit = y.unit
    x_unit = x.unit

    if width_frac is None:
        width_frac = [0.25, 0.65, 0.75, 0.85, 0.95]
    # Sort data values.
    p = 1
    if x[0] > x[1]:
        p = -1
    # Strip units.
    _x = x[::p].value
    _y = np.ma.masked_invalid(y[::p].value)
    # Use channel ranges if provided.
    # Slice the end first to keep the meaning of bchan.
    if echan is not None:
        _x = _x[:echan]
        _y = _y[:echan]
    if bchan is not None:
        _x = _x[bchan:]
        _y = _y[bchan:]
    dx = np.diff(_x)
    ydx = _y[:-1] * dx

    # Find initial guess for central velocity.
    if vc is not None:
        _vc = vc.to(x_unit).value
    else:
        _vc = (_x * _y).sum() / _y.sum()
    vc_idx = np.argmin(abs(_x - _vc))

    # Compute curve of growth.
    b = np.ma.cumsum(ydx[:vc_idx][::-1])  # Blue.
    bp = np.ma.cumsum(ydx[: vc_idx + 1][::-1])
    r = np.ma.cumsum(ydx[vc_idx:])  # Red.
    rp = np.ma.cumsum(ydx[vc_idx + 1 :])
    s = min(len(b), len(r))
    t = b[:s] + r[:s]
    dx = dx[:s]
    xr = _x[vc_idx:] - _vc
    xb = abs(_x[:vc_idx] - _vc)[::-1]
    s = min(len(xb), len(xr))
    xt = xr[:s] + xb[:s]
    # todo: use these to give widths in each direction.
    # vb = x[:vc_idx][::-1]
    # vr = x[vc_idx - 1 :]

    # Find flux.
    # Empirically, fluxes are in error by <3%.
    flux, flux_std, _slope = cog_flux(t, flat_tol)
    flux_std = np.sqrt(flux_std**2 + (1.03 * flux_std) ** 2)
    # Blue.
    flux_b, flux_b_std, slope_b = cog_flux(b, flat_tol)
    flux_bp, _, _ = cog_flux(bp, flat_tol)
    dflux_b = flux_b - flux_bp
    flux_b_std = np.sqrt(flux_b_std**2 + (1.03 * flux_b_std) ** 2 + dflux_b**2)
    # Red.
    flux_r, flux_r_std, slope_r = cog_flux(r, flat_tol)
    flux_rp, _, _ = cog_flux(rp, flat_tol)
    dflux_r = flux_r - flux_rp
    flux_r_std = np.sqrt(flux_r_std**2 + (1.03 * flux_r_std) ** 2 + dflux_r**2)

    # Find line widths.
    if 0.25 not in width_frac:
        width_frac.insert(0, 0.25)
    if 0.85 not in width_frac:
        width_frac.append(0.85)
    widths = dict.fromkeys(width_frac)
    _, _, flat_idx0 = cog_slope(t, flat_tol)
    nt = t[:flat_idx0] / flux
    for f in width_frac:
        idx = np.nanargmin(abs(nt - f))
        widths[f] = xt[idx]

    # Estimate rms from line-free channels.
    if bchan is None:
        _bchan = np.nanargmin(abs(x.value - (_vc - fw * widths[max(width_frac)])))
        if _bchan <= 0:
            _bchan = 0
    else:
        _bchan = bchan
    if echan is None:
        _echan = np.nanargmin(abs(x.value - (_vc + fw * widths[max(width_frac)])))
        if _echan >= len(x):
            _echan = len(x)
    else:
        _echan = echan
    _bchan, _echan = np.sort([_bchan, _echan])

    # Use y values without channel crop.
    rms = np.nanstd(np.hstack((y.value[:_bchan], y.value[_echan:])))

    # Estimate error on line centroid.
    vc_std = 0
    if vc is None:
        fac1 = _vc / (_y * _x).sum() * np.sqrt(np.sum((_x * rms) ** 2))
        fac2 = np.sqrt(len(_x)) * rms * _vc / _y.sum()
        vc_std = np.sqrt(fac1**2 + fac2**2)

    # Estimate error on widths.
    nt_std = np.sqrt((nt / t[:flat_idx0] * rms * dx[:flat_idx0]) ** 2 + (nt / flux * flux_std) ** 2)
    widths_std = dict.fromkeys(width_frac)
    for f in width_frac:
        idx = np.nanargmin(abs(nt - f))
        std = nt_std[idx]
        idx_m = np.nanargmin(abs(nt - std - f))
        idx_p = np.nanargmin(abs(nt + std - f))
        std_m = xt[idx_m] - xt[idx]
        std_p = xt[idx] - xt[idx_p]
        widths_std[f] = np.nanmax((std_m, std_p, dx[idx]))
        widths_std[f] = np.sqrt(
            widths_std[f] ** 2 + (widths[f] * 0.01) ** 2
        )  # Empirically, the widths are in error by <1%.

    # Flux asymmetry.
    a_f = flux_b / flux_r
    if a_f < 1:
        a_f = 1 / a_f

    # Shape asymmetry.
    a_c = slope_b / slope_r
    if a_c < 1:
        a_c = 1 / a_c

    # Concentration.
    c_v = widths[0.85] / widths[0.25]

    flux_unit = y_unit * x_unit

    results = {
        "flux": flux * flux_unit,
        "flux_std": flux_std * flux_unit,
        "flux_r": flux_r * flux_unit,
        "flux_r_std": flux_r_std * flux_unit,
        "flux_b": flux_b * flux_unit,
        "flux_b_std": flux_b_std * flux_unit,
        "width": {k: v * x_unit for k, v in widths.items()},
        "width_std": {k: v * x_unit for k, v in widths_std.items()},
        "A_F": a_f,
        "A_C": a_c,
        "C_V": c_v,
        "rms": rms * y_unit,
        "bchan": _bchan,
        "echan": _echan,
        "vel": _vc * x_unit,
        "vel_std": vc_std * x_unit,
    }

    return results
