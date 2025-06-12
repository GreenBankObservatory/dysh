"""
Core functions for spectral data.
"""

import warnings
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


def sort_spectral_region(spectral_region):
    """
    Sort the elements of a `~specutils.SpectralRegion`.

    Parameters
    ----------
    spectral_region : `~specutils.SpectralRegion`
        `~specutils.SpectralRegion` to be sorted.

    Returns
    -------
    sorted_spectral_region : `~specutils.SpectralRegion`
        Sorted `~specutils.SpectralRegion`.
    """

    unit = spectral_region.lower.unit
    bound_list = np.sort([srb.value for sr in spectral_region.subregions for srb in sr]) * unit
    it = iter(bound_list)
    sorted_spectral_region = SpectralRegion(list(zip(it, it, strict=False)))

    return sorted_spectral_region


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


def exclude_to_spectral_region(exclude, refspec, fix_exclude=True):
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
    fix_exclude: bool
        If True, fix exclude regions that are out of bounds of the spectral axis to be within the spectral axis.
        Default: True

    Returns
    -------
    sr : `~specutils.SpectralRegion`
        A `~specutils.SpectralRegion` corresponding to `exclude`.
    """

    p = refspec  # noqa: F841
    sa = refspec.spectral_axis
    lastchan = len(sa) - 1

    o = 1
    # If the spectral axis is inverted, flip the order of the elements.
    if refspec.spectral_axis_direction == "decreasing":
        o = -1

    if exclude is not None:
        # A single SpectralRegion was given.
        if isinstance(exclude, SpectralRegion):
            sr = exclude
        # `list` of `int` or `Quantity` or `SpectralRegion` was given.
        else:
            # If user provided a single list, we have to
            # turn it into a list of tuples. If SpectralRegion
            # took a list argument, we wouldn't have to do this.
            if type(exclude[0]) is not tuple:
                it = iter(exclude)
                exclude = list(zip(it, it, strict=False))
            try:
                sr = SpectralRegion(exclude)
                # The above will error if the elements are not quantities.
                # In that case use the spectral axis to define the exclusion regions.
            except ValueError:
                # If the user requested to fix the exclude range.
                if fix_exclude:
                    exclude = np.array(exclude)
                    mask = exclude > len(sa)
                    if mask.sum() > 0:
                        msg = f"Setting upper limit to {lastchan}."
                        exclude[exclude > len(sa)] = lastchan
                        warnings.warn(msg)  # noqa: B028
                # If the spectral_axis is decreasing, flip it.
                sr = SpectralRegion(sa[exclude][:, ::o])

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
        o.append((sr[0].quantity, sr[1].quantity))
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


def exclude_to_mask(exclude, refspec):
    # set a mask based on an exclude region
    # mask ~ exclude_to_indices(exclude_to_region())
    pass


def exclude_to_region_list(exclude, spectrum, fix_exclude=True):
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
    fix_exclude : bool
        See `~spectra.core.exclude_to_spectral_region` for details.
        Default: True

    Returns
    -------
    region_list : list of `~specutils.SpectralRegion`
        Regions defined in `exclude` as a list of `~specutils.SpectralRegion`.
    """

    spectral_region = exclude_to_spectral_region(exclude, spectrum, fix_exclude=fix_exclude)
    spectral_region = spectral_region_to_unit(spectral_region, spectrum)
    region_list = spectral_region_to_list(spectral_region)

    return region_list


def baseline(spectrum, order, exclude=None, exclude_region_upper_bounds=True, **kwargs):
    """Fit a baseline to `spectrum`.
    The code uses `~astropy.modeling.fitting.Fitter` and `~astropy.modeling.polynomial` to compute the baseline.
    See the documentation for those modules for details.

    Parameters
    ----------
    spectrum : `~spectra.spectrum.Spectrum`
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

    Returns
    -------
    model : `~specutils.utils.QuantityModel`
        Best fit model.
        See `~specutils.fitting.fit_continuum`.
    """
    kwargs_opts = {
        #'show': False,
        "model": "chebyshev",
        "fitter": LinearLSQFitter(calc_uncertainties=True),
        "fix_exclude": True,
        "exclude_action": "replace",  # {'replace','append', None}
    }
    kwargs_opts.update(kwargs)

    available_models = {
        "chebyshev": Chebyshev1D,
        "hermite": Hermite1D,
        "legendre": Legendre1D,
        "polynomial": Polynomial1D,
    }
    model = minimum_string_match(kwargs_opts["model"], list(available_models.keys()))
    if model == None:  # noqa: E711
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
    # print(f"MODEL {model} FITTER {fitter}")
    p = spectrum
    if np.isnan(p.data).all():
        # @todo handle masks
        return None  # or raise exception
    if exclude is not None:
        regionlist = exclude_to_region_list(exclude, spectrum, fix_exclude=kwargs_opts["fix_exclude"])
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
             The exposure time, typically given in seconds
         delta_freq : `~numpy.ndarray`, float, or `~astropy.units.Quantity`
             The channel width in frequency units
         tsys : `~numpy.ndarray`, float, or `~astropy.units.Quantity`
             The system temperature, typically in K

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
    if type(weight) == u.Quantity:  # noqa: E721
        return weight.value.astype(np.float64)
    else:
        return weight.astype(np.float64)


def mean_data(data, fedge=0.1, nedge=None, median=False):
    """
    special mean of data to exclude the edges like mean_tsys(), with
    an option to use the median instead of the mean.

    Parameters
    ----------
    data : `~numpy.ndarray`
        The spectral data.
    fedge : float, optional
        Fraction of edge channels to exclude at each end, a number between 0 and 1.
        If `nedge` is used, this parameter is not used.
        Default: 0.1, meaning the central 80% bandwidth is used
    nedge : int, optional
        Number of edge channels to exclude. nedge cannot be 0.
        Default: None, meaning use `fedge`
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

    padded = np.empty(newsize, dtype=y.dtype)
    npad = newsize - nch
    nskip = npad // 2
    padded[nskip : nskip + nch] = y
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

    # Pad if requested.
    # Padding reduces aliasing, but is more expensivep.
    if pad:
        padded, nskip = fft_pad(y)
    else:
        padded, nskip = y.copy(), 0

    new_len = len(padded)

    phase_shift = 2.0 * np.pi * shift / new_len

    yf = padded.copy()
    nan_mask = np.isnan(padded)

    if nan_treatment == "fill":
        yf[nan_mask] = fill_value
    else:
        return NotImplemented

    # Take the inverse FFT.
    ifft = np.fft.ifft(yf)

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
    new_y = shifted_y[nskip : nskip + nch].real
    new_nan_mask = nan_mask[nskip : nskip + nch]

    # Shift NaN elements.
    if keep_nan:
        if abs(shift) < 1:
            new_y[new_nan_mask] = np.nan
        else:
            new_y[np.roll(new_nan_mask, int(shift))] = np.nan

    return new_y


def smooth(data, method="hanning", width=1, kernel=None, show=False):
    """
    Smooth or Convolve a spectrum, optionally decimating it.
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
    kernel : numpy array, optional
        A numpy array which is the kernel by which the signal is convolved.
        Use with caution, as it is assumed the kernel is normalized to
        one, and is symmetric. Since width is ill-defined here, the user
        should supply an appropriate number manually.
        NOTE: not implemented yet.
        The default is None.
    show : bool, optional
        If set, the kernel is returned, instead of the convolved array.
        The default is False.

    Raises
    ------
    Exception
        If no valid smoothing method is given.

    Returns
    -------
    s : `~numpy.ndarray`
        The new convolved spectrum.

    """
    # note that these methods always return odd number in the kernel
    available_methods = {
        "boxcar": Box1DKernel,  # (e.g. width=2 gives hanning)
        "hanning": Trapezoid1DKernel,  # only for width=1
        "gaussian": Gaussian1DKernel,
    }
    method = minimum_string_match(method, list(available_methods.keys()))
    if method == None:  # noqa: E711
        raise ValueError(f"Unrecognized input method {method}. Must be one of {list(available_methods.keys())}")
    kernel = available_methods[method](width)
    if show:
        return kernel
    # the boundary='extend' matches  GBTIDL's  /edge_truncate CONVOL() method
    if hasattr(data, "mask"):
        mask = data.mask
    else:
        mask = None  # noqa: F841
    new_data = convolve(data, kernel, boundary="extend")  # , nan_treatment="fill", fill_value=np.nan, mask=mask)
    return new_data


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

    if abs(fshift) > 1:
        raise ValueError("abs(fshift) must be less than one: {fshift}")

    if method == "fft":
        new_y = fft_shift(y, fshift, pad=pad, window=window)
    elif method == "interpolate":
        new_y = ndimage.shift(y, [fshift])

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
