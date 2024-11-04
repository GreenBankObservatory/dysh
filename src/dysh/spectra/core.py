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

from ..coordinates import veltofreq
from ..log import log_function_call, logger
from ..util import minimum_string_match, powerof2


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


def exclude_to_region(exclude, refspec, fix_exclude=False):
    """Convert an exclude list to a list of ~specutuls.SpectralRegion.

    Parameters
    ----------
    exclude : list of 2-tuples of int or `~astropy.units.Quantity`, or `~specutils.SpectralRegion`
        List of region(s) to exclude from the fit.  The tuple(s) represent a range in the form [lower,upper], inclusive.    Examples:

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
        when converting between exclude and axis units (e.g. channels to GHz).
    fix_exclude: bool
        If True, fix exclude regions that are out of bounds of the specctral axis to be within the spectral axis. Default:False

      Returns
      -------
      regionlist : list of `~specutil.SpectralRegion`
        A list of `~specutil.SpectralRegion` corresponding to `exclude` with units of the `refspec.spectral_axis`.
    """
    regionlist = []
    p = refspec
    sa = refspec.spectral_axis
    if exclude is not None:
        regionlist = []
        # a single SpectralRegion was given
        if isinstance(exclude, SpectralRegion):
            b = exclude.bounds
            if b[0] < sa[0] or b[1] > sa[1]:
                msg = f"Exclude limits {b} are not fully within the spectral axis {sa}"
                raise Exception(msg)
            regionlist.append(exclude)
        # list of int or Quantity or SpectralRegion was given
        else:
            # if user provided a single list, we have to
            # add another set of brackets so we an iterate.
            # If SpectralRegion took a list argument, we wouldn't
            # have to do this.
            if len(np.shape(exclude[0])) == 0:
                exclude = [exclude]
            # NB: we are assuming that a SpectralAxis is always [lower...upper].  Is this true???
            for pair in exclude:
                if type(pair[0]) == int:
                    # convert channel to spectral axis units
                    lastchan = len(sa) - 1
                    msg = f"Exclude limits {pair} are not fully within the spectral axis [0,{lastchan}]."
                    if pair[0] < 0 or pair[1] > lastchan:
                        if fix_exclude:
                            msg += f" Setting upper limit to {lastchan}."
                            pair[1] = lastchan
                            warnings.warn(msg)
                        else:
                            raise Exception(msg)
                    pair = [sa[pair[0]], sa[pair[1]]]
                # if it is already a spectral region no additional
                # work is needed
                # @todo we should test that the SpectralRegion is not out of bounds
                if isinstance(pair[0], SpectralRegion):
                    b = pair[0].bounds
                    if b[0] < sa[0] or b[1] > sa[1]:
                        msg = f"Exclude limits {pair} are not fully within the spectral axis {p.spectral_axis}"
                        raise Exception(msg)
                    regionlist.append(pair)
                else:  # it is a Quantity that may need conversion to spectral_axis units
                    q = [pair[0].value, pair[1].value] * pair[0].unit
                    if q.unit.is_equivalent("km/s"):
                        veldef = p.meta.get("VELDEF", None)
                        if veldef is None:
                            raise KeyError("Input spectrum has no VELDEF in header, can't convert to frequency units.")
                        pair = veltofreq(q, p.rest_value, veldef)
                        # offset = p.rest_value - p.radial_velocity.to(sa.unit, equivalencies=p.equivalencies)
                    else:
                        pair[0] = pair[0].to(sa.unit, equivalencies=p.equivalencies)
                        pair[1] = pair[1].to(sa.unit, equivalencies=p.equivalencies)
                    # Ensure test is with sorted [lower,upper]
                    pair = sorted(pair)
                    salimits = sorted([sa[0], sa[-1]])
                    if pair[0] < salimits[0] or pair[1] > salimits[-1]:
                        msg = (
                            f"Exclude limits {pair} are not fully within the spectral axis"
                            f" {[salimits[0],salimits[-1]]}."
                        )
                        if fix_exclude:
                            msg += f" Setting upper limit to {p.spectral_axis[-1]}."
                            pair[1] = sa[-1]
                            warnings.warn(msg)
                        else:
                            raise Exception(msg)
                    sr = SpectralRegion(pair[0], pair[1])
                    regionlist.append(sr)

            return regionlist


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


@log_function_call()
def baseline(spectrum, order, exclude=None, exclude_region_upper_bounds=True, **kwargs):
    """Fit a baseline for a spectrum

    Parameters
    ----------
    spectrum : `~spectra.spectrum.Spectrum`
        The input spectrum
    order : int
        The order of the polynomial series, a.k.a. baseline order
    exclude : list of 2-tuples of int or `~astropy.units.Quantity`, or `~specutils.SpectralRegion`
        List of region(s) to exclude from the fit.  The tuple(s) represent a range in the form [lower,upper], inclusive.  Examples:

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
    fitter : `~astropy.fitting._FitterMeta`
        The fitter to use. Default: `~astropy.fitter.LinearLSQFitter` (with `calc_uncertaintes=True`).  Be care when choosing a different fitter to be sure it is optimized for this problem.
    exclude_region_upper_bounds : bool
        Makes the upper bound of any excision region(s) inclusive. Allows excising channel 0 for lower-sideband data, and the last channel for upper-sideband data.

    Returns
    -------
    models : list of `~astropy.modeling.Model`
        The list of models that contain the fitted model parameters.
        See `~specutuls.fitting.fit_continuum`.

    """
    kwargs_opts = {
        #'show': False,
        "model": "chebyshev",
        "fitter": LinearLSQFitter(calc_uncertainties=True),
        "fix_exclude": False,
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
    if model == None:
        raise ValueError(f'Unrecognized input model {kwargs["model"]}. Must be one of {list(available_models.keys())}')
    selected_model = available_models[model](degree=order)

    _valid_exclude_actions = ["replace", "append", None]
    if kwargs_opts["exclude_action"] not in _valid_exclude_actions:
        raise ValueError(
            f'Unrecognized exclude region action {kwargs["exclude_region"]}. Must be one of {_valid_exclude_actions}'
        )
    fitter = kwargs_opts["fitter"]
    # print(f"MODEL {model} FITTER {fitter}")
    p = spectrum
    if np.isnan(p.data).all():
        # @todo handle masks
        return None  # or raise exception
    if exclude is not None:
        regionlist = exclude_to_region(exclude, spectrum, fix_exclude=kwargs_opts["fix_exclude"])
        if kwargs_opts["exclude_action"] == "replace":
            p._exclude_regions = regionlist
        elif kwargs_opts["exclude_action"] == "append":
            p._exclude_regions.extend(regionlist)
            regionlist = p._exclude_regions
    else:
        # use the spectrum's preset exclude regions if they
        # exist (they will be a list of SpectralRegions or None)
        regionlist = p._exclude_regions
    print(f"EXCLUDING {regionlist}")
    return fit_continuum(
        spectrum=p,
        model=selected_model,
        fitter=fitter,
        exclude_regions=regionlist,
        exclude_region_upper_bounds=exclude_region_upper_bounds,
    )


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
    if type(weight) == u.Quantity:
        return weight.value.astype(np.float64)
    else:
        return weight.astype(np.float64)


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


@log_function_call()
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
    if method == None:
        raise ValueError(f"Unrecognized input method {method}. Must be one of {list(available_methods.keys())}")
    kernel = available_methods[method](width)
    if show:
        return kernel
    # the boundary='extend' matches  GBTIDL's  /edge_truncate CONVOL() method
    if hasattr(data, "mask"):
        mask = data.mask
    else:
        mask = None
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
    Shift `y` by `fshift` channels, where |`fshift`|<1.

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
