#import sys
#import copy
import numpy as np
from astropy.wcs import WCS
from astropy.io import fits
import astropy.units as u
#from astropy.table import Table
from astropy.modeling.polynomial import Polynomial1D,Chebyshev1D
from astropy.modeling.fitting import LevMarLSQFitter,LinearLSQFitter
from specutils import Spectrum1D, SpectrumList,SpectralRegion
from specutils.fitting import fit_continuum
import matplotlib.pyplot as plt
from ..util import uniq
import warnings

#def baseline_all(speclist,order,exclude=None,**kwargs):
#    kwargs_opts = {
#        'remove': False,
#        'show': False,
#        'model':'polynomial',
#        'fitter':  LinearLSQFitter(calc_uncertainties=True),
#    }
#    kwargs_opts.update(kwargs)
#    for p in speclist:
#        p.baseline(order,exclude,**kwargs)

def average(data,axis=0,weights=None):
    """Average a group of spectra or scans.
        TODO: allow data to be SpectrumList or array of Spectrum
       
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
    return np.average(data,axis,weights)

def exclude_to_region(exclude,refspec,fix_exclude=False):
    """Convert an exclude list to a list of ~specutuls.SpectralRegion.

       Parameters
       ----------        

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
                
            refspec: `~spectra.spectrum.Spectrum`
                The reference spectrum whose spectral axis will be used 
                when converting between exclude and axis units (e.g. channels to GHz).
            fix_exclude: bool
                If True, fix exclude regions that are out of bounds of the specctral axis to be within the spectral axis. Default:False
 
      Returns
      ----------        
            regionlist : list of `~specutil.SpectralRegion`
                A list of `~specutil.SpectralRegion` corresponding to `exclude` with units of the `refspec.spectral_axis`.

    """
    regionlist = [] 
    p = refspec
    sa = refspec.spectral_axis
    if exclude is not None:
        regionlist = [] 
        # a single SpectralRegion was given
        if isinstance(exclude,SpectralRegion):
            b = exclude.bounds
            if b[0]<sa[0] or b[1]>sa[1]:
                msg = f"Exclude limits {pair} are not fully within the spectral axis {sa}"
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
            #NB: we are assuming that a SpectralAxis is always [lower...upper].  Is this true???
            for pair in exclude:
                if type(pair[0]) == int:
                # convert channel to spectral axis units
                    lastchan = len(sa)-1
                    msg = f"Exclude limits {pair} are not fully within the spectral axis [0,{lastchan}]." 
                    if pair[0] < 0 or pair[1] > lastchan:
                        if fix_exclude:
                            msg += f" Setting upper limit to {lastchan}."
                            pair[1] = lastchan
                            warnings.warn(msg) 
                        else:
                            raise Exception(msg)
                    pair = [sa[pair[0]],sa[pair[1]]]
                # if it is already a spectral region no additional
                # work is needed
                #@TODO we should test that the SpectralRegion is not out of bounds
                if isinstance(pair[0],SpectralRegion):
                    b = pair[0].bounds
                    if b[0]<sa[0] or b[1]>sa[1]:
                        msg = f"Exclude limits {pair} are not fully within the spectral axis {p.spectral_axis}"
                        raise Exception(msg)
                    regionlist.append(pair)
                else: # it is a Quantity that may need conversion to spectral_axis units
                    if pair[0].unit.is_equivalent("km/s"):
                        offset = p.rest_value - p.radial_velocity.to(sa.unit,equivalencies = p.equivalencies)
                    else:
                        offset = 0
                    pair[0] = offset + pair[0].to(sa.unit,equivalencies = p.equivalencies)
                    pair[1] = offset + pair[1].to(sa.unit,equivalencies = p.equivalencies)
                    # Ensure test is with sorted [lower,upper]
                    pair = sorted(pair) 
                    salimits = sorted([sa[0],sa[-1]])
                    if pair[0] < salimits[0] or pair[1] > salimits[-1]:
                        msg = f"Exclude limits {pair} are not fully within the spectral axis {[salimits[0],salimits[-1]]}."
                        if fix_exclude:
                            msg += f" Setting upper limit to {p.spectral_axis[-1]}."
                            pair[1] = sa[-1]
                            warnings.warn(msg) 
                        else:
                            raise Exception(msg)
                    sr = SpectralRegion(pair[0],pair[1])
                    regionlist.append(sr)

            return regionlist

def region_to_axis_indices(region,refspec):
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
    #@TODO needs to work for multiple spectral regions? or just loop outside this call
    p = refspec
    sa = refspec.spectral_axis
    if region.lower.unit != sa.unit:
        #@todo if they are conformable, then allow it and convert
        raise Exception(f"Axis units of region [{region.lower.unit}] and refspec [{sa.unit}] not identical")
    b = [x.value for x in region.bounds]
    indices = np.abs(np.subtract.outer(sa.value, b)).argmin(0)
    return indices

def exclude_to_mask(exclude,refspec):
    # set a mask based on an exclude region
    # mask ~ exclude_to_indices(exclude_to_region())
    pass

def baseline(spectrum,order,exclude=None,**kwargs):
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
           The fitter to use. Default: `~astropy.fitter.LinearLSQFitter` (with `calc_uncertaintes=True).  Be care when choosing a different fitter to be sure it is optimized for this problem.

       Returns
       -------
       models : list of `~astropy.modeling.Model`
           The list of models that contain the fitted model parameters.
           See `~specutuls.fitting.fit_continuum`.
           
    """
    kwargs_opts = {
        #'show': False,
        'model':'polynomial',
        'fitter':  LinearLSQFitter(calc_uncertainties=True),
        'fix_exclude': False,
        'exclude_action': 'replace', # {'replace','append',None}
    }
    kwargs_opts.update(kwargs)

    _valid_models = ["polynomial", "chebyshev"]
    _valid_exclude_actions = ['replace','append',None]
    # @todo replace with minimum_string_match
    if kwargs_opts["model"] not in _valid_models:
        raise ValueError(f'Unrecognized input model {kwargs["model"]}. Must be one of {_valid_models}')
    if kwargs_opts['model'] == "polynomial":
        model = Polynomial1D(degree=order)
    elif kwargs_opts['model'] == "chebyshev":
        model = Chebyshev1D(degree=order)
    else:
        # should never get here, unless we someday allow user to input a astropy.model
        raise ValueError(f'Unrecognized input model {kwargs["model"]}. Must be one of {_valid_models}')

    if kwargs_opts['exclude_action'] not in _valid_exclude_actions:
        raise ValueError(f'Unrecognized exclude region action {kwargs["exclude_region"]}. Must be one of {_valid_exclude_actions}')
    fitter = kwargs_opts['fitter']
    #print(f"MODEL {model} FITTER {fitter}")
    p = spectrum
    if np.isnan(p.data).all():
        #@Todo handle masks
        return None # or raise exception
    if exclude is not None:
        regionlist = exclude_to_region(exclude,spectrum,fix_exclude=kwargs_opts['fix_exclude'])
        if kwargs_opts['exclude_action'] == 'replace':
            p._exclude_regions = regionlist
        elif kwargs_opts['exclude_action'] == 'append':
            p._exclude_regions.extend(regionlist)
            regionlist = p._exclude_regions
    else:
        # use the spectrum's preset exclude regions if they
        # exist (they will be a list of SpectralRegions or None)
        regionlist = p._exclude_regions
    print(f"EXCLUDING {regionlist}")
    return fit_continuum(spectrum=p,
                model=model,
                fitter=fitter,
                exclude_regions=regionlist)

def mean_tsys(calon, caloff, tcal, mode=0, fedge=10, nedge=None):
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
            TODO: Ask PJT why the options?

        fedge : int
            Fraction of edge channels to exclude at each end, in percent. Default: 10, meaning the central 80% bandwidth is used

        nedge : int
            Number of edge channels to exclude. Default: None, meaning use `fedge`

    Returns
    -------
        meanTsys : `~numpy.ndarray`-like 
            The mean system temperature
    """
    #@todo Pedro thinks about a version that takes a spectrum with multiple SpectralRegions to exclude.
    nchan = len(calon)
    if nedge == None:
        nedge = nchan // fedge    # 10 %
    # Python uses exclusive array ranges while GBTIDL uses inclusive ones.
    # Therefore we have to add a channel to the upper edge of the range
    # below in order to reproduce exactly what GBTIDL gets for Tsys.  
    # See github issue #28.
    # Define the channel range once.
    chrng = slice(nedge,-(nedge-1),1)

    # Make them doubles. Probably not worth it.
    caloff = caloff.astype('d')
    calon = calon.astype('d')

    if mode == 0:  #mode = 0 matches GBTIDL output for Tsys values
        meanoff = np.nanmean(caloff[chrng])
        meandiff = np.nanmean(calon[chrng]) - np.nanmean(caloff[chrng])
        if False:
            if meandiff < 0  :
                print(f"moff {meanoff}, mdif {meandiff}, tc {tcal}")
                print(f"CALON: {calon[nedge:-(nedge-1)]}")
                print(f"CALOF: {caloff[nedge:-(nedge-1)]}")
                print(f"DIFF: {calon[nedge:-(nedge-1)]-caloff[nedge:-(nedge-1)]}")
                print(f"CALOF: {caloff[nedge:-(nedge-1)]}")
        meanTsys = ( meanoff / meandiff * tcal + tcal/2.0 )
    else:
        meanTsys = np.mean( caloff[chrng] / (calon[chrng] - caloff[chrng]) )
        meanTsys = meanTsys * tcal + tcal/2.0

    # meandiff can sometimes be negative, which makes Tsys negative!
    # GBTIDL also takes abs(Tsys) because it does sqrt(Tsys^2)
    return np.abs(meanTsys) 


def veldef_to_convention(veldef):
    """given a VELDEF, return the velocity convention expected by Spectrum(1D)

        Parameters
        ----------
            veldef : str
                velocity definition from FITS header, e.g., 'OPTI-HELO', 'VELO-LSR'
        
        Returns
        -------
            convention : str
            velocity convention string, one of {'radio', 'optical', 'relativistic'}  or None if `velframe` can't be parsed
    """

    #@TODO GBT defines these wrong.  Need to sort out and have special version for GBT
    prefix = veldef[0:4].lower()
    if prefix == "opti":
        return 'optical'
    if prefix == "velo" or prefix == "radi":
        return 'radio'
    if prefix == "rela":
        return 'relativistic'
    return None


def tsys_weight(exposure,delta_freq,tsys):
    r"""Compute the system temperature based weight(s) using the expression:
        
        :math:`w = t_{exp} \times \delta_\nu / T_{sys}^2,`

       where :math:`t_{exp}` is the exposure time, :math:`\delta_\nu` is the frequency resolution, and :math:`T_{sys}` is the system temperature.

       If `exposure`, `delta_freq`, or `tsys` parameters are given as `~astropy.units.Quantity`,
       they will be converted to seconds, Hz, K, respectively for the calculation.
       (Not that this really matters since weights are relative to each other)
    
       **Note:** The parameters cannot be given as arrays of `~astropy.units.Quantity`, e.g., 

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
    weight = abs(delta_freq)*exposure*np.power(tsys,-2)
    if type(weight) == u.Quantity:
        return weight.value.astype(np.longdouble)
    else:
        return weight.astype(np.longdouble)
