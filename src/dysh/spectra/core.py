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

def baseline_all(speclist,order,exclude=None,**kwargs):
    kwargs_opts = {
        'remove': False,
        'show': False,
        'model':'polynomial',
        'fitter':  LinearLSQFitter(calc_uncertainties=True),
    }
    kwargs_opts.update(kwargs)
    i=0
    bad = 0
    for p in speclist:
        p.baseline(order,exclude,**kwargs)

def baseline(spectrum,order,exclude=None,**kwargs):
    """Fit a baseline for a spectrum

       Parameters
       ----------        
            spectrum : ~Spectrum
                The input spectrum
            order : int
                The order of the polynomial series, a.k.a. baseline order
            exclude : list of 2-tuples of int or ~astropy.units.Quantity, or ~specutils.SpectralRegion
                List of region(s) to exclude from the fit.  The tuple(s) represent a range in the form [lower,upper], inclusive.  
in channel units.  

                Examples: One channel-based region: [11,51], Two channel-based regions: [(11,51),(99,123)]. One ~astropy.units.Quantity region: [110.198*u.GHz,110.204*u.GHz]. One compound ~specutils.SpectralRegion: SpectralRegion([(110.198*u.GHz,110.204*u.GHz),(110.196*u.GHz,110.197*u.GHz)]).

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
    }
    kwargs_opts.update(kwargs)

    _valid_models = ["polynomial", "chebyshev"]
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
    fitter = kwargs_opts['fitter']
    #print(f"MODEL {model} FITTER {fitter}")
    p = spectrum
    if np.isnan(p.data).all():
        print('returning none')
        return None # or raise exception
    if exclude is not None:
        regionlist = []
        # a single SpectralRegion was given
        if isinstance(exclude,SpectralRegion):
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
                    lastchan = len(p.spectral_axis)-1
                    msg = f"Exclude limits {pair} are not fully within the spectral axis [0,{lastchan}]." 
                    if pair[0] < 0 or pair[1] > lastchan:
                        if kwargs_opts['fix_exclude']:
                            msg += f" Setting upper limit to {lastchan}."
                            pair[1] = lastchan
                            warnings.warn(msg) 
                        else:
                            raise Exception(msg)
                    pair = [p.spectral_axis[pair[0]],p.spectral_axis[pair[1]]]
                # if it is already a spectral region no additional
                # work is needed
                #@TODO we should test that the SpectralRegion is not out of bounds
                if isinstance(pair[0],SpectralRegion):
                    regionlist.append(pair)
                else: # it is a Quantity that may need conversion to spectral_axis units
                    if pair[0].unit.is_equivalent("km/s"):
                        offset = p.rest_value - p.radial_velocity.to(p.spectral_axis.unit,equivalencies = p.equivalencies)
                    else:
                        offset = 0
                    pair[0] = offset + pair[0].to(p.spectral_axis.unit,equivalencies = p.equivalencies)
                    pair[1] = offset + pair[1].to(p.spectral_axis.unit,equivalencies = p.equivalencies)
                    pair = sorted(pair) # SpectralRegion requires sorted [lower,upper]
                    if pair[0] < p.spectral_axis[0] or pair[1] > p.spectral_axis[-1]:
                        msg = f"Exclude limits {pair} are not fully within the spectral axis {[p.spectral_axis[0],p.spectral_axis[-1]]}."
                        if kwargs_opts['fix_exclude']:
                            msg += f" Setting upper limit to {p.spectral_axis[-1]}."
                            pair[1] = p.spectral_axis[-1]
                            warnings.warn(msg) 
                        else:
                            raise Exception(msg)
                    sr = SpectralRegion(pair[0],pair[1])
                    print(f"EXCLUDING {sr}")
                    regionlist.append(sr)
        return fit_continuum(spectrum=p,
                model=model,
                fitter=fitter,
                exclude_regions=regionlist)
    else:
        return fit_continuum(spectrum=p,
                model=model,
                fitter=fitter)

def baseline_old(speclist,order,exclude=None,**kwargs):
    kwargs_opts = {
        'remove': True,
        'show': False,
        'model':'polynomial',
        'fitter':  LinearLSQFitter(calc_uncertainties=True),
    }
    kwargs_opts.update(kwargs)

    _valid_models = ["polynomial", "chebyshev"]
    # @todo replace with minimum_string_match
    if kwargs["model"] not in _valid_models:
        raise ValueError(f'Unrecognized input model {kwargs["model"]}. Must be one of {_valid_models}')
    print(f"BL {order} for {len(speclist)} spectra")
    i=0
    bad = 0
    if kwargs_opts['model'] == "polynomial":
        model = Polynomial1D(degree=order)
    elif kwargs_opts['model'] == "chebyshev":
        model = Chebyshev1D(degree=order)
    else:
        # should never get here
        raise ValueError(f'Unrecognized input model {kwargs["model"]}. Must be one of {_valid_models}')
    fitter = kwargs_opts['fitter']
    print(f"MODEL {model} FITTER {fitter}")
    try:
        if exclude is not None:
            for p in speclist:
                if np.isnan(p.data).all():
                    bad+=1
                    continue
                fc = fit_continuum(spectrum=p,
                        model=model,
                        fitter=fitter,
                        exclude_regions=[exclude])
                i=i+1
        else:
            for p in speclist:
                if np.isnan(p.data).all():
                    bad+=1
                    continue
                fc = fit_continuum(spectrum=p,
                        model=model,
                        fitter=fitter)
                i=i+1
    except Exception as e:
        print(f"At spectrum {i}, Exception was {e}")
        print(p)
        print("DATA MEAN: ",np.nanmean(p.data))
        return p
    if plot:
        fig,ax = plt.subplots()
        ax.plot(x,p.flux)
        ax.plot(x,fc(x))
        plt.show()
    print(f"NUMBER OF BAD SPECTRA: {bad}")
    return None

def dcmeantsys(calon, caloff, tcal, mode=0, fedge=10, nedge=None):
    """
    Following the GBTIDL routine with same name, get the system temperature from 
    the neighboring calon and caloff, which reflect the state of the noise diode.
    We define an extra way to set the edge size, nedge, if you prefer to use 
    number of edge channels instead of the inverse fraction.
    
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
    nchan = len(calon)
    if nedge == None:
        nedge = nchan // fedge     # 10 %
    if mode == 0:
        meanoff = np.mean(caloff[nedge:-nedge])
        meandiff = np.mean(calon[nedge:-nedge] - caloff[nedge:-nedge])
        meanTsys = ( meanoff / meandiff * tcal + tcal/2.0 )
    else:
        meanTsys = np.mean( caloff[nedge:-nedge] / (calon[nedge:-nedge] - caloff[nedge:-nedge]) )
        meanTsys = meanTsys * tcal + tcal/2.0
    return meanTsys


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
