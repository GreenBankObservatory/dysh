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
            exclude : list of 2-tuples
                List of regions to exclude from the fit, in channel units.  The tuple(s) are in the form [lower,upper], inclusive.  Examples: One region: [11,51], Two regions: [(11,51),(99,123)]. Default: no exclude region
            model : str
                One of 'polynomial' or 'chebyshev', Default: 'polynomial'
            fitter  :  `~astropy.fitting._FitterMeta`
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
        unit = p._spectral_axis.unit
        if type(exclude[0] == int): #ugh, klugy
            exclude = [exclude]
        for pair in exclude:
            regionlist.append(SpectralRegion(pair[0]*unit,pair[1]*unit))
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

