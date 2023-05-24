import sys
import copy
import numpy as np
from astropy.wcs import WCS
from astropy.units import cds
from astropy.modeling.models import Gaussian1D
from astropy.modeling.polynomial import Polynomial1D
from astropy.io import fits
from astropy.modeling import models, fitting
import astropy.units as u
from astropy.table import Table
from astropy.modeling.fitting import LevMarLSQFitter,LinearLSQFitter
from specutils import Spectrum1D, SpectrumList,SpectralRegion
from specutils.fitting import fit_continuum
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from .spectrum import Spectrum

class Obsblock():
    '''Class that holds a series of spectra on which bulk operations can be performed'''
    def __init__(self,speclist,index):
        self._speclist = speclist
        self._index = index # pandas dataframe

    def __getitem__(self,i):
        return self._speclist[i]

    def __len__(self):
        return len(self._speclist)

    def __op__(self,opname):
        pass

    def baseline(self,order,exclude=None,**kwargs):
        """compute and optionally remove a baseline"""
        pass

