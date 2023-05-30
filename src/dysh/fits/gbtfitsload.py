"""Load SDFITS files produced by the Green Bank Telescope"""
import sys
import copy
from astropy.wcs import WCS
from astropy.units import cds
from astropy.io import fits
from astropy.modeling import models, fitting
import astropy.units as u
from astropy.table import Table
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from ..spectra.spectrum import Spectrum
from ..spectra.obsblock import Obsblock
from ..spectra import dcmeantsys
from .sdfitsload import SDFITSLoad
from ..util import uniq

# from GBT IDL users guide Table 6.7
_PROCEDURES = ["Track", "OnOff", "OffOn", "OffOnSameHA", "Nod", "SubBeamNod"] 

class GBTFITSLoad(SDFITSLoad):
    """GBT-specific container for bintables from selected HDU(s)"""
    def __init__(self, filename, source=None,hdu=None):
        SDFITSLoad.__init__(self,filename,source,hdu)#,fix=False)
        print("==GBTLoad %s" % filename)

        self._compute_proc()
        self.ushow(0,'OBJECT')
        self.ushow(0,'SCAN')
        self.ushow(0,'SAMPLER')
        #ushow('PLNUM')
        #ushow('IFNUM')
        self.ushow(0,'SIG')
        self.ushow(0,'CAL')
        self.ushow(0,'PROCSEQN')
        self.ushow(0,'PROCSIZE')
        self.ushow(0,'OBSMODE')  
        self.ushow(0,'SIDEBAND')

    def _compute_proc(self):
        """Compute the procedure string from obsmode and add to index"""
        for i in range(len(self._ptable)):
            df = self._ptable[i]["OBSMODE"].str.split(':',expand=True)
            self._ptable[i]["PROC"] = df[0]
            # Assign these to something that might be usefule later, since we have them
            self._ptable[i]["_OBSTYPE"] = df[1]
            self._ptable[i]["_SUBOBSMODE"] = df[2]

    def summary(self, scans=None, verbose=True):
#  From GBTIDL 
#Intended to work with un-calibrated GBT data and is
# likely to give confusing results for other data.  For other data,
# list is usually more useful.
        """Create a summary list of the input dataset.  

        Parameters
        ----------
            scans : 2-tuple
                The beginning and ending scan  to use. Default: show all scans

            verbose: bool
                If True, list every record, otherwise list compact summary
                TODO : make compact work ala gbtidl
            concat : bool
                If true concatenate summaries of multiple HDUs, otherwise return list of summaries.

        Returns
        -------
            summary - list of `~pandas.DataFrame`
                Summary of the data as a DataFrame, one per HDU
            
        """
        show = ["SCAN", "OBJECT", "VELOCITY", "PROC", "PROCSEQN", 
                "RESTFREQ", "IFNUM","FEED", "AZIMUTH", "ELEVATIO", "FDNUM"] 
        summary = []
        #alternative use df.concatenate and make one big df for all HDUs
        if self._ptable is None:
            self._create_index()
        for df in self._ptable:
            _df = df[df.columns & show]
            # Don't use /= here because SettingWithCopyWarning: 
            #    A value is trying to be set on a copy of a slice from a DataFrame.
            #    Try using .loc[row_indexer,col_indexer] = value instead
            _df.loc[:,"VELOCITY"] = _df["VELOCITY"]/1E3   # convert to km/s
            _df["RESTFREQ"] = _df["RESTFREQ"]/1.0E9 # convert to GHz
            if scans is not None:
                _df = _df[(_df["SCAN"]>=scans[0]) & ( _df["SCAN"] <= scans[1])].reindex(columns=show)
                summary.append(_df)
            else:
                summary.append(_df.reindex(columns=show))

        return summary

