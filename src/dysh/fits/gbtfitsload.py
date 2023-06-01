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
            # make a copy here because we can't guarantee if this is a 
            # view or a copy without it. See https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#returning-a-view-versus-a-copy
            _df = df[df.columns & show].copy()
            _df.loc[:,"VELOCITY"] /= 1E3   # convert to km/s
            _df["RESTFREQ"] = _df["RESTFREQ"]/1.0E9 # convert to GHz
            if scans is not None:
                _df = self.select_scans(scans,_df).reindex(columns=show)
                #_df = _df[(_df["SCAN"]>=scans[0]) & ( _df["SCAN"] <= scans[1])].reindex(columns=show)
                summary.append(_df)
            else:
                summary.append(_df.reindex(columns=show))

        return summary

    def velocity_convention(self,veldef,velframe):
        # GBT uses VELDEF and VELFRAME incorrectly. 
        return "doppler_radio"

    def select_scans(self,scans,df):
        return df[(df["SCAN"]>=scans[0]) & ( df["SCAN"] <= scans[1])]

    def select_onoff(self,df):
        return df[(df["PROC"]=="OnOff") | ( df["PROC"] == "OffOn")]

    def select(self,key,value,df):
        return df[(df[key]==value) | ( df[key] == value)]

    def getps(self,scans=None,**kwargs):
        '''Get the rows that contain position-switched data.  These include ONs and OFFs.

        Parameters
        ----------
            scans : 2-tuple
                The beginning and ending scan to use. Default: show all scans.
        '''
        if self._ptable is None:
            self._create_index()
        #for df in self._ptable:

    def sonoff2(self):
        s = {"ON": [], "OFF" :[]}
        for df in self._ptable:
            #OnOff lowest scan number is on
            #dfonoff = self.select(df,"PROC","OnOff"))
            #OffOn lowest scan number is off
            #dfoffon = self.select(df,"PROC","OffOn"))
            dfon  = self.select("_OBSTYPE","PSWITCHON",df)
            dfoff = self.select("_OBSTYPE","PSWITCHOFF",df)
            onscans = uniq(list(dfon["SCAN"]))
            #print("ON: ",onscans)
            s["ON"].extend(onscans)
            offscans = uniq(list(dfoff["SCAN"]))
            #print("OFF: ",offscans)
            s["OFF"].extend(offscans)

        s["ON"] = uniq(s["ON"])
        s["OFF"] = uniq(s["OFF"])
        return s

    def sonoffrows(self,scans=None,bintable=0): #@TODO deal with mulitple bintables
        rows = {"ON": [], "OFF" :[]}
        if not scans:
            scans = self.sonoff2()
        df = self._ptable[bintable]
        for k in scans:
            rows[k] = list(df[df["SCAN"].isin(scans[k])].index)
        return rows
        
    def getscans(self,scans,bintable=0):
        df = self._ptable[bintable]
        df[df["SCAN"].isin(scans)]

# @This only works if the FITS FILE contains ONLY ons and offs
def sonoff(scan, procseqn):
    """
    return the list of On and Off scan numbers
    there must be a more elegant python way to do this....
    """
    sp = {}
    for (i,j) in zip(scan, procseqn):
        sp[i] = j
    
    us1 = uniq(scan)
    up1 = uniq(procseqn)
    
    sd = {}
    for i in up1:
        sd[i] = []
        
    for s in us1:
        sd[sp[s]].append(s)

    return sd

