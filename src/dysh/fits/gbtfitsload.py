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
from ..spectra.scan import GBTPSScan
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

    def summary(self, scans=None, verbose=True, bintable=0):
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

            bintable : int or list
                Which bintable to summarize, zero-based. Default:0
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

    def _create_index_if_needed(self):
        if self._ptable is None:
            self._create_index()

    def getps(self,scans=None,bintable=0,**kwargs):
        '''Get the rows that contain position-switched data.  These include ONs and OFFs.

            kwargs: plnum, feed, ifnum, integration, calibrate=T/F, average=T/F, tsys, weights
            [sampler], ap_eff [if requested units are Jy]
            
        Parameters
        ----------
            scans : int or 2-tuple
                Single scan number or list of scan numbers to use. Default: all scans.
                Scan numbers can be Ons or Offs
        TODO: figure how to allow [startscan, endscan]
        Returns 
        -------
                ? ScanBlock
        '''
        # all ON/OFF scans
        kwargs_opts = {
            'ifnum': 0,
            'plnum' : 0, # I prefer "pol"
            'fdnum' : 0,
            'calibrate': True,
            'average': False,
            'tsys': None,
            'weights': None,
        }
        kwargs_opts.update(kwargs)

        ifnum = kwargs_opts['ifnum']
        plnum = kwargs_opts['plnum']
        scanlist = self.onoff_scan_list(scans,ifnum=ifnum,plnum=plnum,bintable=bintable)
        # add ifnum,plnum
        rows = self.onoff_rows(scans,ifnum=ifnum,plnum=plnum,bintable=bintable)
        # do not pass scan list here. We need all the cal rows. They will 
        # be intersected with scan rows in GBTPSScan
        # add ifnum,plnum
        calrows = self.calonoff_rows(scans=None,bintable=bintable)
        g = GBTPSScan(self,scanlist,rows,calrows,bintable)
        return g


    def onoff_scan_list(self,scans=None,ifnum=0,plnum=0,bintable=0):
        self._create_index_if_needed()
        #print(f"onoff_scan_list(scans={scans},if={ifnum},pl={plnum})")
        s = {"ON": [], "OFF" :[]}
        if type(scans) == int:
            scans = [scans]
        if False: # this does all bintables.
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

        df    = self._ptable[bintable]
        df = df[(df["PLNUM"] == plnum) & (df["IFNUM"] == ifnum)]
        dfon  = self.select("_OBSTYPE","PSWITCHON",df)
        dfoff = self.select("_OBSTYPE","PSWITCHOFF",df)
        onscans = uniq(list(dfon["SCAN"])) # wouldn't set() do this too?
        offscans = uniq(list(dfoff["SCAN"]))
        if scans is not None:
        # The companion scan will always be +/- 1 depending if procseqn is 1(ON) or 2(OFF)
        # First check the requested scan number(s) are even in the ONs or OFFs of this bintable
            seton = set(onscans)
            setoff = set(offscans)
            onrequested = seton.intersection(scans)
            #print("ON REQUESTED ",onrequested)
            offrequested = setoff.intersection(scans)
            #print("OFF REQUESTED ",offrequested)
            if len(onrequested) == 0 and len(offrequested) == 0:
                raise ValueError(f"Scans {scans} not found in ONs or OFFs of bintable {bintable}")
        # Then check that for each requested ON/OFF there is a matching OFF/ON
        # and build the final matched list of ONs and OFfs.
            sons = list(onrequested.copy())
            soffs = list(offrequested.copy())
            missingoff = []
            missingon = []
            for i in onrequested:
                expectedoff = i+1
                #print(f"DOING ONQUESTED {i}, looking for off {expectedoff}")
                if len(setoff.intersection([expectedoff])) == 0:
                    missingoff.append(expectedoff)
                else:
                    soffs.append(expectedoff)
            for i in offrequested:
                expectedon = i-1
                #print(f"DOING OFFEQUESTED {i}, looking for on {expectedon}")
                if len(seton.intersection([expectedon])) == 0:
                    missingon.append(expectedon)
                else:
                    sons.append(expectedon)
            if len(missingoff) > 0:
                raise ValueError(f"For the requested ON scans {onrequested}, the OFF scans {missingoff} were not present in bintable {bintable}")
            if len(missingon) > 0:
                raise ValueError(f"For the requested OFF scans {offrequested}, the ON scans {missingon} were not present in bintable {bintable}")
            #print("ON",sorted(sons))
            #print("OFF",sorted(soffs))
            s["ON"] = sorted(set(sons))
            s["OFF"] = sorted(set(soffs))
        else:
            s["ON"] = uniq(list(dfon["SCAN"]))
            s["OFF"] = uniq(list(dfoff["SCAN"]))

        return s

    def calonoff_rows(self,scans=None,bintable=0):
        self._create_index_if_needed()
        s = {"ON": [], "OFF" :[]}
        if type(scans) == int:
            scans = [scans]
        df    = self._ptable[bintable]
        if scans is not None:
            df = df[df["SCAN"].isin(scans)]
        dfon  = self.select("CAL","T",df)
        dfoff = self.select("CAL","F",df)
        s["ON"]  = list(dfon.index)
        s["OFF"] = list(dfoff.index)
        return s

    def onoff_rows(self,scans=None,ifnum=0,plnum=0,bintable=0): 
    #@TODO deal with mulitple bintables
    #@TODO rename this sigref_rows?
    # keep the bintable keyword and allow iteration over bintables if requested (bintable=None) 
        #print(f"onoff_rows(scans={scans},if={ifnum},pl={plnum})")
        self._create_index_if_needed()
        rows = {"ON": [], "OFF" :[]}
        if type(scans) is int:
            scans = [scans]
        if scans is not None:
            scans = self.onoff_scan_list(None,ifnum,plnum,bintable)
        else:
            scans = self.onoff_scan_list(scans,ifnum,plnum,bintable)
        #scans is now a dict of "ON" "OFF
        for key in scans: 
            rows[key] = self.scan_rows(scans[key],ifnum,plnum,bintable)
        return rows
        
    def scan_rows(self,scans,ifnum=0,plnum=0,bintable=0):
        #scans is a list
        #print(f"scan_rows(scans={scans},if={ifnum},pl={plnum})")
        self._create_index_if_needed()
        if scans is None:
            raise ValueError("Parameter 'scans' cannot be None. It must be int or list of int")
        df = self._ptable[bintable]
        df = df[df["SCAN"].isin(scans) & (df["IFNUM"] == ifnum) & (df["PLNUM"] == plnum)]
        rows = list(df.index)
        if len(rows) == 0:
            raise Exception(f"Scans {scans} not found in bintable {bintable}")
        return rows
