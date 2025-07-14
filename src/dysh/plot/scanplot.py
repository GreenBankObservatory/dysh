"""
Plot spectrograms from a ScanBlock using matplotlib
"""

import datetime as dt
from copy import deepcopy

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.utils.masked import Masked
from matplotlib.patches import Rectangle
from matplotlib.widgets import Button, SpanSelector

from ..coordinates import (
    Observatory,
    crval4_to_pol,
    decode_veldef,
    frame_to_label,
    ra2ha,
)

_KMS = u.km / u.s


class ScanPlot:
    """hello"""


    def __init__(self, scanblock_or_scan, **kwargs):
        self.reset()
        self._plot_kwargs.update(kwargs)
        self._figure = None
        acceptable_types = ["PSScan", "TPScan", "NodScan", "FSScan", "SubBeamNodScan"]

        #determine if input is a ScanBlock or a ScanBase (raise exception if neither)
        self._type = type(scanblock_or_scan).split(".")[-1][-2]
        if self._type == "ScanBlock":
            self._scanblock = scanblock_or_scan
            self._num_scans = len(scanblock)
        elif self._type in acceptable_types:
            self._scan = scanblock_or_scan
        else:
            raise Exception(f"Plotter input {type(scanblock_or_scan)} does not appear to be a valid input object type")

        # handle scanblocks
        if self._type == "ScanBlock":
            self._scan_nos = [] # scan numbers in the scan block
            self._nint_nos = [] # number of integrations in each scan
            self._timestamps = [] # 0-indexed timestamps in sec for every integration
            self._spectral_axis = self._scanblock[0].timeaverage().spectral_axis
            for i,scan in self._scanblock:
                if i==0:
                    self.spectrogram = scan._calibrated
                else:
                    self.spectrogram = np.append(self.spectrogram, scan._calibrated, axis=0)
                self._scan_nos.append(scan.scan)
                self._nint_nos.append(scan.nint)
                # TODO: figure out how to deal with generating a "time" axis
                # agnostic of scan proctype (pos sw, etc will have gaps between scans due to OFF)
                # self._timestamps.append(scan.)

        # handle scans
        elif self._type in acceptable_types:
            self.spectrogram = self._calibrated
            self._scan_nos = self._scan.scan





    def plot(self):
        """hi hello"""

        self.__init__(self._scanblock_or_scan, **kwargs)


        #self._set_xaxis_info()



