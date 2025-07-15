"""
Plot spectrograms from a ScanBlock using matplotlib
"""

import datetime as dt
from copy import deepcopy

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np




_KMS = u.km / u.s


class ScanPlot:
    """hello"""


    def __init__(self, scanblock_or_scan, **kwargs):
        self.reset()
        self._scanblock_or_scan = scanblock_or_scan
        self._plot_kwargs.update(kwargs)
        self._plt = plt
        self._figure = None
        self._axis = None
        self._title = self._plot_kwargs["title"]
        acceptable_types = ["PSScan", "TPScan", "NodScan", "FSScan", "SubBeamNodScan"]

        #determine if input is a ScanBlock or a ScanBase (raise exception if neither)
        self._type = str(type(scanblock_or_scan)).split(".")[-1][:-2]
        if self._type == "ScanBlock":
            self._scanblock = scanblock_or_scan
            self._num_scans = len(self._scanblock)
        elif self._type in acceptable_types:
            self._scan = scanblock_or_scan
        else:
            raise Exception(f"Plotter input {self._type} does not appear to be a valid input object type")

        # handle scanblocks
        if self._type == "ScanBlock":
            self._scan_nos = [] # scan numbers in the scan block
            self._nint_nos = [] # number of integrations in each scan
            self._timestamps = [] # 0-indexed timestamps in sec for every integration
            self._spectral_axis = self._scanblock[0].timeaverage().spectral_axis
            for i,scan in enumerate(self._scanblock):
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





    def plot(self, **kwargs):
        """hi hello"""

        self.__init__(self._scanblock_or_scan, **kwargs)
        plt.ion()

        #self._set_xaxis_info()
        this_plot_kwargs = deepcopy(self._plot_kwargs)
        this_plot_kwargs.update(kwargs)

        if True:
            self._figure, self._axis = self._plt.subplots(figsize=(10,6))
        
        #ax_lw = 3
        #self._axis.tick_params(axis='both',direction='in',width=2,length=8,top=True,right=True,pad=2)
        print(self.spectrogram)
        self._axis.imshow(self.spectrogram)
    




    def reset(self):
        self._plot_kwargs = {
            "title": None,
        }



