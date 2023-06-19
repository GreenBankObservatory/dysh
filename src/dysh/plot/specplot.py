from ..spectra.spectrum import Spectrum
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
import copy

class SpectrumPlot():
    def __init__(self,spectrum,**kwargs):
        kwargs_opts = {
            'xmin': None,
            'xmax': None,
            'ymin': None,
            'ymax': None,
            'xlabel':None,
            'ylabel':None,
            'grid' :False,
            'figsize':None,
            #'capsize':3,
            'linewidth': 2.0,
            'linestyle': 'steps-mid',
            'markersize': 8,
            'color':None,
            #'axis':None,
            #'label':None,
            'aspect': 'auto',
            'bbox_to_anchor':None,
            'loc':"best",
            'legend':None,
            'test': False
        }
        kwargs_opts.update(kwargs)
        self._plot_kwargs = kwargs_opts
        self._plt = plt
        self._figure = None
        self._axis = None
        self._spectrum = spectrum

# def __call__ (see pyspeckit)

    @property
    def axis(self):
        return self._axis
    @property
    def figure(self):
        return self._figure
    @property
    def spectrum(self):
        return self._spectrum

    def plot(self,**kwargs):
        """Plot the spectrum"""
        self._plot_kwargs.update(kwargs)
        if self._figure is None:
            self._figure,self._axis = self._plt.subplots(figsize=self._plot_kwargs['figsize'])

        s = self._spectrum
        lw =  self._plot_kwargs['linewidth']
        self._axis.plot(s.spectral_axis,s.flux,color=self._plot_kwargs['color'],
                        lw = lw)
        self._axis.set_xlim(self._plot_kwargs['xmin'],self._plot_kwargs['xmax'])
        self._axis.set_ylim(self._plot_kwargs['ymin'],self._plot_kwargs['ymax'])
        if self._plot_kwargs['grid']:
            self._axis.grid(visible=True,which='major',axis='both',lw=lw/2,
                                color='k',alpha=0.33)
            self._axis.grid(visible=True,which='minor',axis='both',lw=lw/2,
                                color='k',alpha=0.22,linestyle='--')

        self._set_labels()

    def _set_labels(self):
        """set x and y labels according to spectral units"""
        pass


