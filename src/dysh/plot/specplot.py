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
            'xaxis_unit': None,
            'yaxis_unit': None,
            'grid' :False,
            'figsize':None,
            #'capsize':3,
            'linewidth': 2.0,
            'linestyle': 'steps-mid',
            'markersize': 8,
            'color':None,
            'title':None,
            #'axis':None,
            #'label':None,
            'aspect': 'auto',
            'bbox_to_anchor':None,
            'loc':"best",
            'legend':None,
            'show_baseline':True,
            'test': False
        }
        kwargs_opts.update(kwargs)
        self._plot_kwargs = kwargs_opts
        self._plt = plt
        self._figure = None
        self._axis = None
        self._spectrum = spectrum
        self._title = self._plot_kwargs['title']

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
        # xtype = 'velocity, 'frequency', 'wavelength'
        self._plot_kwargs.update(kwargs)
        if self._figure is None:
            self._figure,self._axis = self._plt.subplots(figsize=self._plot_kwargs['figsize'])
        #else:
        #    self._axis.cla()

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

        self._set_labels(self._plot_kwargs['title'],self._plot_kwargs['xlabel'],self._plot_kwargs['ylabel'])
        self.refresh()

    def _set_labels(self,title=None,xlabel=None,ylabel=None,**kwargs):
        """set x and y labels according to spectral units"""
        if title is not None:
            self._title = title
        if hasattr(self.spectrum.wcs,'wcs'):
           ctype = self.spectrum.wcs.wcs.ctype
        elif self.spectrum.meta is not None:
           ctype = []
           ctype.append(self.spectrum.meta.get('CTYPE1',None))
           ctype.append(self.spectrum.meta.get('CTYPE2',None))
           ctype.append(self.spectrum.meta.get('CTYPE3',None))
        print('ctype is ',ctype)
        xunit = self.spectrum.spectral_axis.unit
        yunit = self.spectrum.unit
        if xlabel is not None:
            self.axis.set_xlabel(xlabel)
        elif ctype[0] in ['FREQ']:
            xlabel = f'Frequency ({xunit})'
            self.axis.set_xlabel(xlabel)
        if ylabel is not None:
            self.axis.set_ylabel(ylabel)
        elif yunit in ['K']:
            self.axis.set_ylabel(r"Brightness Temperature $T$ (K)")
        elif self.Spectrum.unit == 'mJy':
            self.axis.set_ylabel(r"$S_\\nu$ (mJy)")
        elif self.Spectrum.unit == 'Jy':
            self.axis.set_ylabel(r"$S_\\nu$ (Jy)")


    def refresh(self):
        if self.axis is not None:
            self.axis.figure.canvas.draw()
            print('redrawing')
            #self.axis.figure.canvas.draw_idle()
            self._plt.show()


