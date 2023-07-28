import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
import copy
"""Plot a spectrum using matplotlib"""
class SpectrumPlot():
    '''
    The SpectrumPlot class is for simple plotting of a ~spectrum.Spectrum
    using matplotlib functions.   Plots attributes are modified using keywords
    (\*\*kwargs) described below SpectrumPlot will attempt to make smart default
    choices for the plot if no additional keywords are given.
    The attributes are "sticky" meaning that an attribute set via
    instantiation or by the `plot()` method will stay set until changed
    or reset using the `reset()` method.  

    Parameters
    ----------
    spectrum : `~spectra.spectrum.Spectrum`
        The spectrum to plot
    \**kwargs : dict
        Plot attribute keyword arguments, see below.

    Keyword Arguments
    -----------------
        xaxis_unit : str or ~astrpy.unit.Unit
            The units to use on the x-axis, e.g. "km/s" to plot velocity 
        yaxis_unit : str or ~astrpy.unit.Unit
            The units to use on the y-axis
        xmin : float
            Minimum x-axis value
        xmax : float
            Maximum x-axis value
        ymin : float
            Minimum y-axis value
        ymax : float
            Maximum y-axis value
        xlabel : str
            x-axis label
        ylabel : str
            y-axis label
        grid : bool
            Show a plot grid or not
        figsize : tuple
            Figure size (see matplotlib)
        linewidth : float
            Line width, default: 2.0.  lw also works
        linestyle : str
            Line style, default 'steps-mid'.  ls also works
        color : str
            Line color, c also works
        title : str
            Plot title
        aspect : str
            plot aspect ratio, default: 'auto'
        show_baseline : bool
            show the baseline - not yet implemented
    '''
    # loc, legend, bbox_to_anchor

    def __init__(self,spectrum,**kwargs):
        self.reset()
        self._plot_kwargs.update(kwargs)
        self._plt = plt
        self._figure = None
        self._axis = None
        self._spectrum = spectrum
        self._title = self._plot_kwargs['title']

# def __call__ (see pyspeckit)

    @property
    def axis(self):
        """The underlying :class:`~matplotlib.Axes` object"""
        return self._axis
    @property
    def figure(self):
        """The underlying :class:`~matplotlib.Figure` object"""
        return self._figure
    @property
    def spectrum(self):
        """The underlying `~spectra.spectrum.Spectrum`"""
        return self._spectrum

    def plot(self,**kwargs):
        """
        Plot the spectrum.

        Parameters
        ----------
        \**kwargs : various
            keyword=value arguments (need to describe these in a central place)
        """
        # xtype = 'velocity, 'frequency', 'wavelength'
        self._plot_kwargs.update(kwargs)
        #if self._figure is None:
        if True:
            self._figure,self._axis = self._plt.subplots(figsize=self._plot_kwargs['figsize'])
        #else:
        #    self._axis.cla()

        s = self._spectrum
        sa = s.spectral_axis
        lw =  self._plot_kwargs['linewidth']
        xunit = self._plot_kwargs["xaxis_unit"]
        yunit = self._plot_kwargs["yaxis_unit"]
        if xunit is not None:
            if "chan" in xunit:
                sa = np.arange(len(sa))
                self._plot_kwargs['xlabel'] = "Channel"
            else:
                # convert the x axis to the requested 
                #print(f"EQUIV {equiv} doppler_rest {sa.doppler_rest} [{rfq}] convention {convention}")
                #sa = s.spectral_axis.to( self._plot_kwargs["xaxis_unit"], equivalencies=equiv,doppler_rest=rfq, doppler_convention=convention)
                sa = s.velocity.to( self._plot_kwargs["xaxis_unit"], equivalencies=s.equivalencies)
                self._plot_kwargs['xlabel'] = f'Velocity ({xunit})'
        sf = s.flux
        if yunit is not None:
            sf = s.flux.to(yunit)
        self._axis.plot(sa,sf,color=self._plot_kwargs['color'], lw = lw)
        self._axis.set_xlim(self._plot_kwargs['xmin'],self._plot_kwargs['xmax'])
        self._axis.set_ylim(self._plot_kwargs['ymin'],self._plot_kwargs['ymax'])
        self._axis.tick_params(axis='both',which='both',bottom=True,top=True,left=True,right=True,direction='in')
        if self._plot_kwargs['grid']:
            self._axis.grid(visible=True,which='major',axis='both',lw=lw/2,
                                color='k',alpha=0.33)
            self._axis.grid(visible=True,which='minor',axis='both',lw=lw/2,
                                color='k',alpha=0.22,linestyle='--')

        self._set_labels(**self._plot_kwargs)
        #self._axis.axhline(y=0,color='red',lw=2)
        self.refresh()

    def reset(self):
        """Reset the plot keyword arguments to their defaults."""
        self._plot_kwargs = {
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

    def _set_labels(self,title=None,xlabel=None,ylabel=None,**kwargs):
        """Set x and y labels according to spectral units

           Parameters
           ----------
           title : str
               plot title
           xlabel : str
               x-axis label
           ylabel : str
               x-axis label
           \*\*kwargs : various
               other keyword=value arguments
        """
        if title is not None:
            self._title = title
        if hasattr(self.spectrum.wcs,'wcs'):
           ctype = self.spectrum.wcs.wcs.ctype
        elif self.spectrum.meta is not None:
           ctype = []
           ctype.append(self.spectrum.meta.get('CTYPE1',None))
           ctype.append(self.spectrum.meta.get('CTYPE2',None))
           ctype.append(self.spectrum.meta.get('CTYPE3',None))
        #print('ctype is ',ctype)
        if kwargs.get("xaxis_unit",None) is not None:
            xunit = kwargs["xaxis_unit"]
        else:
            xunit = self.spectrum.spectral_axis.unit
        if kwargs.get("yaxis_unit",None) is not None:
            yunit = u.Unit(kwargs["yaxis_unit"])
        else:
            yunit = self.spectrum.unit
        if xlabel is not None:
            self.axis.set_xlabel(xlabel)
        elif ctype[0] in ['FREQ']:
            xlabel = f'Frequency ({xunit})'
            self.axis.set_xlabel(xlabel)
        elif ctype[0] in ['VELO', 'VRAD', 'VOPT']:
            xlabel = f'Velocity ({xunit})'
            self.axis.set_xlabel(xlabel)
        elif ctype[0] in ['WAVE', 'AWAV']:
            xlabel = f'Wavelength({xunit})'
            self.axis.set_xlabel(xlabel)
        elif xunit is not None:
            xlabel = xunit
            self.axis.set_xlabel(xlabel)
        #print(f"ylabel {ylabel} yunit {yunit} sunit {self.spectrum.unit}")
        if ylabel is not None:
            self.axis.set_ylabel(ylabel)
        elif yunit.is_equivalent(u.K):
            self.axis.set_ylabel(f"$T_A$ ({yunit})")
        elif self.spectrum.unit.is_equivalent(u.Jy):
            snu = r"$S_{\nu}$"
            self.axis.set_ylabel(f"{snu} ({yunit})")

    def _show_exclude(self,**kwargs):
        '''TODO: Method to show the exclude array on the plot'''
        kwargs_opts = {
            'loc': "bottom", #top,bottom ?
            'color' : 'silver',
        }
        kwargs_opts.update(kwargs)
        #if kwargs_opts['loc'] == 'bottom':
        #    self._ax.axhline
        


    def refresh(self):
        """Refresh the plot"""
        if self.axis is not None:
            self.axis.figure.canvas.draw()
            #print('redrawing')
            #self.axis.figure.canvas.draw_idle()
            self._plt.show()

    def savefig(self,file,**kwargs):
        """Save the plot"""
        self.figure.savefig(file,*kwargs)

