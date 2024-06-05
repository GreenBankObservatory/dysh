import matplotlib.pyplot as plt
#from dysh.fits.gbtfitsload import GBTFITSLoad

class SpectrumPlot:
    """Single plot"""

    def __init__(self):
        self._plt = plt
        self._fig = None
        self._axes = None
        self._spectrum = None
        self.reset_kwargs()
        self.make_canvas()
        self.plot_spectrum()

    @property
    def fig(self):
        """The underlying :class:`~matplotlib.Axes` object"""
        return self._fig
    
    @property
    def axes(self):
        """The underlying :class:`~matplotlib.Axes` object"""
        return self._axes
    
    @property
    def spectrum(self):
        """The underlying `~spectra.spectrum.Spectrum`"""
        return self._spectrum
    
    def set_spectrum(self, spectrum):
        self.spectrum = spectrum

    def make_canvas(self):
        self._fig, self._axes = plt.subplots(1)
        self.update_kwargs()

    def plot_spectrum(self):
        # s = self._spectrum
        # sa = s.spectral_axes
        # lw = self._plot_kwargs["linewidth"]
        # xunit = self._plot_kwargs["xaxis_unit"]
        # yunit = self._plot_kwargs["yaxis_unit"]
        self._plt.plot([1,2,3,4], [1,2,3,4], color=self._axes_kwargs["line_color"])

    def update_kwargs(self):
        """Update the axis configuration"""
        self._fig.set_facecolor(self._axes_kwargs["facecolor"])
        self._axes.set_facecolor(self._axes_kwargs["facecolor"])

        self._axes.title.set_color(self._axes_kwargs["title_color"])
        self._axes.xaxis.label.set_color(self._axes_kwargs["xlabel_color"])
        self._axes.yaxis.label.set_color(self._axes_kwargs["ylabel_color"])

        self._axes.tick_params(axis='x', colors=self._axes_kwargs["tick_color"])
        self._axes.tick_params(axis='y', colors=self._axes_kwargs["tick_color"])

        self._axes.spines['left'].set_color(self._axes_kwargs["spine_color"])
        self._axes.spines['right'].set_color(self._axes_kwargs["spine_color"])
        self._axes.spines['top'].set_color(self._axes_kwargs["spine_color"])
        self._axes.spines['bottom'].set_color(self._axes_kwargs["spine_color"])

    def reset_kwargs(self):
        """Reset the plot keyword arguments to their defaults."""
        self._axes_kwargs = {
            "xmin": None,
            "xmax": None,
            "ymin": None,
            "ymax": None,
            "xlabel": None,
            "ylabel": None,
            "xlabel_color": "white",
            "ylabel_color": "white",
            "xaxis_unit": None,
            "yaxis_unit": None,
            "tick_color": "white",
            "spine_color": "white",
            
            "linewidth": 2.0,
            "linestyle": "steps-mid",
            "markersize": 8,
            "line_color": "red",
            "title": None,
            "title_color": "white",

            "aspect": "auto",
            "bbox_to_anchor": None,
            "loc": "best",
            "legend": None,
            "show_baseline": True,
            "test": False,

            "titlecolor": "white",
            "facecolor": "black",
            "edgecolor": "none",
        }
    
    def show(self):
        """Show the plot"""
        self._plt.show()
    
    def savefig(self, file, **kwargs):
        """Save the plot"""
        self.figure.savefig(file, facecolor=self.figure.get_facecolor(), edgecolor=self.figure.get_edgecolor())

class GridPlot:
    """Plot that can have multiple subplots"""

    def __init__(self, spectrum=None, **kwargs):
        self.reset_kwargs()
        self._plt = plt
        self._figure = None
        self._axes = None
        self._axes_mosaic = None
        self._axes_dict = None
        self.make_canvas()
        # self._title = self._plot_kwargs["title"]

    @property
    def axis(self):
        """The underlying :class:`~matplotlib.Axes` object"""
        return self._axes

    @property
    def figure(self):
        """The underlying :class:`~matplotlib.Figure` object"""
        return self._figure
    
    def make_canvas(self, **kwargs):
        """Make a blank figure"""
        self._figure = self._plt.figure()

        self.reset_mosaic()
        self.set_mosaic(self._axes_mosaic)
        #self.make_canvas()
        self.update_kwargs()

    def reset_mosaic(self):
        """Reset the mosaic to a single plot"""
        self._axes_mosaic = """
            A
            """
        
    def identify_axes(ax_dict, fontsize=48):
        """
        Helper to identify the Axes in the examples below.

        Draws the label in a large font in the center of the Axes.

        From https://matplotlib.org/stable/users/explain/axes/mosaic.html#mosaic

        Parameters
        ----------
        ax_dict : dict[str, Axes]
            Mapping between the title / label and the Axes.
        fontsize : int, optional
            How big the label should be.
        """
        kw = dict(ha="center", va="center", fontsize=fontsize, color="darkgrey")
        for k, ax in ax_dict.items():
            ax.text(0.5, 0.5, k, transform=ax.transAxes, **kw)

    def set_mosaic(self, mosaic):
        """Set the figure mosaic"""
        self._axes_mosaic = mosaic
        self._axes_dict = self._figure.subplot_mosaic(mosaic)

    def update_kwargs(self):
        """Update the figure configuration"""
        self._figure.set_facecolor(self._plot_kwargs["facecolor"])
        self._figure.set_edgecolor(self._plot_kwargs["edgecolor"])

    def set_title(self, title: str):
        """Add or update the title of the figure"""
        self._figure.suptitle(title, color=self._plot_kwargs["titlecolor"])

    def reset_kwargs(self):
        """Reset the plot keyword arguments to their defaults."""
        self._plot_kwargs = {
            "xmin": None,
            "xmax": None,
            "ymin": None,
            "ymax": None,
            "xlabel": None,
            "ylabel": None,
            "xaxis_unit": None,
            "yaxis_unit": None,
            "grid": False,
            "figsize": None,
            #'capsize':3,
            "linewidth": 2.0,
            "linestyle": "steps-mid",
            "markersize": 8,
            "color": None,
            "title": None,
            #'axis':None,
            #'label':None,
            "aspect": "auto",
            "bbox_to_anchor": None,
            "loc": "best",
            "legend": None,
            "show_baseline": True,
            "test": False,

            "titlecolor": "white",
            "facecolor": "black",
            "edgecolor": "none",
        }

    def show(self):
        """Show the plot"""
        self._figure.show()
    
    def savefig(self, file, **kwargs):
        """Save the plot"""
        self.figure.savefig(file, facecolor=self.figure.get_facecolor(), edgecolor=self.figure.get_edgecolor()) # *kwargs)