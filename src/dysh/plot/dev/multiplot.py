import matplotlib.pyplot as plt
from matplotlib.axes import Axes

# from patchworklib import Brick


class SinglePlot:
    """Single plot"""

    def __init__(self):
        self._plt = plt
        self.reset_kwargs()
        self.make_canvas()
        self._plt.close()

    @property
    def fig(self):
        """The underlying :class:`~matplotlib.Axes` object"""
        return self._fig

    @property
    def axis(self):
        """The underlying :class:`~matplotlib.Axes` object"""
        return self._axis

    def make_canvas(self):
        self._fig, self._axis = plt.subplots(1)
        self.update_kwargs()

    def update_kwargs(self):
        """Update the axis configuration"""
        self._axis.set_facecolor(self._axis_kwargs["facecolor"])

        self._axis.title.set_color(self._axis_kwargs["title_color"])
        self._axis.xaxis.label.set_color(self._axis_kwargs["xlabel_color"])
        self._axis.yaxis.label.set_color(self._axis_kwargs["ylabel_color"])

        self._axis.tick_params(axis="x", colors=self._axis_kwargs["tick_color"])
        self._axis.tick_params(axis="y", colors=self._axis_kwargs["tick_color"])

        self._axis.spines["left"].set_color(self._axis_kwargs["spine_color"])
        self._axis.spines["right"].set_color(self._axis_kwargs["spine_color"])
        self._axis.spines["top"].set_color(self._axis_kwargs["spine_color"])
        self._axis.spines["bottom"].set_color(self._axis_kwargs["spine_color"])

    def reset_kwargs(self):
        """Reset the plot keyword arguments to their defaults."""
        self._axis_kwargs = {
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
            "color": None,
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


class GridPlot:
    """Plot that can have multiple subplots"""

    def __init__(self, spectrum=None, **kwargs):
        self.reset_kwargs()
        # self._spectrum = spectrum
        # self._plot_kwargs["vel_convention"] = spectrum.velocity_convention
        # self._plot_kwargs["vel_frame"] = spectrum.velocity_frame
        # self._plot_kwargs.update(kwargs)
        self._plt = plt
        self._figure = None
        self._axis = None
        self._axis_mosaic = None
        self._axis_dict = None
        self.make_canvas()
        # self._title = self._plot_kwargs["title"]

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

    def make_canvas(self, **kwargs):
        """Make a blank figure"""
        self._figure = self._plt.figure()

        self.reset_mosaic()
        self.set_mosaic(self._axis_mosaic)
        # self.make_canvas()
        self.update_kwargs()

    def reset_mosaic(self):
        """Reset the mosaic to a single plot"""
        self._axis_mosaic = """
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
        self._axis_mosaic = mosaic
        self._axis_dict = self._figure.subplot_mosaic(mosaic)

    def add_subplot(self, single_plot: SinglePlot, loc: str):
        """Add a subplot to the figure"""
        # ax = self._figure.add_subplot(111)
        # axis = SinglePlot()
        # self._axis_grid.append(axis)
        self._figure.add_subfigure(single_plot.fig)
        # self._axis_dict[loc] = single_plot.axis

    def update_kwargs(self):
        """Update the figure configuration"""
        self._figure.set_facecolor(self._plot_kwargs["facecolor"])
        self._figure.set_edgecolor(self._plot_kwargs["edgecolor"])

    def plot_spectrum(self):
        s = self._spectrum
        sa = s.spectral_axis
        lw = self._plot_kwargs["linewidth"]
        xunit = self._plot_kwargs["xaxis_unit"]
        yunit = self._plot_kwargs["yaxis_unit"]

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
        self._plt.show()

    def savefig(self, file, **kwargs):
        """Save the plot"""
        self.figure.savefig(
            file, facecolor=self.figure.get_facecolor(), edgecolor=self.figure.get_edgecolor()
        )  # *kwargs)


if __name__ == "__main__":
    test_grid = GridPlot()
    # test_grid.make_canvas()
    # test_grid.update_kwargs()
    test_grid.set_title("Suptitle")

    # test_mosaic = "AB;CD"
    # test_grid.set_mosaic(test_mosaic)

    test_plot_0 = SinglePlot()
    # test_plot_1 = SinglePlot()
    # test_plot_2 = SinglePlot()
    # test_plot_3 = SinglePlot()

    test_grid.add_subplot(test_plot_0, "A")
    # test_grid.add_subplot(test_plot_1, "B")
    # test_grid.add_subplot(test_plot_2, "C")
    # test_grid.add_subplot(test_plot_3, "D")

    test_grid.show()
