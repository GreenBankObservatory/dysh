"""
Plot a spectrum using matplotlib
"""

from copy import deepcopy

import astropy.units as u
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

matplotlib.use("QtAgg")
from ..coordinates import frame_to_label

_KMS = u.km / u.s


class SpectrumPlot:
    # @todo make xaxis_unit='chan[nel]' work
    r"""
    The SpectrumPlot class is for simple plotting of a `~spectrum.Spectrum`
    using matplotlib functions. Plots attributes are modified using keywords
    (\*\*kwargs) described below SpectrumPlot will attempt to make smart default
    choices for the plot if no additional keywords are given.

    Parameters
    ----------
    spectrum : `~spectra.spectrum.Spectrum`
        The spectrum to plot
    **kwargs : dict
        Plot attribute keyword arguments, see below.

    Other Parameters
    ----------------
    xaxis_unit : str or `~astropy.unit.Unit`
        The units to use on the x-axis, e.g. "km/s" to plot velocity
    yaxis_unit : str or `~astropy.unit.Unit`
        The units to use on the y-axis
    xmin : float
        Minimum x-axis value, in `xaxis_unit`
    xmax : float
        Maximum x-axis value, in `yaxis_unit`
    ymin : float
        Minimum y-axis value, in `xaxis_unit`
    ymax : float
        Maximum y-axis value, in `yaxis_unit`
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
    vel_frame : str
        The velocity frame (see VELDEF FITS Keyword)
    doppler_convention: str
        The velocity convention (see VELDEF FITS Keyword)
    """

    # loc, legend, bbox_to_anchor

    def __init__(self, spectrum, **kwargs):
        self.reset_kwargs()
        self._spectrum = spectrum
        self._set_xaxis_info()
        self._plt = plt
        self._figure = None
        self._axis = None
        self._title = self._plot_kwargs["title"]
        self._show_masked = True

        # self.update_kwargs()
        self._plot_kwargs.update(kwargs)

    # def __call__ (see pyspeckit)

    def _set_xaxis_info(self):
        """Ensure the xaxis info is up to date if say, the spectrum frame has changed."""
        self._plot_kwargs["doppler_convention"] = self._spectrum.doppler_convention
        self._plot_kwargs["vel_frame"] = self._spectrum.velocity_frame
        self._plot_kwargs["xaxis_unit"] = self._spectrum.spectral_axis.unit
        self._plot_kwargs["yaxis_unit"] = self._spectrum.unit

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

    def plot(self, **kwargs):
        # @todo document kwargs here
        r"""
        Plot the spectrum.

        Parameters
        ----------
        **kwargs : various
            keyword=value arguments (need to describe these in a central place)
        """
        # xtype = 'velocity, 'frequency', 'wavelength'
        # if self._figure is None:

        self._set_xaxis_info()
        # plot arguments for this call of plot(). i.e. non-sticky plot attributes
        this_plot_kwargs = deepcopy(self._plot_kwargs)
        this_plot_kwargs.update(kwargs)
        if True:  # @todo deal with plot reuse (notebook vs script)
            self._figure, self._axis = self._plt.subplots(figsize=this_plot_kwargs["figsize"])
        # else:
        #    self._axis.cla()

        s = self._spectrum
        sa = s.spectral_axis
        lw = this_plot_kwargs["linewidth"]
        xunit = this_plot_kwargs["xaxis_unit"]
        yunit = this_plot_kwargs["yaxis_unit"]
        if "vel_frame" not in this_plot_kwargs:
            this_plot_kwargs["vel_frame"] = s.velocity_frame
        if xunit is None:
            xunit = str(sa.unit)
        if "chan" in str(xunit).lower():
            sa = np.arange(len(sa))
            this_plot_kwargs["xlabel"] = "Channel"
        else:
            # convert the x axis to the requested
            # print(f"EQUIV {equiv} doppler_rest {sa.doppler_rest} [{rfq}] convention {convention}")
            # sa = s.spectral_axis.to( self._plot_kwargs["xaxis_unit"], equivalencies=equiv,doppler_rest=rfq, doppler_convention=convention)
            sa = s.velocity_axis_to(
                unit=xunit,
                toframe=this_plot_kwargs["vel_frame"],
                doppler_convention=this_plot_kwargs["doppler_convention"],
            )
        sf = s.flux
        if yunit is not None:
            sf = s.flux.to(yunit)
        plot_indx = self._spectrum.mask == False
        self._axis.plot(sa[plot_indx], sf[plot_indx], color=this_plot_kwargs["line_color"], lw=lw)
        if self._show_masked:
            self._axis.plot(sa[~plot_indx], sf[~plot_indx], color=this_plot_kwargs["line_color_masked"], lw=lw)
        self._axis.set_xlim(this_plot_kwargs["xmin"], this_plot_kwargs["xmax"])
        self._axis.set_ylim(this_plot_kwargs["ymin"], this_plot_kwargs["ymax"])
        self._axis.tick_params(axis="both", which="both", bottom=True, top=True, left=True, right=True, direction="in")
        if this_plot_kwargs["grid"]:
            self._axis.grid(visible=True, which="major", axis="both", lw=lw / 2, color="k", alpha=0.33)
            self._axis.grid(visible=True, which="minor", axis="both", lw=lw / 2, color="k", alpha=0.22, linestyle="--")

        self._set_labels(**this_plot_kwargs)
        # self._axis.axhline(y=0,color='red',lw=2)
        if self._title is not None:
            self._axis.set_title(self._title, color=this_plot_kwargs["titlecolor"])

        self._figure.set_facecolor(this_plot_kwargs["facecolor"])
        self._figure.set_edgecolor(this_plot_kwargs["edgecolor"])
        self._axis.set_facecolor(this_plot_kwargs["facecolor"])

        self._axis.title.set_color(this_plot_kwargs["title_color"])
        self._axis.xaxis.label.set_color(this_plot_kwargs["xlabel_color"])
        self._axis.yaxis.label.set_color(this_plot_kwargs["ylabel_color"])

        self._axis.tick_params(axis="x", colors=this_plot_kwargs["tick_color"])
        self._axis.tick_params(axis="y", colors=this_plot_kwargs["tick_color"])

        self._axis.spines["left"].set_color(this_plot_kwargs["spine_color"])
        self._axis.spines["right"].set_color(this_plot_kwargs["spine_color"])
        self._axis.spines["top"].set_color(this_plot_kwargs["spine_color"])
        self._axis.spines["bottom"].set_color(this_plot_kwargs["spine_color"])

        self.refresh()

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
            "xlabel_color": "white",
            "ylabel_color": "white",
            "tick_color": "white",
            "spine_color": "white",
            #'capsize':3,
            "linewidth": 2.0,
            "linestyle": "steps-mid",
            "markersize": 8,
            "line_color": "white",
            "line_color_masked": "red",
            "title": None,
            "title_color": "white",
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

    def update_kwargs(self):
        """Update the figure configuration"""
        self._figure.set_facecolor(self._plot_kwargs["facecolor"])
        self._figure.set_edgecolor(self._plot_kwargs["edgecolor"])

    def set_title(self, title: str):
        """Add or update the title of the figure"""
        self._figure.suptitle(title, color=self._plot_kwargs["titlecolor"])

    def plot_selection(self):
        if self._figure is not None:
            pass

    def _compose_xlabel(self, **kwargs):
        """Create a sensible spectral axis label given units, velframe, and doppler convention"""
        xlabel = kwargs.get("xlabel", None)
        if xlabel is not None:
            return xlabel
        if kwargs["doppler_convention"] == "radio":
            subscript = "_{rad}$"
        elif kwargs["doppler_convention"] == "optical":
            subscript = "_{opt}$"
        elif kwargs["doppler_convention"] == "relativistic":
            subscript = "_{rel}$"
        else:  # should never happen
            subscript = ""
        if kwargs.get("xaxis_unit", None) is not None:
            xunit = u.Unit(kwargs["xaxis_unit"])
        else:
            xunit = self.spectrum.spectral_axis.unit
        if xunit.is_equivalent(u.Hz):
            xname = r"$\nu" + subscript
        elif xunit.is_equivalent(_KMS):
            xname = r"V$" + subscript
        elif xunit.is_equivalent(u.angstrom):
            xname = r"$\lambda" + subscript
        # Channel is handled in plot() with kwargs['xlabel']
        else:
            raise ValueError(f"Unrecognized spectral axis unit: {xunit}")
        xlabel = f"{frame_to_label[kwargs['vel_frame']]} {xname} ({xunit})"
        return xlabel

    def _set_labels(self, **kwargs):
        r"""Set x and y labels according to spectral units

        Parameters
        ----------

        **kwargs : various
            title : str
                plot title
            xlabel : str
                x-axis label
            ylabel : str
                x-axis label

            and other keyword=value arguments
        """
        title = kwargs.get("title", None)
        xlabel = kwargs.get("xlabel", None)
        ylabel = kwargs.get("ylabel", None)
        if title is not None:
            self._title = title
        if kwargs.get("yaxis_unit", None) is not None:
            yunit = u.Unit(kwargs["yaxis_unit"])
        else:
            yunit = self.spectrum.unit
        self.axis.set_xlabel(self._compose_xlabel(**kwargs))
        if ylabel is not None:
            self.axis.set_ylabel(ylabel)
        elif yunit.is_equivalent(u.K):
            self.axis.set_ylabel(f"$T_A$ ({yunit})")
        elif self.spectrum.unit.is_equivalent(u.Jy):
            snu = r"$S_{\nu}$"
            self.axis.set_ylabel(f"{snu} ({yunit})")

    def _show_exclude(self, **kwargs):
        """TODO: Method to show the exclude array on the plot"""
        kwargs_opts = {
            "loc": "bottom",  # top,bottom ?
            "color": "silver",
        }
        kwargs_opts.update(kwargs)
        # if kwargs_opts['loc'] == 'bottom':
        #    self._ax.axhline

    def refresh(self):
        """Refresh the plot"""
        if self.axis is not None:
            self.axis.figure.canvas.draw()
            # self.axis.figure.canvas.draw_idle()
            self._plt.show()

    def show(self):
        """Show the plot"""
        self._figure.show()

    def savefig(self, file, **kwargs):
        r"""Save the plot

        Parameters
        ----------
        file - str
            The output file name
        **kwargs : dict or key=value pairs
            Other arguments to pass to `~matplotlib.pyplot.savefig`

        """
        self.figure.savefig(file, *kwargs)
