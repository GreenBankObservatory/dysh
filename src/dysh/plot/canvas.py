import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg

from dysh.coordinates import frame_to_label

_KMS = u.km / u.s


class InteractiveFigure1D(FigureCanvasQTAgg):
    """Figure object for interactive plotting"""

    def __init__(self, parent=None):
        self._plt = plt
        self._figure, self._axis = plt.subplots()
        # super().__init__(self._figure)
        if parent is not None:
            self.setParent(parent)
        self._plt.ion()

    @property
    def plt(self):
        return self._plt

    @property
    def figure(self):
        return self._figure

    @property
    def axis(self):
        return self._axis


class StaticFigure1D(FigureCanvasAgg):
    """Figure object for static plotting"""

    def __init__(self):
        self._plt = plt
        self._figure, self._axis = plt.subplots()
        self._plt.ioff()

    @property
    def plt(self):
        return self._plt

    @property
    def figure(self):
        return self._figure

    @property
    def axis(self):
        return self._axis


class SpectrumPlotCanvas:
    """Plot canvas for SpectrumPlot"""

    def __init__(self, interactive=False):
        self._interactive = interactive
        self._set_backend()

    def _set_backend(self):
        """Set the backend based on if the plot is interactive or not"""
        if self.interactive:
            self._backend = "InteractiveFigure1D"
            InteractiveFigure1D.__init__(self)
        else:
            self._backend = "StaticFigure1D"
            StaticFigure1D.__init__(self)

    @property
    def interactive(self):
        """Boolean representing if the plot is interactive or not"""
        return self._interactive

    @property
    def backend(self):
        """The backend for the figure canvas"""
        return self._backend

    @backend.setter
    def set_backend(self, value):
        """Set the backend value"""
        if isinstance(value, str):
            self._backend = value
        else:
            raise ValueError(f"Backend must be a string. Instead received {value} of type {type(value)}.")


class SpectrumPlot:
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

    def __init__(self, spectrum, interactive=False, **kwargs):
        self.reset()
        self._interactive = interactive
        self._canvas = None
        self._set_canvas()
        self._spectrum = spectrum
        self._set_xaxis_info()
        self._set_yaxis_info()
        self._plot_kwargs.update(kwargs)
        self._title = self._plot_kwargs["title"]
        self.plot()

    def _set_canvas(self):
        self._canvas = SpectrumPlotCanvas(self.interactive)
        self._plt = self._canvas._plt
        self._figure = self._canvas._figure
        self._axis = self._canvas._axis

    def _set_xaxis_info(self):
        """Ensure the xaxis info is up to date if say, the spectrum frame has changed."""
        self._plot_kwargs["doppler_convention"] = self._spectrum.doppler_convention
        self._plot_kwargs["vel_frame"] = self._spectrum.velocity_frame
        self._plot_kwargs["xaxis_unit"] = self._spectrum.spectral_axis.unit

    def _set_yaxis_info(self):
        self._plot_kwargs["yaxis_unit"] = self._spectrum.unit

    @property
    def canvas(self):
        """The plotting canvas"""
        return self._canvas

    @property
    def interactive(self):
        """Boolean representing if the plot is interactive or not"""
        return self._interactive

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
        r"""
        Plot the spectrum.

        Parameters
        ----------
        **kwargs : various
            keyword=value arguments (need to describe these in a central place)
        """

        self._set_xaxis_info()

        s = self._spectrum
        sa = s.spectral_axis
        lw = self._plot_kwargs["linewidth"]
        xunit = self._plot_kwargs["xaxis_unit"]
        yunit = self._plot_kwargs["yaxis_unit"]
        if "vel_frame" not in self._plot_kwargs:
            self._plot_kwargs["vel_frame"] = s.velocity_frame
        if xunit is None:
            xunit = str(sa.unit)
        if "chan" in str(xunit).lower():
            sa = np.arange(len(sa))
            self._plot_kwargs["xlabel"] = "Channel"
        else:
            # convert the x axis to the requested
            # print(f"EQUIV {equiv} doppler_rest {sa.doppler_rest} [{rfq}] convention {convention}")
            # sa = s.spectral_axis.to( self._plot_kwargs["xaxis_unit"], equivalencies=equiv,doppler_rest=rfq, doppler_convention=convention)
            sa = s.velocity_axis_to(
                unit=xunit,
                toframe=self._plot_kwargs["vel_frame"],
                doppler_convention=self._plot_kwargs["doppler_convention"],
            )
        sf = s.flux
        if yunit is not None:
            sf = s.flux.to(yunit)
        self._axis.plot(sa, sf, color=self._plot_kwargs["color"], lw=lw)
        self._axis.set_xlim(self._plot_kwargs["xmin"], self._plot_kwargs["xmax"])
        self._axis.set_ylim(self._plot_kwargs["ymin"], self._plot_kwargs["ymax"])
        self._axis.tick_params(axis="both", which="both", bottom=True, top=True, left=True, right=True, direction="in")
        if self._plot_kwargs["grid"]:
            self._axis.grid(visible=True, which="major", axis="both", lw=lw / 2, color="k", alpha=0.33)
            self._axis.grid(visible=True, which="minor", axis="both", lw=lw / 2, color="k", alpha=0.22, linestyle="--")

        self._set_labels(**self._plot_kwargs)
        # self._axis.axhline(y=0,color='red',lw=2)
        if self._title is not None:
            self._axis.set_title(self._title)

    def reset(self):
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
        }

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
            self._plt.show()

    def show(self):
        self.refresh()

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
