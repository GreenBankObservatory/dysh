"""
Plot a spectrum using matplotlib
"""

from copy import deepcopy

import astropy.units as u
import matplotlib as mpl
import numpy as np
from astropy.units.quantity import Quantity
from astropy.utils.masked import Masked
from matplotlib.widgets import SpanSelector

from dysh.log import logger

from ..coordinates import (
    decode_veldef,
    frame_to_label,
)
from ..util.docstring_manip import docstring_parameter
from .core import check_kwargs, parse_html
from .plotbase import PlotBase

_KMS = u.km / u.s


kwargs_docstring = """xaxis_unit : str or `~astropy.unit.Unit`
    The units to use on the x-axis, e.g. "km/s" to plot velocity.
yaxis_unit : str or `~astropy.unit.Unit`
    The units to use on the y-axis.
xmin : float
    Minimum x-axis value, in `xaxis_unit`.
xmax : float
    Maximum x-axis value, in `yaxis_unit`.
ymin : float
    Minimum y-axis value, in `xaxis_unit`.
ymax : float
    Maximum y-axis value, in `yaxis_unit`.
xlabel : str
    x-axis label.
ylabel : str
    y-axis label.
label : str
    Label for legend.
alpha : float
    Alpha value for the plot. Between 0 and 1.
grid : bool
    Show a plot grid or not.
figsize : tuple
    Figure size (see matplotlib).
linewidth : float
    Line width, default: 2.0.
drawstyle : str
    Line style, default 'default'.
color : str
    Line color, c also works.
title : str
    Plot title.
vel_frame : str
    The velocity frame (see VELDEF FITS Keyword).
doppler_convention: str
    The velocity convention (see VELDEF FITS Keyword).
"""


@docstring_parameter(kwargs_docstring)
class SpectrumPlot(PlotBase):
    r"""
    The SpectrumPlot class is for simple plotting of a `~dysh.spectra.spectrum.Spectrum`
    using matplotlib functions. Plots attributes are modified using keywords
    (\*\*kwargs) described below. SpectrumPlot will attempt to make smart default
    choices for the plot if no additional keywords are given.

    Parameters
    ----------
    spectrum : `~dysh.spectra.spectrum.Spectrum`
        The spectrum to plot

    Other Parameters
    ----------------
    {0}
    """

    def __init__(self, spectrum, select=True, **kwargs):
        super().__init__()
        self.reset()
        self._spectrum = spectrum
        self._sa = spectrum._spectral_axis
        self._set_xaxis_info()
        self._plot_kwargs.update(kwargs)
        self._title = self._plot_kwargs["title"]
        self._select = select
        self._selector: MultiSpanSelector = None
        self._scan_numbers = np.array([self._spectrum.meta["SCAN"]])

        self._init_selector()

        self._set_frontend()

    def _init_plot(self):
        super()._init_plot()

    def _init_selector(self):
        if self._select:
            if not hasattr(self, "_selector") or self._selector is None:
                self._selector = MultiSpanSelector(self.axes, minspan=1)
            elif hasattr(self, "_selector"):
                self._selector.clear()
        # If select is False, and there is a selector, close it.
        elif self.has_selector():
            self._selector.close()
            self._selector = None

    def _set_xaxis_info(self):
        """Ensure the xaxis info is up to date if say, the spectrum frame has changed."""
        self._plot_kwargs["doppler_convention"] = self._spectrum.doppler_convention
        self._plot_kwargs["vel_frame"] = self._spectrum.velocity_frame
        self._plot_kwargs["xaxis_unit"] = self._spectrum.spectral_axis.unit
        self._plot_kwargs["yaxis_unit"] = self._spectrum.unit

    @property
    def spectrum(self):
        """The underlying `~dysh.spectra.spectrum.Spectrum`"""
        return self._spectrum

    def default_plot_kwargs(self):
        return {
            "xmin": None,
            "xmax": None,
            "ymin": None,
            "ymax": None,
            "xlabel": None,
            "ylabel": None,
            "xaxis_unit": u.GHz,
            "yaxis_unit": None,
            "grid": False,
            "label": None,
            "alpha": 1.0,
            "figsize": None,
            "linewidth": 2.0,
            "drawstyle": "default",
            "color": None,
            "title": None,
            "doppler_convention": None,
            "vel_frame": None,
        }

    def reset(self):
        """Reset the plot keyword arguments to their defaults."""
        self._plot_kwargs = self.default_plot_kwargs()

    @docstring_parameter(kwargs_docstring)
    def plot(self, show_header=True, select=True, oshow=None, oshow_kwargs=None, **kwargs):
        """
        Plot the spectrum.

        Parameters
        ----------
        show_header : bool
            Show informational header.
        select : bool
            Allow selecting regions via click and drag.
        oshow : list or `~dysh.spectra.spectrum.Spectrum`
            Spectra to overlay in the plot.
        oshow_kwargs : dict
            Dictionary with parameters for `SpectrumPlot.oshow`.
            These include color, linestyle, label, and alpha.

        Other Parameters
        ----------------
        {0}
        """

        check_kwargs(self.default_plot_kwargs(), kwargs)

        self._select = select

        self._set_xaxis_info()
        # Plot arguments for this call of plot(). i.e. non-sticky plot attributes
        this_plot_kwargs = deepcopy(self._plot_kwargs)
        this_plot_kwargs.update(kwargs)

        if not self.has_figure():
            self._init_plot()

        if self.has_axes():
            self.axes.clear()
            self._init_selector()

        if not self._frontend.is_window_alive():
            self._set_frontend()

        if self._select:
            self._frontend.connect_buttons(self)

        # TODO: procedurally generate subplot params based on show header/buttons args.
        # ideally place left/right params right here, then top gets determined below.

        s = self._spectrum

        lw = this_plot_kwargs["linewidth"]
        if "vel_frame" not in this_plot_kwargs:
            if u.Unit(self._plot_kwargs["xaxis_unit"]).is_equivalent("km/s") and "VELDEF" in s.meta:
                # If the user specified velocity units, default to
                # the velframe the data were taken in.  This we can
                # get from VELDEF keyword.  See issue #303
                this_plot_kwargs["vel_frame"] = decode_veldef(s.meta["VELDEF"])[1].lower()
            else:
                this_plot_kwargs["vel_frame"] = s.velocity_frame
        if "chan" in str(this_plot_kwargs["xaxis_unit"]).lower():
            self._sa = u.Quantity(np.arange(len(self._sa)))
            this_plot_kwargs["xlabel"] = "Channel"
        else:
            # convert the x axis to the requested velocity frame and Doppler convention.
            self._sa = s.velocity_axis_to(
                unit=this_plot_kwargs["xaxis_unit"],
                toframe=this_plot_kwargs["vel_frame"],
                doppler_convention=this_plot_kwargs["doppler_convention"],
            )

        sf = s.flux
        if this_plot_kwargs["yaxis_unit"] is not None:
            sf = s.flux.to(this_plot_kwargs["yaxis_unit"])
        sf = Masked(sf, s.mask)

        lines = self.axes.plot(
            self._sa,
            sf,
            color=this_plot_kwargs["color"],
            lw=lw,
            drawstyle=this_plot_kwargs["drawstyle"],
            label=this_plot_kwargs["label"],
            alpha=this_plot_kwargs["alpha"],
        )
        self._line = lines[0]

        if this_plot_kwargs["label"] is not None:
            self.axes.legend()

        if this_plot_kwargs["xmin"] is not None or this_plot_kwargs["xmax"] is not None:
            self.axes.set_xlim(this_plot_kwargs["xmin"], this_plot_kwargs["xmax"])
        if this_plot_kwargs["ymin"] is not None or this_plot_kwargs["ymax"] is not None:
            self.axes.set_ylim(this_plot_kwargs["ymin"], this_plot_kwargs["ymax"])

        self.axes.tick_params(axis="both", which="both", bottom=True, top=True, left=True, right=True, direction="in")
        if this_plot_kwargs["grid"]:
            self.axes.grid(visible=True, which="major", axis="both", lw=lw / 2, color="k", alpha=0.33)
            self.axes.grid(visible=True, which="minor", axis="both", lw=lw / 2, color="k", alpha=0.22, linestyle="--")
        self._set_labels(**this_plot_kwargs)
        if self._title is not None:
            self.axes.set_title(self._title)

        if show_header:
            self.figure.subplots_adjust(top=0.79, left=0.09, right=0.95)
            self._set_header(s)
        if self._select:
            self._selector.set_minspan(self._get_dx())
            self._selector.set_xy0((self._sa.value.mean(), 0))
            self._selector.init_first_span()
            self._spectrum._selection = self._selector.get_selected_regions()
        if oshow is not None:
            if isinstance(oshow, type(self._spectrum)):
                oshow = [oshow]
            if type(oshow) is not list:
                raise TypeError(f"oshow ({oshow}) must be a list or Spectrum")
            for i, sp in enumerate(oshow):
                if not isinstance(sp, type(self._spectrum)):
                    raise TypeError(f"Element {i} of oshow ({oshow}) is not a Spectrum")
            if oshow_kwargs is None:
                oshow_kwargs = {}
            self.oshow(oshow, **oshow_kwargs)

        self.show()
        self.figure.canvas.draw_idle()

    def _get_dx(self) -> float:
        """Return the channel separation as a float."""
        return abs(self._sa[0].value - self._sa[1].value)

    def _compose_xlabel(self, **kwargs):
        """Create a sensible spectral axis label given units, velframe, and doppler convention"""
        xlabel = kwargs.get("xlabel", None)
        if xlabel is not None:
            return xlabel
        if kwargs["doppler_convention"] == "radio":
            subscript = "_{rad}"
        elif kwargs["doppler_convention"] == "optical":
            subscript = "_{opt}"
        elif kwargs["doppler_convention"] == "relativistic":
            subscript = "_{rel}"
        else:  # should never happen
            subscript = ""
        if kwargs.get("xaxis_unit", None) is not None:
            xunit = u.Unit(kwargs["xaxis_unit"])
        else:
            xunit = self.spectrum.spectral_axis.unit
        if xunit.is_equivalent(u.Hz):
            xname = r"\nu"
        elif xunit.is_equivalent(_KMS):
            xname = r"V" + subscript
        elif xunit.is_equivalent(u.angstrom):
            xname = r"\lambda"
        # Channel is handled in plot() with kwargs['xlabel']
        else:
            raise ValueError(f"Unrecognized spectral axis unit: {xunit}")
        _xunit = xunit.to_string(format="latex_inline")
        xlabel = f"{frame_to_label[kwargs['vel_frame']]} ${xname}$ ({_xunit})"
        return xlabel

    def _set_labels(self, **kwargs):
        r"""Set x and y labels according to spectral units

        Parameters
        ----------
        title : str
            Plot title.
        xlabel : str
            x-axis label.
        ylabel : str
            x-axis label.
        doppler_convention : str
            Doppler convention for x-axis.
        xaxis_unit : str
            Units for x-axis.
        yaxis_unit : str
            Units for y-axis.
        """
        title = kwargs.get("title", None)
        ylabel = kwargs.get("ylabel", None)
        if title is not None:
            self._title = title
        if kwargs.get("yaxis_unit", None) is not None:
            yunit = u.Unit(kwargs["yaxis_unit"])
        else:
            yunit = self.spectrum.unit
        self.axes.set_xlabel(self._compose_xlabel(**kwargs))
        if ylabel is not None:
            self.axes.set_ylabel(ylabel)
        else:
            if "TSCALE" in self.spectrum.meta:
                ylabel = self.spectrum.meta["TSCALE"]
            elif "TUNIT7" in self.spectrum.meta:
                tunit7 = self.spectrum.meta["TUNIT7"]
                if tunit7 == "Ta":  # what about Ta*
                    ylabel = "Ta"
                    yunit = "K"
                elif tunit7 == "Ta*":
                    ylabel = "Ta*"
                    yunit = "K"
                elif tunit7 == "Jy":
                    ylabel = "Flux"
                    yunit = "Jy"
                else:
                    ylabel = "Unknown"
                    yunit = "()"
                logger.info(f"Missing TSCALE: patching Y-axis as '{ylabel} ({yunit})'")
            self.axes.set_ylabel(f"{ylabel} ({yunit})")

    def _show_exclude(self, **kwargs):
        """TODO: Method to show the exclude array on the plot"""
        kwargs_opts = {
            "loc": "bottom",  # top,bottom ?
            "color": "silver",
        }
        kwargs_opts.update(kwargs)

    def get_selected_regions(self):
        """ """
        regions = self._selector.get_selected_regions()
        return [tuple(np.sort([np.argmin(abs(p - self._sa.value)) for p in r])) for r in regions]

    def clear_overlays(self, blines=True, oshows=True, catalog=True):
        """Clear Overlays from the plot.

        Parameters
        ----------
        blines : bool
            Remove baseline models overlaid on the plot. Default: True
        oshows : bool
            Remove other spectra overlaid on the plot. Default: True
        catalog : bool
            Remove catalog spectral lines overlaid on the plot. Default: True
        """
        if blines:
            self._clear_overlay_objects("lines", "baseline")
        if oshows:
            self._clear_overlay_objects("lines", "oshow")
        if catalog:
            self._clear_overlay_objects("lines", "catalogline")
            self._clear_overlay_objects("texts", "catalogtext")

    def clear_lines(self, gid):
        self._clear_overlay_objects("lines", gid)

    def _clear_overlay_objects(self, otype, gid):
        """
        Clears lines with `gid` from the plot.

        Parameters
        ----------
        otype : str
            Type of overlay. Can be "lines" or "texts".
        gid : str
            Group id for the lines to be cleared.
        """
        if otype == "lines":
            tgt_list = self.axes.lines
        elif otype == "texts":
            tgt_list = self.axes.texts

        for b in tgt_list:
            if b.get_gid() == gid:
                b.remove()

        self.figure.canvas.draw_idle()

    def oshow(self, spectra, color=None, linestyle=None, label=None, alpha=None):
        """
        Add `spectra` to the current plot.

        Parameters
        ----------
        spectra : list of `dysh.spectra.spectrum.Spectrum` or `dysh.spectra.spectrum.Spectrum`
            Spectra to add to the plot.
        color : list of valid `matplotlib` colors or `matplotlib` color
            Colors for the spectra. There must be one element per spectra.
        linestyle : list of valid `matplotlib` linestyles or `matplotlib` linestyle
            Linestyles for the spectra. There must be one element per spectra.
        label : list of str
            Labels for the spectra. There must be one element per spectra.
        alpha : list of float
            Alpha values for the spectra, between 0 and 1. There must be one element per spectra.
        """

        # If a single Spectrum is the input, make everything a list.
        if isinstance(spectra, type(self._spectrum)):
            spectra = [spectra]
            if color is not None:
                color = [color]
            if linestyle is not None:
                linestyle = [linestyle]
            if label is not None:
                label = [label]
            if alpha is not None:
                alpha = [alpha]

        for i, s in enumerate(spectra):
            if not isinstance(s, type(self._spectrum)):
                raise TypeError(f"Element {i} of spectra ({s}) is not a Spectrum.")

        # Pack args together, and check that we have enough of each.
        zargs = (spectra,)
        if color is not None:
            # Check that we have enough colors.
            if len(color) != len(spectra):
                raise ValueError(f"How do I color {len(spectra)} spectra with {len(color)} colors?")
            zargs += (color,)
        else:
            zargs += ([None] * len(spectra),)
        if linestyle is not None:
            # Check that we have enough linestyles.
            if len(linestyle) != len(spectra):
                raise ValueError(f"How do I style {len(spectra)} spectra with {len(linestyle)} linestyles?")
            zargs += (linestyle,)
        else:
            zargs += ([None] * len(spectra),)
        if label is not None:
            if len(label) != len(spectra):
                raise ValueError(f"How do I label {len(spectra)} spectra with {len(label)} labels?")
            zargs += (label,)
        else:
            zargs += ([None] * len(spectra),)
        if alpha is not None:
            if len(alpha) != len(spectra):
                raise ValueError(f"How do I set alpha for {len(spectra)} spectra with {len(label)} alpha values?")
            zargs += (alpha,)
        else:
            zargs += ([None] * len(spectra),)

        for s, c, ls, l, a in zip(*zargs, strict=True):
            self._oshow(s, color=c, linestyle=ls, label=l, alpha=a)

    def _oshow(self, oshow_spectrum, color=None, linestyle=None, label=None, alpha=None):
        sf = oshow_spectrum.flux.to(self._spectrum.unit)
        sa = oshow_spectrum.velocity_axis_to(
            unit=self._plot_kwargs["xaxis_unit"],
            toframe=self._plot_kwargs["vel_frame"],
            doppler_convention=self._plot_kwargs["doppler_convention"],
        )

        self.axes.plot(sa, sf, color=color, linestyle=linestyle, label=label, alpha=alpha, gid="oshow")
        if label is not None:
            self.axes.legend()
        self.update_limits()
        self.figure.canvas.draw_idle()

    def show_catalog_lines(self, rotation=0, **kwargs):
        """
        Overlay spectral lines from various catalogs on the plot, with annotations.

        Parameters
        ----------
        rotation : float, degrees
            Rotate the annotation text CCW to aid in readability. Default 0.
        **kwargs
            All other kwargs get passed to `dysh.line.query_lines`.
        """

        self.sl_tbl = self._spectrum.query_lines(**kwargs)

        fsize = 9  # font size
        num_vsteps = 7  # number of vertical steps of annotations
        rot_factor = (rotation / 90) * 0.3 / num_vsteps  # adjust ylocs to avoid rotated text running into each other
        fracstep = 0.04 + rot_factor
        ystart = 0.86 - (num_vsteps * fracstep)

        for i, line in enumerate(self.sl_tbl):
            line_name = parse_html(line["name"])
            line_freq = (
                (line["obs_frequency"] * u.MHz)
                .to(self._plot_kwargs["xaxis_unit"], equivalencies=self.spectrum.equivalencies)
                .value
            )

            vloc = ystart + (i % num_vsteps) * fracstep

            self.axes.axvline(line_freq, c="k", linewidth=1, gid="catalogline")
            self.axes.annotate(
                line_name,
                (line_freq, vloc),
                xycoords=("data", "axes fraction"),
                size=fsize,
                gid="catalogtext",
                rotation=rotation,
            )

        self.figure.canvas.draw_idle()

    def annotate_vline(self, xval, text="", rotation=0):
        """
        Add a single annotated vline to the plot. Can be cleared with the "catalog" gid.

        Parameters
        ----------
        xval : float
            X value of the line, in the same units as the plot.
        text : str
            Associated text for the vline. Defaults to an empty string.
        rotation : float
            Rotate the text CCW degrees. Default 0.
        """
        fsize = 9
        self.axes.axvline(xval, c="k", linewidth=1, gid="catalogline")
        self.axes.annotate(
            text,
            (xval, 0.7),
            xycoords=("data", "axes fraction"),
            size=fsize,
            gid="catalogtext",
            rotation=rotation,
        )
        self.figure.canvas.draw_idle()

    def update_limits(self):
        """
        Recompute the data limits based on the current artists (`~matplotlib.axes.Axes.relim`),
        and then autoscale the view limits using the data limits (`~matplotlib.axes.Axes.autoscale_view`).
        """
        self.axes.relim(visible_only=True)
        self.axes.autoscale_view()
        self.figure.canvas.draw_idle()

    def set_yaxis_unit(self, yunit: str | Quantity) -> None:
        """
        Set the units of the ordinate axis of the plot.

        Parameters
        ----------
        yunit : str or `~astropy.units.quantity.Quantity`
            New units for the y-axis. It must be a valid unit for the data.
        """

        if isinstance(yunit, str):
            yunit = u.Unit(yunit)

        y_valid_physical_types = [
            "spectral flux density",  # Flux
            "temperature",  # Kelvin
            "unknown",  # counts
        ]

        if u.get_physical_type(yunit) not in y_valid_physical_types:
            # This would be an error downstream anyway.
            pass

        ## Define the equivalency between K and Jy.
        ## Won't use for now, as we require the beam information.
        # freq = self._spectrum.spectral_axis.to("Hz").mean()
        # lmbd = ac.c/freq
        # beam = 1.06*lmbd/(100*u.m)*u.rad # 100m hardcoded...
        # beam_area = np.pi*(beam)**2/(4*np.log(2))
        # equiv = u.brightness_temperature(freq, beam_area=beam_area)
        ## Then use the equivalency with: self._spectrum.flux.to(yunit, equivalencies=equiv)
        # Here we could use one of our own exceptions, but we have none.
        try:
            yvals = self._spectrum.flux.to(yunit)
        except u.errors.UnitConversionError as e:
            raise (e)
        yvals = Masked(yvals, self._spectrum.mask)

        self.set_ydata(yvals, keep_unit=True)

        self._plot_kwargs["yaxis_unit"] = yunit
        self._set_labels(**self._plot_kwargs)

        self.update_limits()

    def set_xaxis_unit(self, xunit: str | Quantity) -> None:
        """
        Set the units of the abscissa axis of the plot.

        Parameters
        ----------
        xunit : str or `~astropy.units.quantity.Quantity`
            New units for the x-axis. It must be a valid unit for the data.
        """

        if isinstance(xunit, str) and "ch" not in xunit:
            xunit = u.Unit(xunit)

        x_valid_physical_types = [
            "speed",
            "length",
            "frequency",
            "dimensionless",  # channels
        ]

        if isinstance(xunit, Quantity):
            if u.get_physical_type(xunit) not in x_valid_physical_types:
                # This would be an error downstream anyway.
                pass

        if isinstance(xunit, str) and "ch" in xunit:
            self._sa = u.Quantity(np.arange(len(self._spectrum.data)))
            self._plot_kwargs["xlabel"] = "Channels"
        else:
            self._sa = self._spectrum.velocity_axis_to(
                unit=xunit,
                toframe=self._plot_kwargs["vel_frame"],
                doppler_convention=self._plot_kwargs["doppler_convention"],
            )
            self._plot_kwargs["xlabel"] = None

        self._line.set_xdata(self._sa.value)

        if self._selector is not None:
            self._selector.set_xy0((self._sa.value.mean(), 0))
            self._selector.set_minspan(self._get_dx())
            # For now, just clear any selected regions.
            # In the future, these should be transformed to the new
            # coordinates.
            self._selector.clear()
            # Can't add another region selector until the bounds have been
            # updated, because the ToolLineHandles uses the axes limits
            # as starting positions.

        self._plot_kwargs["xaxis_unit"] = xunit
        self._set_labels(**self._plot_kwargs)

        self.update_limits()

        # Now it should be safe to add a blank span.
        if self._selector is not None:
            self._selector.init_first_span()

    def set_ydata(self, ydata, keep_unit=False) -> None:
        """
        Set the ordinate data values for the plot line.
        That is the line that representes the underlying Spectrum.

        Parameters
        ----------
        ydata : array
            Ordinate axes values to display.
            These can be floats or quantities.
        keep_unit : bool
            Only used if `ydata` is a `~astropy.units.quantity.Quantity`.
            If True, keep the units of `ydata`.
            If False, scale to the units of the current plot line.

        """

        if isinstance(ydata, Quantity):
            if not keep_unit:
                ydata = ydata.to(self._plot_kwargs["yaxis_unit"])

        self._line.set_ydata(ydata)

        self.update_limits()

    def show_baseline(self, y, color="k") -> None:
        """
        Show a baseline model with values `y`.

        Parameters
        ----------
        y : array
            Ordinate axes values of the baseline model.
        color :
            Color for the baseline data.
        """

        self.axes.plot(self._sa, y.to(self._plot_kwargs["yaxis_unit"]), c=color, gid="baseline")

        self.figure.canvas.draw_idle()

    def _xaxis_unit_options(self):
        return ["GHz", "MHz", "kHz", "Hz", "km/s", "m/s", "channels"]

    def _yaxis_unit_options(self):
        yunit = self._spectrum.flux.unit
        yunit_type = u.get_physical_type(yunit)
        if yunit_type in ["spectral flux density"]:
            opts = ["Jy", "mJy"]
        elif yunit_type in ["temperature"]:
            opts = ["K", "mK"]
        elif yunit_type == "unknown":
            opts = ["counts"]
        else:
            raise ValueError(f"Can't handle this type of units ({yunit} with physical type: {yunit_type}).")
        return opts


class MultiSpanSelector:
    def __init__(self, ax, minspan=1, xy0=(0, 0)):
        self.ax = ax
        self.canvas = ax.figure.canvas
        self.spans = []
        self.minspan = minspan
        self.xy0 = xy0
        self.colors = {
            "edge": (0, 0, 0, 1),
            "face": (0, 0, 0, 0.3),
            "edge_selected": (*mpl.colors.to_rgb("#6c3483"), 1.0),
        }
        # Register callbacks before creating any spans.
        self.cid_press = self.canvas.mpl_connect("button_press_event", self.on_press)

        self.active_span = None
        self.selected_span = None

    def init_first_span(self):
        """
        Add the first span selector to this object.
        """
        if len(self.spans) > 0:
            return
        else:
            self.spans.extend(self.init_span())
            self.active_span = self.spans[0]

    def init_span(self):
        span = SpanSelector(
            self.ax,
            self.on_select,
            direction="horizontal",
            useblit=False,  # blitting causes spans to blink on press.
            interactive=True,
            drag_from_anywhere=True,
            ignore_event_outside=True,
            props=dict(facecolor=self.colors["face"], alpha=0.3, xy=self.xy0),
            minspan=self.minspan,
        )
        return [span]

    def on_select(self, vmin, vmax):
        span = vmax - vmin
        if span > self.minspan and np.all(np.diff(self.get_selected_regions(False)) > self.minspan):
            if self.active_span is not None:
                self.active_span.set_active(False)
                self.active_span = None
            self.spans.extend(self.init_span())
            self.active_span = self.spans[-1]
        elif self.active_span is not None:
            self.active_span.set_active(False)
            self.active_span = None
        return

    def on_press(self, event):
        # Do nothing if outside the axes.
        if event.inaxes != self.ax:
            return

        # Do nothing if another widget is enabled.
        if self.ax.get_navigate_mode() is not None:
            return

        # Determine if the event is in a span.
        if len(self.spans) > 1:
            got_one = False
            for span in self.spans:
                # Select only a single span at a time.
                if span._contains(event) and not got_one and np.diff(span.extents) > self.minspan:
                    got_one = True
                    self.active_span = span
                    self.active_span.set_active(True)
                    self.active_span.set_props(**{"linewidth": 20, "edgecolor": self.colors["edge_selected"]})
                    self.selected_span = span
                else:
                    span.set_active(False)
                    props = {"linewidth": 1, "edgecolor": self.colors["edge"]}
                    span.set_props(**props)

            if not got_one:
                self.active_span = None
                self.selected_span = None
                for span in self.spans:
                    # Determine if there's a span that needs to be completed
                    # and activate it.
                    if np.diff(span.extents) <= self.minspan:
                        self.active_span = span
                        self.active_span.set_active(True)

    def disconnect(self):
        """Disconnect all event handlers to prevent memory leaks and dangling references."""
        if hasattr(self, "cid_press") and self.cid_press is not None:
            self.canvas.mpl_disconnect(self.cid_press)
            self.cid_press = None

    def clear(self):
        """Clear all regions."""
        for span in self.spans:
            span.clear()
            span.disconnect_events()
            del span
        self.spans.clear()
        self.spans = []
        self.selected_span = None
        self.active_span = None

    def clear_region(self, event=None):
        """Clear a single region."""
        if not self.selected_span:
            return
        self.selected_span.clear()
        self.selected_span.disconnect_events()
        self.spans.remove(self.selected_span)
        del self.selected_span
        self.selected_span = None
        self.active_span = None
        # If there's only one more span
        # activate it.
        if len(self.spans) == 1:
            self.spans[0].set_active(True)
            self.active_span = self.spans[0]

    def clear_regions(self, event=None):
        """Clear all regions and start a new empty region."""
        self.clear()
        # Add a new span.
        if self.ax is not None:
            self.spans.extend(self.init_span())
            self.active_span = self.spans[-1]

    def close(self):
        self.disconnect()
        self.clear()
        self.ax = None
        self.canvas = None

    def get_selected_regions(self, ignore_incomplete=True):
        """
        Parameters
        ----------
        ignore_complete : bool
            If True ignore spans that are smaller than `self.minspan`.

        Returns
        -------
        regions : list of tuples
            List with edges of the spans as tuples.
        """
        if ignore_incomplete:
            regions = [span.extents for span in self.spans if np.diff(span.extents) > self.minspan]
        else:
            regions = [span.extents for span in self.spans]
        return regions

    def set_minspan(self, minspan: float) -> None:
        """Set the minimum span in data units for the spans."""
        self.minspan = minspan
        for span in self.spans:
            span.minspan = minspan

    def set_xy0(self, xy0: tuple) -> None:
        """Set the bottom left coordinates in data units for the spans and their handles."""
        self.xy0 = xy0
        for span in self.spans:
            if np.diff(span.extents) < self.minspan:
                span._selection_artist.update({"xy": self.xy0})
                span._edge_handles.set_data([self.xy0[0], self.xy0[0]])
