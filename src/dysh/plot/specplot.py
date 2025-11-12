"""
Plot a spectrum using matplotlib
"""

from copy import deepcopy

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import EarthLocation
from astropy.utils.masked import Masked
from matplotlib.patches import Rectangle
from matplotlib.widgets import Button, CheckButtons, SpanSelector

from ..coordinates import (
    decode_veldef,
    frame_to_label,
)
from ..util.docstring_manip import docstring_parameter
from . import PlotBase, check_kwargs

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

    def __init__(self, spectrum, **kwargs):
        super().__init__()
        self._spectrum = spectrum
        self._sa = spectrum._spectral_axis
        self._set_xaxis_info()
        self._plot_kwargs.update(kwargs)
        self._title = self._plot_kwargs["title"]
        self._selector: InteractiveSpanSelector = None
        self._freezey = (self._plot_kwargs["ymin"] is not None) or (self._plot_kwargs["ymax"] is not None)
        self._freezex = (self._plot_kwargs["xmin"] is not None) or (self._plot_kwargs["xmax"] is not None)
        self._scan_numbers = np.array([self._spectrum.meta["SCAN"]])

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
            "xaxis_unit": None,
            "yaxis_unit": None,
            "grid": False,
            "figsize": None,
            "linewidth": 2.0,
            "drawstyle": "default",
            "color": None,
            "title": None,
        }

    def reset(self):
        """Reset the plot keyword arguments to their defaults."""
        self._plot_kwargs = self.default_plot_kwargs()

    @docstring_parameter(kwargs_docstring)
    def plot(self, show_header=True, select=True, oshow=None, hist=False, **kwargs):
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
        hist : bool
            Set the plot type to histogram. Default: False

        Other Parameters
        ----------------
        {0}
        """

        check_kwargs(self.default_plot_kwargs(), kwargs)

        self._set_xaxis_info()
        # Plot arguments for this call of plot(). i.e. non-sticky plot attributes
        this_plot_kwargs = deepcopy(self._plot_kwargs)
        this_plot_kwargs.update(kwargs)
        if True:  # @todo deal with plot reuse (notebook vs script)
            self._figure, self._axis = self._plt.subplots(figsize=(10, 6))

        # TODO: procedurally generate subplot params based on show header/buttons args.
        # ideally place left/right params right here, then top gets determined below.

        s = self._spectrum

        lw = this_plot_kwargs["linewidth"]
        self._xunit = this_plot_kwargs["xaxis_unit"]  # need to kick back a ref to xunit for baseline overlays
        self._yunit = this_plot_kwargs["yaxis_unit"]
        if self._xunit is None:
            self._xunit = str(sa.unit)  # noqa: F821
        if "vel_frame" not in this_plot_kwargs:
            if u.Unit(self._xunit).is_equivalent("km/s") and "VELDEF" in s.meta:
                # If the user specified velocity units, default to
                # the velframe the data were taken in.  This we can
                # get from VELDEF keyword.  See issue #303
                this_plot_kwargs["vel_frame"] = decode_veldef(s.meta["VELDEF"])[1].lower()
            else:
                this_plot_kwargs["vel_frame"] = s.velocity_frame
        if "chan" in str(self._xunit).lower():
            self._sa = u.Quantity(np.arange(len(self._sa)))
            this_plot_kwargs["xlabel"] = "Channel"
        else:
            # convert the x axis to the requested velocity frame and Doppler convention.
            self._sa = s.velocity_axis_to(
                unit=self._xunit,
                toframe=this_plot_kwargs["vel_frame"],
                doppler_convention=this_plot_kwargs["doppler_convention"],
            )

        sf = s.flux
        if self._yunit is not None:
            sf = s.flux.to(self._yunit)
        sf = Masked(sf, s.mask)

        lines = self._axis.plot(
            self._sa, sf, color=this_plot_kwargs["color"], lw=lw, drawstyle=this_plot_kwargs["drawstyle"]
        )
        self._line = lines[0]

        if not this_plot_kwargs["xmin"] and not this_plot_kwargs["xmax"]:
            self._axis.set_xlim(np.min(self._sa).value, np.max(self._sa).value)
        else:
            self._axis.set_xlim(this_plot_kwargs["xmin"], this_plot_kwargs["xmax"])

        if self._freezey:
            self._axis.autoscale(enable=False)
        else:
            self._axis.autoscale(axis="y", enable=True)
        self._axis.set_ylim(this_plot_kwargs["ymin"], this_plot_kwargs["ymax"])

        self._axis.tick_params(axis="both", which="both", bottom=True, top=True, left=True, right=True, direction="in")
        if this_plot_kwargs["grid"]:
            self._axis.grid(visible=True, which="major", axis="both", lw=lw / 2, color="k", alpha=0.33)
            self._axis.grid(visible=True, which="minor", axis="both", lw=lw / 2, color="k", alpha=0.22, linestyle="--")
        self._set_labels(**this_plot_kwargs)
        if self._title is not None:
            self._axis.set_title(self._title)

        # persistent zline reference
        self._zline = self._axis.axhline(0, c="k", visible=False, linewidth=1)

        if show_header:
            self._figure.subplots_adjust(top=0.7, left=0.09, right=0.95)
            self._set_header(s)
            self._menu = Menu(self)
        if select:
            self._selector = InteractiveSpanSelector(self._axis)
            self._spectrum._selection = self._selector.get_selected_regions()
        if oshow is not None:
            if isinstance(oshow, type(self._spectrum)):
                oshow = [oshow]
            if type(oshow) is not list:
                raise TypeError(f"oshow ({oshow}) must be a list or Spectrum")
            for i, sp in enumerate(oshow):
                if not isinstance(sp, type(self._spectrum)):
                    raise TypeError(f"Element {i} of oshow ({oshow}) is not a Spectrum")
                self._oshow(sp)

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
        self.axis.set_xlabel(self._compose_xlabel(**kwargs))
        if ylabel is not None:
            self.axis.set_ylabel(ylabel)
        else:
            # @todo It would be nice if yunit could be latex. e.g. T_A^* instead of Ta*
            self.axis.set_ylabel(f"{self.spectrum.meta['TSCALE']} ({yunit})")

    def _show_exclude(self, **kwargs):
        """TODO: Method to show the exclude array on the plot"""
        kwargs_opts = {
            "loc": "bottom",  # top,bottom ?
            "color": "silver",
        }
        kwargs_opts.update(kwargs)
        # if kwargs_opts['loc'] == 'bottom':
        #    self._ax.axhline

    def get_selected_regions(self):
        """ """
        regions = self._selector.get_selected_regions()
        return [tuple(np.sort([np.argmin(abs(p - self._sa.value)) for p in r])) for r in regions]

    def freex(self):
        """Free the X-axis if limits have been set. Resets the limits to be the span of the spectrum."""
        self._freezex = False
        mins = []
        maxs = []
        for line in self._axis.lines:
            mins.append(line._x.min())
            maxs.append(line._x.max())
        self._axis.set_xlim((min(mins), max(maxs)))

    def freey(self):
        """Free the Y-axis if limits have been set. Autoscales the Y-axis according to your matplotlib configuration."""
        self._freezey = False
        self._axis.relim()
        self._axis.autoscale(axis="y", enable=True)
        self._axis.autoscale_view()

    def freexy(self):
        r"""Free the X and Y axes simultaneously. See `freex` and `freey` for more details."""
        self.freex()
        self.freey()

    def clear_overlays(self, blines=True, oshows=True):
        """Clear Overlays from the plot.

        Parameters
        ----------
        blines : bool
            Remove only baseline models overlaid on the plot. Default: True
        """
        if blines:
            self.clear_lines("baseline")
        if oshows:
            self.clear_lines("oshow")

    def clear_lines(self, gid):
        """
        Clears lines with `gid` from the plot.

        Parameters
        ----------
        gid : str
            Group id for the lines to be cleared.
        """

        for b in self._axis.lines:
            if b.get_gid() == gid:
                b.remove()

    def oshow(self, spectra, color=None, linestyle=None):
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
        """

        # If a single Spectrum is the input, make everything a list.
        if isinstance(spectra, type(self._spectrum)):
            spectra = [spectra]
            if color is not None:
                color = [color]
            if linestyle is not None:
                linestyle = [linestyle]

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

        for s, c, ls in zip(*zargs, strict=True):
            self._oshow(s, color=c, linestyle=ls)

    def _oshow(self, oshow_spectrum, color=None, linestyle=None):
        this_plot_kwargs = deepcopy(self._plot_kwargs)
        sf = oshow_spectrum.flux.to(self._spectrum.unit)
        sa = oshow_spectrum.velocity_axis_to(
            unit=self._xunit,
            toframe=this_plot_kwargs["vel_frame"],
            doppler_convention=this_plot_kwargs["doppler_convention"],
        )

        self._axis.plot(sa, sf, color=color, linestyle=linestyle, gid="oshow")

        self.freexy()


class InteractiveSpanSelector:
    def __init__(self, ax):
        self.ax = ax
        self.canvas = ax.figure.canvas
        self.regions = []
        self.active_patch = None
        self.press = None
        self.dragging_edge = None
        self.edge_threshold = 0.03  # in axes fraction
        self.colors = {
            "edge": (0, 0, 0, 1),
            "face": (0, 0, 0, 0.3),
            "edge_selected": plt.matplotlib.colors.to_rgb("#6c3483") + (1.0,),  # noqa: RUF005
        }

        # SpanSelector for creating new regions.
        self.span = SpanSelector(
            ax,
            self.onselect,
            direction="horizontal",
            useblit=True,
            interactive=False,
            props=dict(facecolor=self.colors["face"], alpha=0.3),
        )

        # Button to clear all selections.
        self.region_clear_button_ax = self.canvas.figure.add_axes([0.1, 0.025, 0.12, 0.04])
        self.region_clear_button = Button(self.region_clear_button_ax, "Clear Regions")
        self.region_clear_button.on_clicked(self.clear_regions)

        # Button to clear a single region.
        self.region_del_button_ax = self.canvas.figure.add_axes([0.24, 0.025, 0.12, 0.04])
        self.region_del_button = Button(self.region_del_button_ax, "Delete Region")
        self.region_del_button.on_clicked(self.clear_region)

        # Connect interaction events for dragging/resizing
        self.cid_press = self.canvas.mpl_connect("button_press_event", self.on_press)
        self.cid_release = self.canvas.mpl_connect("button_release_event", self.on_release)
        self.cid_motion = self.canvas.mpl_connect("motion_notify_event", self.on_motion)
        self.cid_key = plt.gcf().canvas.mpl_connect("key_press_event", self.on_key_press)

    def onselect(self, vmin, vmax):
        if abs(vmax - vmin) < 1e-6:
            return  # ignore tiny selections
        rect = Rectangle(
            (vmin, 0),
            vmax - vmin,
            1,
            transform=self.ax.get_xaxis_transform(),
            facecolor=self.colors["face"],
            edgecolor=self.colors["edge"],
        )
        self.ax.add_patch(rect)
        self.regions.append(rect)
        self.canvas.draw()

    def on_press(self, event):
        if event.inaxes != self.ax:
            return
        # Do nothing if another widget is enabled.
        if self.ax.get_navigate_mode() is not None:
            return
        got_one = False
        for patch in self.regions:
            contains, _attr = patch.contains(event)
            if contains and not got_one:
                self.span.set_active(False)
                x0 = patch.get_x()
                x1 = x0 + patch.get_width()
                xtol = self.edge_threshold * (self.ax.get_xlim()[1] - self.ax.get_xlim()[0])
                if abs(event.xdata - x0) <= xtol or abs(event.xdata + x0) <= xtol:
                    self.dragging_edge = "left"
                elif abs(event.xdata - x1) <= xtol or abs(event.xdata + x1) <= xtol:
                    self.dragging_edge = "right"
                else:
                    self.dragging_edge = "move"
                self.active_patch = patch
                self.press = event.xdata, x0, x1
                self.active_patch.set_edgecolor(self.colors["edge_selected"])
                self.active_patch.set_linewidth(2.0)
                got_one = True
            else:
                patch.set_edgecolor(self.colors["edge"])
                patch.set_linewidth(1.0)

    def on_motion(self, event):
        if not self.active_patch or event.inaxes != self.ax or self.press is None:
            return
        xdata, x0, x1 = self.press
        dx = event.xdata - xdata
        if self.dragging_edge == "move":
            new_x = x0 + dx
            self.active_patch.set_x(new_x)
        elif self.dragging_edge == "left":
            new_x0 = x0 + dx
            if new_x0 < x1:
                self.active_patch.set_x(new_x0)
                self.active_patch.set_width(x1 - new_x0)
        elif self.dragging_edge == "right":
            new_x1 = x1 + dx
            if new_x1 > x0:
                self.active_patch.set_width(new_x1 - x0)
        self.canvas.draw_idle()

    def on_release(self, event):
        self.press = None
        self.dragging_edge = None
        self.span.set_active(True)

    def on_key_press(self, event):
        if event.key == "d":
            self.clear_region(event)
        if event.key == "D":
            self.clear_regions(event)

    def clear_regions(self, event=None):
        for patch in self.regions:
            patch.remove()
        self.regions.clear()
        self.canvas.draw_idle()

    def clear_region(self, event=None):
        if not self.active_patch:
            return
        idx = self.regions.index(self.active_patch)  # noqa: F841
        self.regions.remove(self.active_patch)
        self.active_patch.remove()
        self.active_patch = None
        self.canvas.draw_idle()

    def get_selected_regions(self):
        return [(patch.get_x(), patch.get_x() + patch.get_width()) for patch in self.regions]


class Menu:
    def __init__(self, specplot):
        self.specplot = specplot
        self.canvas = self.specplot._axis.figure.canvas
        self.regionshow = True
        self.gbtloc = EarthLocation.of_site("Green Bank Telescope")

        hcoords = [0.08, 0.21, 0.34, 0.47, 0.60, 0.73]
        vcoords = [0.93, 0.88]

        hsize = 0.12
        vsize = 0.04

        # Button to write an ASCII file.
        # self.writeascii_button_ax = self.canvas.figure.add_axes([hcoords[0], vcoords[0], hsize, vsize])
        # self.writeascii_button = Button(self.writeascii_button_ax, "Write ASCII")
        # self.writeascii_button.on_clicked(self._writeascii)

        # CheckButton to show zline
        self.zline_button_ax = self.canvas.figure.add_axes([hcoords[1], vcoords[0], hsize, vsize])
        self.zline_button = CheckButtons(
            ax=self.zline_button_ax,
            labels=["Zline"],
        )
        self.zline_button.on_clicked(self._zline_enable)

        # CheckButton to plot hist style
        self.hist_button_ax = self.canvas.figure.add_axes([hcoords[1], vcoords[1], hsize, vsize])
        self.hist_button = CheckButtons(
            ax=self.hist_button_ax,
            labels=["Hist"],
        )
        self.hist_button.on_clicked(self._hist_enable)

        # CheckButton to show selected regions
        self.regionshow_button_ax = self.canvas.figure.add_axes([hcoords[2], vcoords[0], hsize, vsize])
        self.regionshow_button = CheckButtons(
            ax=self.regionshow_button_ax,
            labels=["Show Regs"],
        )
        self.regionshow_button.on_clicked(self._regionshow_enable)

        # Button to clear overlays
        self.clearoverlay_button_ax = self.canvas.figure.add_axes([hcoords[2], vcoords[1], hsize, vsize])
        self.clearoverlay_button = Button(self.clearoverlay_button_ax, "Clear Overlays")
        self.clearoverlay_button.on_clicked(self._clearoverlays)

        # Button/Radio combo to select leftclick
        # self.leftclick_button_ax = self.canvas.figure.add_axes([hcoords[3], vcoords[0], hsize, vsize])
        # self.leftclick_button = Button(self.leftclick_button_ax, "Left Click")
        # self.leftclick_button.on_clicked(self.open_leftclick_radio)
        # self.leftclick_radio_ax = self.canvas.figure.add_axes([hcoords[3], vcoords[1]-vsize*4, hsize*2, vsize*4],
        # zorder=100)
        # self.leftclick_radio_ax.set_visible(False)

        # Button to cycle xunit
        self.xunit_cycle = {
            "chan": "chan",
            "Hz": u.Hz,
            "kHz": u.kHz,
            "MHz": u.MHz,
            "GHz": u.GHz,
            "m/s": u.m / u.s,
            "km/s": u.km / u.s,
        }
        try:
            self.xunit_ind = list(self.xunit_cycle.values()).index(self.specplot._xunit)
        except:
            self.xunit_ind = 1

        self.xunit_button_ax = self.canvas.figure.add_axes([hcoords[3], vcoords[1], hsize, vsize])
        self.xunit_button = Button(self.xunit_button_ax, list(self.xunit_cycle)[self.xunit_ind])
        self.xunit_button.on_clicked(self.choose_xunit)

        # Button to cycle vframe
        self.vframe_cycle = {
            "TOPO": "topo",
            "LSRK": "lsrk",
            "LSRD": "lsrd",
            "GEO": "gcrs",
            "HEL": "hcrs",
            "BARY": "icrs",
            "GAL": "galactocentric",
        }
        try:
            self.vframe_ind = list(self.vframe_cycle.values()).index(self.specplot._plot_kwargs["vel_frame"])
        except:
            self.vframe_ind = 0
        self.vframe_button_ax = self.canvas.figure.add_axes([hcoords[4], vcoords[0], hsize, vsize])
        self.vframe_button = Button(self.vframe_button_ax, list(self.vframe_cycle)[self.vframe_ind])
        self.vframe_button.on_clicked(self.choose_vframe)
        # self.vframe_radio_ax = self.canvas.figure.add_axes([hcoords[4], vcoords[1]-vsize*7, hsize*2, vsize*7],
        # zorder=100)
        # self.vframe_radio_ax.set_visible(False)

        # Button to cycle voffset
        # self.voffset_button_ax = self.canvas.figure.add_axes([hcoords[4], vcoords[1], hsize, vsize])
        # self.voffset_button = Button(self.voffset_button_ax, "Voffset")
        # self.voffset_button.on_clicked(self.open_voffset_radio)
        # self.voffset_radio_ax = self.canvas.figure.add_axes([hcoords[4], vcoords[1]-vsize*2, hsize*2, vsize*2],
        # zorder=100)
        # self.voffset_radio_ax.set_visible(False)

        # Button to cycle vdef
        self.vdef_cycle = {"RADIO": "radio", "OPTI": "optical", "RELA": "relativistic"}
        try:
            self.vdef_ind = list(self.vdef_cycle.values()).index(self.specplot._plot_kwargs["doppler_convention"])
        except:
            self.vdef_ind = 0
        self.vdef_button_ax = self.canvas.figure.add_axes([hcoords[5], vcoords[0], hsize, vsize])
        self.vdef_button = Button(self.vdef_button_ax, list(self.vdef_cycle)[self.vdef_ind])
        self.vdef_button.on_clicked(self.choose_vdef)
        # self.vdef_radio_ax = self.canvas.figure.add_axes([hcoords[5], vcoords[1]-vsize*3, hsize*2, vsize*3],
        # zorder=100)
        # self.vdef_radio_ax.set_visible(False)

        # Button/Radio combo to select center freq
        # self.centfreq_button_ax = self.canvas.figure.add_axes([hcoords[5], vcoords[1], hsize, vsize])
        # self.centfreq_button = Button(self.centfreq_button_ax, "CentFreq")
        # self.centfreq_button.on_clicked(self.open_centfreq_radio)
        # self.centfreq_radio_ax = self.canvas.figure.add_axes([hcoords[5], vcoords[1]-vsize*2, hsize*2, vsize*2],
        # zorder=100)
        # self.centfreq_radio_ax.set_visible(False)

    # def _writeascii(self, event=None):
    #     print("pass")

    def _zline_enable(self, event=None):
        print(f"zline is now {not self.specplot._zline.get_visible()}")
        self.specplot._zline.set_visible(not self.specplot._zline.get_visible())
        # self.zline_button.set_active(not self.zline_button.get_active())
        # self.specplot._axis.figure.canvas.draw_idle()
        # self.specplot._zline._axis.canvas.draw_idle()

    def _hist_enable(self, event=None):
        hist = self.specplot._line.get_drawstyle() == "steps-mid"
        print(f"hist is now {not hist}")
        if hist:
            self.specplot._line.set_drawstyle("default")
        else:
            self.specplot._line.set_drawstyle("steps-mid")
        # self.specplot._histline.set_visible(not self.specplot._histline.get_visible())
        # self.specplot._zline._axis.canvas.draw_idle()

    def _regionshow_enable(self, event=None):
        # need to have a persistent bool here in case users disable showing regions then make more
        self.regionshow = not self.regionshow
        print(f"regionshow is now {self.regionshow}")
        if self.specplot._selector is not None:
            for region in self.specplot._selector.regions:
                region.set_visible(self.regionshow)

    def _clearoverlays(self, event=None):
        print("clear overlays")
        self.specplot.clear_overlays()

    def choose_leftclick(self, event=None):
        self.xunit_ind += 1
        i = self.xunit_ind % len(self.xunit_cycle)
        self.xunit_button.label.set_text(self.xunit_cycle[i])
        # print(choice)
        # self.leftclick_radio = None
        # self.leftclick_radio_ax.set_visible(False)
        self.specplot._plt.draw()
        self.specplot._axis.figure.canvas.draw_idle()

    def choose_xunit(self, c):
        self.xunit_ind += 1
        i = self.xunit_ind % len(self.xunit_cycle)
        # self.xunit_button.label.set_text(self.xunit_cycle[i][0])
        self.xunit_button.label.set_text(list(self.xunit_cycle)[i])
        self.specplot._xunit = list(self.xunit_cycle.values())[i]
        # TODO light amount of refactoring here between choose_xunit, chose_vframe, choose_vdef etc

        if i == 0:
            self.specplot._sa = u.Quantity(np.arange(len(self.specplot._sa)))
            self.specplot._plot_kwargs["xlabel"] = "Channel"
            self.specplot.axis.set_xlabel(self.specplot._plot_kwargs["xlabel"])
            self.specplot._plot_kwargs["xaxis_unit"] = None
        else:
            self.specplot._sa = self.specplot._spectrum.velocity_axis_to(
                unit=self.specplot._xunit,
                toframe=self.specplot._plot_kwargs["vel_frame"],
                doppler_convention=self.specplot._plot_kwargs["doppler_convention"],
            )
            self.specplot._spectrum._spectral_axis = self.specplot._sa
            self.specplot._set_xaxis_info()
            self.specplot._plot_kwargs["xlabel"] = self.specplot._compose_xlabel(**self.specplot._plot_kwargs)

            self.specplot._set_labels(**self.specplot._plot_kwargs)

        self.specplot._line.set_xdata(self.specplot._sa)
        self.specplot._axis.set_xlim(np.min(self.specplot._sa).value, np.max(self.specplot._sa).value)
        self.specplot._plt.draw()
        self.specplot._axis.figure.canvas.draw_idle()

    def choose_vframe(self, choice):
        self.vframe_ind += 1
        i = self.vframe_ind % len(self.vframe_cycle)
        self.vframe_button.label.set_text(list(self.vframe_cycle)[i])
        self.specplot._plot_kwargs["vel_frame"] = list(self.vframe_cycle.values())[i]

        # if i == 0:
        #     self.specplot.spectrum.set_frame(self.gbtloc)
        # else:
        self.specplot.spectrum.set_frame(self.specplot._plot_kwargs["vel_frame"])

        self.specplot._sa = self.specplot.spectrum.spectral_axis
        self.specplot._set_xaxis_info()
        self.specplot._plot_kwargs["xlabel"] = self.specplot._compose_xlabel(**self.specplot._plot_kwargs)
        self.specplot._set_labels(**self.specplot._plot_kwargs)
        self.specplot._line.set_xdata(self.specplot._sa)
        self.specplot._axis.set_xlim(np.min(self.specplot._sa).value, np.max(self.specplot._sa).value)

        self.specplot._plt.draw()
        self.specplot._axis.figure.canvas.draw_idle()

    # def open_voffset_radio(self, event=None):
    #     print("choose voffset")
    #     self.voffset_radio_ax.set_visible(True)
    #     self.voffset_radio = RadioButtons(self.voffset_radio_ax, ("vsrc", "0"))
    #     self.voffset_radio.on_clicked(self.choose_voffset)
    #     self.specplot._plt.draw()
    #     self.specplot._axis.figure.canvas.draw_idle()

    # def choose_voffset(self, choice):
    #     print(choice)
    #     self.voffset_radio = None
    #     self.voffset_radio_ax.set_visible(False)
    #     self.specplot._plt.draw()
    #     self.specplot._axis.figure.canvas.draw_idle()

    def choose_vdef(self, choice):
        self.vdef_ind += 1
        i = self.vdef_ind % len(self.vdef_cycle)
        self.vdef_button.label.set_text(list(self.vdef_cycle)[i])
        self.specplot._plot_kwargs["doppler_convention"] = list(self.vdef_cycle.values())[i]

        self.specplot.spectrum.set_frame(self.specplot._plot_kwargs["vel_frame"])

        self.specplot._sa = self.specplot._spectrum.velocity_axis_to(
            unit=self.specplot._xunit,
            toframe=self.specplot._plot_kwargs["vel_frame"],
            doppler_convention=self.specplot._plot_kwargs["doppler_convention"],
        )
        self.specplot._spectrum._spectral_axis = self.specplot._sa
        self.specplot._set_xaxis_info()
        self.specplot._plot_kwargs["xlabel"] = self.specplot._compose_xlabel(**self.specplot._plot_kwargs)

        self.specplot._set_labels(**self.specplot._plot_kwargs)

        self.specplot._line.set_xdata(self.specplot._sa)
        self.specplot._axis.set_xlim(np.min(self.specplot._sa).value, np.max(self.specplot._sa).value)
        self.specplot._plt.draw()
        self.specplot._axis.figure.canvas.draw_idle()

    # def open_centfreq_radio(self, event=None):
    #     print("choose centfreq")
    #     self.centfreq_radio_ax.set_visible(True)
    #     self.centfreq_radio = RadioButtons(self.centfreq_radio_ax, ("abs", "rel"))
    #     self.centfreq_radio.on_clicked(self.choose_centfreq)
    #     self.specplot._plt.draw()
    #     self.specplot._axis.figure.canvas.draw_idle()

    # def choose_centfreq(self, choice):
    #     print(choice)
    #     self.centfreq_radio = None
    #     self.centfreq_radio_ax.set_visible(False)
    #     self.specplot._plt.draw()
    #     self.specplot._axis.figure.canvas.draw_idle()
