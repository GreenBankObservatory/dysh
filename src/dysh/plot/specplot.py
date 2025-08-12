"""
Plot a spectrum using matplotlib
"""

import datetime as dt
import warnings
from copy import deepcopy

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.utils.masked import Masked
from matplotlib.patches import Rectangle
from matplotlib.widgets import Button, SpanSelector

from ..coordinates import (
    Observatory,
    crval4_to_pol,
    decode_veldef,
    frame_to_label,
    ra2ha,
)

_KMS = u.km / u.s


class PlotBase:
    """This class describes describes the common interface to all Plot classes."""

    def __init__(self, **kwargs):
        self.reset()
        self._figure = None
        self._axis = None
        self._plt = plt
        self._plt.rcParams["font.family"] = "monospace"

    def _plot_type(self):
        """The plot object"""
        return self.__class__.__name__

    @property
    def axis(self):
        """The underlying :class:`~matplotlib.Axes` object"""
        return self._axis

    @property
    def figure(self):
        """The underlying :class:`~matplotlib.Figure` object"""
        return self._figure

    def reset(self):
        """Reset the plot keyword arguments to their defaults."""
        if self._plot_type() == "SpectrumPlot":
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
        elif self._plot_type() == "ScanPlot":
            self._plot_kwargs = {
                "title": None,
                "cmap": "inferno",
                "interpolation": "nearest",
            }

    def _set_header(self, s):
        move_vcoords_bool = 0.1 * (self._plot_type() == "ScanPlot")
        fsize_small = 9
        fsize_large = 14
        xyc = "figure fraction"

        hcoords = np.array([0.05, 0.21, 0.41, 0.59, 0.77])
        vcoords = np.array([0.84, 0.8, 0.76]) + move_vcoords_bool

        def time_formatter(time_sec):
            hh = int(time_sec // 3600)
            mm = int((time_sec - 3600 * hh) // 60)
            ss = np.around((time_sec - 3600 * hh - 60 * mm), 1)
            return f"{str(hh).zfill(2)} {str(mm).zfill(2)} {str(ss).zfill(3)}"

        def coord_formatter(s):
            sc = SkyCoord(
                s.meta["CRVAL2"],
                s.meta["CRVAL3"],
                unit="deg",
                frame=s.meta["RADESYS"].lower(),
                obstime=Time(s.meta["DATE-OBS"]),
                location=Observatory.get_earth_location(s.meta["SITELONG"], s.meta["SITELAT"], s.meta["SITEELEV"]),
            )
            out_str = sc.transform_to("fk5").to_string("hmsdms", sep=" ", precision=2)[:-1]
            out_ra = out_str[:11]
            out_dec = out_str[12:]
            return out_ra, out_dec

        # col 1
        self._axis.annotate(f"Scan     {s.meta['SCAN']}", (hcoords[0], vcoords[0]), xycoords=xyc, size=fsize_small)
        self._axis.annotate(f"{s.meta['DATE-OBS'][:10]}", (hcoords[0], vcoords[1]), xycoords=xyc, size=fsize_small)
        self._axis.annotate(f"{s.meta['OBSERVER']}", (hcoords[0], vcoords[2]), xycoords=xyc, size=fsize_small)

        # col 2
        ljust = 4  # Left justify the column names by this amount.
        velo = s.meta["VELOCITY"] * 1e-3  # * u.km / u.s # GBTIDL doesn't say km/s so neither will I (saves space)
        self._axis.annotate(
            f"{'V':<{ljust}}: {velo} {s.meta['VELDEF']}", (hcoords[1], vcoords[0]), xycoords=xyc, size=fsize_small
        )
        self._axis.annotate(
            f"{'Int':<{ljust}}: {time_formatter(s.meta['EXPOSURE'])}",
            (hcoords[1], vcoords[1]),
            xycoords=xyc,
            size=fsize_small,
        )
        self._axis.annotate(
            f"{'LST':<{ljust}}: {time_formatter(s.meta['LST'])}",
            (hcoords[1], vcoords[2]),
            xycoords=xyc,
            size=fsize_small,
        )

        # col 3
        # TODO: need to understand frequencies to assign correct title
        # instead of just forcing to GHz with 5 decimal points
        ljust = 5
        f0 = np.around(s.meta["RESTFREQ"] * 1e-9, 5) * u.GHz
        self._axis.annotate(f"{'F0':<{ljust}}:  {f0}", (hcoords[2], vcoords[0]), xycoords=xyc, size=fsize_small)
        fsky = np.around(s.meta["OBSFREQ"] * 1e-9, 5) * u.GHz  # or CRVAL1?
        self._axis.annotate(f"{'Fsky':<{ljust}}:  {fsky}", (hcoords[2], vcoords[1]), xycoords=xyc, size=fsize_small)
        bw = np.around(s.meta["BANDWID"] * 1e-6, 4) * u.MHz
        self._axis.annotate(f"{'BW':<{ljust}}:  {bw}", (hcoords[2], vcoords[2]), xycoords=xyc, size=fsize_small)

        # col 4
        ljust = 4
        self._axis.annotate(
            f"{'Pol':<{ljust}}:   {crval4_to_pol[s.meta['CRVAL4']]}",
            (hcoords[3], vcoords[0]),
            xycoords=xyc,
            size=fsize_small,
        )
        self._axis.annotate(
            f"{'IF':<{ljust}}:    {s.meta['IFNUM']}", (hcoords[3], vcoords[1]), xycoords=xyc, size=fsize_small
        )
        self._axis.annotate(f"{s.meta['PROJID']}", (hcoords[3], vcoords[2]), xycoords=xyc, size=fsize_small)

        # col 5
        _tsys = np.around(s.meta["TSYS"], 2)
        self._axis.annotate(f"Tsys   :  {_tsys}", (hcoords[4], vcoords[0]), xycoords=xyc, size=fsize_small)
        self._axis.annotate(
            f"Tcal   :  {np.around(s.meta['TCAL'], 2)}", (hcoords[4], vcoords[1]), xycoords=xyc, size=fsize_small
        )
        self._axis.annotate(f"{s.meta['PROC']}", (hcoords[4], vcoords[2]), xycoords=xyc, size=fsize_small)

        # bottom row
        vcoord_bot = 0.72 + move_vcoords_bool
        hcoord_bot = 0.95
        ra, dec = coord_formatter(s)
        self._axis.annotate(f"{ra}  {dec}", (hcoords[0], vcoord_bot), xycoords=xyc, size=fsize_small)
        if self._axis.get_title() == "":
            self._axis.annotate(
                f"{s.meta['OBJECT']}", (0.5, vcoord_bot), xycoords=xyc, size=fsize_large, horizontalalignment="center"
            )
        az = np.around(s.meta["AZIMUTH"], 1)
        el = np.around(s.meta["ELEVATIO"], 1)
        ha = ra2ha(s.meta["LST"], s.meta["CRVAL2"])
        self._axis.annotate(
            f"Az: {az}  El: {el}  HA: {ha}",
            (hcoord_bot, vcoord_bot),
            xycoords=xyc,
            size=fsize_small,
            horizontalalignment="right",
        )

        # last corner -- current date time.
        ts = str(dt.datetime.now())[:19]
        self._axis.annotate(
            f"{ts}", (hcoord_bot - 0.1, 0.01), xycoords=xyc, size=fsize_small, horizontalalignment="right"
        )

    def refresh(self):
        """Refresh the plot"""
        if self.axis is not None:
            self.axis.figure.canvas.draw_idle()

    def savefig(self, file, **kwargs):
        r"""Save the plot

        Parameters
        ----------
        file - str
            The output file name
        **kwargs : dict or key=value pairs
            Other arguments to pass to `~matplotlib.pyplot.savefig`

        """
        # TODO: add clause about cutting off the top of the figure where the interactive buttons are
        # bbox_inches = matplotlib.transforms.Bbox((0,0,10,hgt)) (warn: 10 is hardcoded in specplot)
        # or, set_visible to False
        # buttons are currently listed in the _localaxes, but this includes the plot window at index 0
        # so if the plot window ever goes missing, check the order in this list
        # there has to be a better way to do this
        # TODO: put buttons in a sub/different axes so we only have to hide the axes object instead of
        # a list of all the buttons and plots
        for button in self.figure._localaxes[1:]:
            button.set_visible(False)
        self.figure.savefig(file, *kwargs)
        for button in self.figure._localaxes[1:]:
            button.set_visible(True)


class SpectrumPlot(PlotBase):
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
        super().__init__()
        self._spectrum = spectrum
        self._sa = spectrum._spectral_axis
        self._set_xaxis_info()
        self._plot_kwargs.update(kwargs)
        self._title = self._plot_kwargs["title"]
        self._selector: InteractiveSpanSelector = None
        self._freezey = (self._plot_kwargs["ymin"] is not None) or (self._plot_kwargs["ymax"] is not None)
        self._freezex = (self._plot_kwargs["xmin"] is not None) or (self._plot_kwargs["xmax"] is not None)

    # def __call__ (see pyspeckit)

    def _set_xaxis_info(self):
        """Ensure the xaxis info is up to date if say, the spectrum frame has changed."""
        self._plot_kwargs["doppler_convention"] = self._spectrum.doppler_convention
        self._plot_kwargs["vel_frame"] = self._spectrum.velocity_frame
        self._plot_kwargs["xaxis_unit"] = self._spectrum.spectral_axis.unit
        self._plot_kwargs["yaxis_unit"] = self._spectrum.unit

    @property
    def spectrum(self):
        """The underlying `~spectra.spectrum.Spectrum`"""
        return self._spectrum

    def plot(self, show_header=True, select=True, **kwargs):
        # @todo document kwargs here
        r"""
        Plot the spectrum.

        Parameters
        ----------
        show_header : bool
            Show informational header in the style of GBTIDL, default: True.
        select : bool
            Allow selecting regions via click and drag for baseline computation, default: True
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
            self._figure, self._axis = self._plt.subplots(figsize=(10, 6))

        # TODO: procedurally generate subplot params based on show header/buttons args.
        # ideally place left/right params right here, then top gets determined below.

        s = self._spectrum

        lw = this_plot_kwargs["linewidth"]
        self._xunit = this_plot_kwargs["xaxis_unit"]  # need to kick back a ref to xunit for baseline overlays
        yunit = this_plot_kwargs["yaxis_unit"]
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
            # convert the x axis to the requested
            # print(f"EQUIV {equiv} doppler_rest {sa.doppler_rest} [{rfq}] convention {convention}")
            # sa = s.spectral_axis.to( self._plot_kwargs["xaxis_unit"],
            #   equivalencies=equiv,doppler_rest=rfq, doppler_convention=convention)
            self._sa = s.velocity_axis_to(
                unit=self._xunit,
                toframe=this_plot_kwargs["vel_frame"],
                doppler_convention=this_plot_kwargs["doppler_convention"],
            )

        sf = s.flux
        if yunit is not None:
            sf = s.flux.to(yunit)
        sf = Masked(sf, s.mask)
        lines = self._axis.plot(self._sa, sf, color=this_plot_kwargs["color"], lw=lw)
        self._line = lines[0]

        if not this_plot_kwargs["xmin"] and not this_plot_kwargs["xmax"]:
            self._axis.set_xlim(np.min(self._sa).value, np.max(self._sa).value)
        else:
            self._axis.set_xlim(this_plot_kwargs["xmin"], this_plot_kwargs["xmax"])
        self._axis.set_ylim(this_plot_kwargs["ymin"], this_plot_kwargs["ymax"])
        if self._freezey:
            self._axis.autoscale(enable=False)
        else:
            self._axis.autoscale(axis="y", enable=True)
        self._axis.tick_params(axis="both", which="both", bottom=True, top=True, left=True, right=True, direction="in")
        if this_plot_kwargs["grid"]:
            self._axis.grid(visible=True, which="major", axis="both", lw=lw / 2, color="k", alpha=0.33)
            self._axis.grid(visible=True, which="minor", axis="both", lw=lw / 2, color="k", alpha=0.22, linestyle="--")

        self._set_labels(**this_plot_kwargs)
        # self._axis.axhline(y=0,color='red',lw=2)
        if self._title is not None:
            self._axis.set_title(self._title)

        if show_header:
            self._figure.subplots_adjust(top=0.7, left=0.09, right=0.95)
            self._set_header(s)

        if select:
            self._selector = InteractiveSpanSelector(self._axis)
            self._spectrum._selection = self._selector.get_selected_regions()

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
        xlabel = kwargs.get("xlabel", None)  # noqa: F841
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

    def savefig(self, file, **kwargs):
        r"""Save the plot

        Parameters
        ----------
        file - str
            The output file name
        **kwargs : dict or key=value pairs
            Other arguments to pass to `~matplotlib.pyplot.savefig`

        """
        # TODO: add clause about cutting off the top of the figure where the interactive buttons are
        # bbox_inches = matplotlib.transforms.Bbox((0,0,10,hgt)) (warn: 10 is hardcoded in specplot)
        # or, set_visible to False
        # buttons are currently listed in the _localaxes, but this includes the plot window at index 0
        # so if the plot window ever goes missing, check the order in this list
        # there has to be a better way to do this
        # TODO: put buttons in a sub/different axes so we only have to hide the axes object instead of
        # a list of all the buttons and plots
        for button in self.figure._localaxes[1:]:
            button.set_visible(False)
        self.figure.savefig(file, *kwargs)
        for button in self.figure._localaxes[1:]:
            button.set_visible(True)

    def get_selected_regions(self):
        """ """
        regions = self._selector.get_selected_regions()
        return [tuple(np.sort([np.argmin(abs(p - self._sa.value)) for p in r])) for r in regions]

    def freex(self):
        """ "Free the X-axis if limits have been set. Resets the limits to be the span of the spectrum."""
        self._freezex = False
        # This line (and the other in specplot.py) will have to be addressed when we
        # implement multiple IF windows in the same plot
        self._axis.set_xlim(self._sa.min.value, self._sa.max.value)

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

    def clear_overlays(self, blines=True):
        """Clear Overlays from the plot.

        Parameters
        ----------
        blines : bool
            Remove only baseline models overlaid on the plot. Default: True
        """
        # clear baseline models
        if blines:
            for b in self._axis.lines:
                if b.get_gid() == "baseline":
                    b.remove()


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
        got_one = False
        for patch in self.regions:
            contains, attr = patch.contains(event)
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


class ScanPlot(PlotBase):
    r"""
    The ScanPlot class is for simple plotting of a `~scan.Scan` or `~scan.ScanBlock`
    using matplotlib functions. Plots attributes are modified using keywords
    (\*\*kwargs) described below SpectrumPlot will attempt to make smart default
    choices for the plot if no additional keywords are given.

    Parameters
    ----------
    scanblock_or_scan : `~spectra.scan.Scan` or `~spectra.scan.ScanBlock`
        The scan or scanblock to plot.
    **kwargs : dict
        Plot attribute keyword arguments, see below.

    Other Parameters
    ----------------
    spectral_unit : str
        The units to use on the frequency axis. Can be 'MHz' or 'GHz'.
    """

    def __init__(self, scanblock_or_scan, **kwargs):
        super().__init__()
        self._scanblock_or_scan = scanblock_or_scan
        self._plot_kwargs.update(kwargs)
        self._axis2 = None
        # self._title = self._plot_kwargs["title"]# todo: deal with when refactoring
        acceptable_types = ["PSScan", "TPScan", "NodScan", "FSScan", "SubBeamNodScan"]
        self._spectrum = self._scanblock_or_scan.timeaverage()
        self._sa = self._spectrum.spectral_axis

        # determine if input is a ScanBlock or a ScanBase (raise exception if neither)
        self._type = scanblock_or_scan.__class__.__name__
        if self._type == "ScanBlock":
            self._scanblock = scanblock_or_scan
            self._num_scans = len(self._scanblock)
        elif self._type in acceptable_types:
            self._scan = scanblock_or_scan
        else:
            raise Exception(f"Plotter input {self._type} does not appear to be a valid input object type")

        # handle scanblocks
        if self._type == "ScanBlock":
            self._scan_nos = []  # scan numbers in the scan block
            self._nint_nos = []  # number of integrations in each scan
            self._timestamps = []  # 0-indexed timestamps in sec for every integration
            for i, scan in enumerate(self._scanblock):
                if i == 0:
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
            self.spectrogram = self._scan._calibrated
            self._scan_nos = self._scan.scan

        self.spectrogram = self.spectrogram.T

    def plot(self, spectral_unit=None, **kwargs):
        r"""
        Plot the scan.

        Parameters
        ----------
        spectral_unit : `~astropy.unit.Unit`
            The units to use on the frequency axis. Default: MHz if below 1 GHz, GHz if above.
        **kwargs : various
            keyword=value arguments (need to describe these in a central place)
        """

        this_plot_kwargs = deepcopy(self._plot_kwargs)
        this_plot_kwargs.update(kwargs)

        cmap = kwargs.get("cmap", "inferno")
        interpolation = kwargs.get("interpolation", "nearest")

        if True:
            self._figure, self._axis = self._plt.subplots(figsize=(10, 6))
            self._axis2 = self._axis.twinx()

        self._figure.subplots_adjust(top=0.79, left=0.1, right=1.05)
        self._set_header(self._spectrum)

        # self._axis.tick_params(axis='both',direction='inout',length=8,top=False,right=False,pad=2)
        # self._axis.yaxis.set_label_position('left')
        # self._axis.yaxis.set_ticks_position('left')

        self.im = self._axis.imshow(self.spectrogram, aspect="auto", cmap=cmap, interpolation=interpolation)

        # second "plot" to get different scales on x2, y2 axes
        # self._axis2.tick_params(axis='both',direction='inout',
        #     length=8,bottom=False,left=False,top=True,right=True,pad=2)
        # self._axis2.yaxis.set_label_position('right')
        # self._axis2.yaxis.set_ticks_position('right')
        if spectral_unit is not None:
            self._sa = self._sa.to(spectral_unit)
        else:
            if self._sa[0] / (u.GHz) < 1:
                self._sa = self._sa.to(u.MHz)
            else:
                self._sa = self._sa.to(u.GHz)
        stop = self.spectrogram.shape[1]
        step = self.spectrogram.shape[1] / self.spectrogram.shape[0]
        im2 = self._axis2.plot(np.arange(0, stop, step), self._sa, linewidth=0)  # noqa: F841
        self._axis2.set_ylim((np.min(self._sa).value, np.max(self._sa).value))

        self._axis.set_xlim(0, stop - 0.5)
        self._set_labels()

    def _set_labels(self):
        # x1: bottom
        # x2: top
        # y1: left
        # y2: right
        # z: colorbar

        x1_label = "Integration"
        self._axis.set_xlabel(x1_label)

        y1_label = "Channel"
        self._axis.set_ylabel(y1_label)

        y2_unit = self._sa.unit
        if y2_unit.is_equivalent(u.Hz):
            nu = r"$\nu$"
            y2_label = f"{nu} ({y2_unit})"
        self._axis2.set_ylabel(y2_label)

        z_unit = self._spectrum.unit
        if z_unit.is_equivalent(u.K):
            z_label = f"$T_A$ ({z_unit})"
        elif z_unit.is_equivalent(u.Jy):
            snu = r"$S_{\nu}$"
            z_label = f"{snu} ({z_unit})"
        elif z_unit.is_equivalent(u.ct):
            z_label = "Counts"
        else:
            warnings.warn("Flux units are unknown", stacklevel=2)
            z_label = ""
        self._colorbar = self._figure.colorbar(self.im, label=z_label, pad=0.1)
        # matplotlib won't set this before the Figure is drawn.
        self._figure.draw_without_rendering()
        # If there's an offset, add it to the label and make the offset invisible.
        if self._colorbar.ax.yaxis.offsetText.get_text() != "":
            off = self._colorbar.ax.yaxis.offsetText.get_text()
            e = off.split("e")[1]
            self._colorbar.set_label(z_label + rf"($\times10^{{{e}}}$)")
            self._colorbar.ax.yaxis.offsetText.set_visible(False)

    def set_clim(self, vmin, vmax):
        """
        Set the vmin and vmax parameters of the image.

        Parameters
        ----------
        vmin : float
            The minimum value of the color scale.
        vmax : float
            The maximum value of the color scale.
        """
        self.im.set_clim(vmin=vmin, vmax=vmax)

    def set_interpolation(self, interpolation="nearest"):
        """
        Set the interpolation of the image.

        Parameters
        ----------
        interpolation : str
            Interpolation method. Default: "nearest".
        """
        self.im.set_interpolation(interpolation)

    def set_cmap(self, cmap="inferno"):
        """
        Set the cmap of the image.

        Parameters
        ----------
        cmap : str
            cmap used for the color scale. Default: "inferno".
        """
        self.im.set_cmap(cmap)
