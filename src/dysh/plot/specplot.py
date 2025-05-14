"""
Plot a spectrum using matplotlib
"""

import datetime as dt
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
        self.reset()
        self._spectrum = spectrum
        self._set_xaxis_info()
        self._plot_kwargs.update(kwargs)
        self._plt = plt
        self._figure = None
        self._axis = None
        self._title = self._plot_kwargs["title"]
        self._selector: InteractiveSpanSelector = None

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

    def plot(self, show_header=True, select=True, show=True, **kwargs):
        # @todo document kwargs here
        r"""
        Plot the spectrum.

        Parameters
        ----------
        show_header : bool
            Show informational header in the style of GBTIDL, default: True.
        **kwargs : various
            keyword=value arguments (need to describe these in a central place)
        """
        if show:
            plt.ion()
        plt.rcParams["font.family"] = "monospace"
        # plt.rcParams['axes.formatter.useoffset'] = False # Disable use of offset.

        # xtype = 'velocity, 'frequency', 'wavelength'
        # if self._figure is None:
        self._set_xaxis_info()
        # plot arguments for this call of plot(). i.e. non-sticky plot attributes
        this_plot_kwargs = deepcopy(self._plot_kwargs)
        this_plot_kwargs.update(kwargs)
        if True:  # @todo deal with plot reuse (notebook vs script)
            self._figure, self._axis = self._plt.subplots(figsize=(10, 6))

        # else:
        #    self._axis.cla()
        def apply_region_selection(x, y):  # or list of start/stop values?
            """Apply region selected using Selection"""
            # sdf.Select(blah blah blah)
            print(x, y)

        # TODO: procedurally generate subplot params based on show header/buttons args.
        # ideally place left/right params right here, then top gets determined below.

        s = self._spectrum
        if show_header:
            self._figure.subplots_adjust(top=0.7, left=0.09, right=0.95)
            self._set_header(s)

            # callback = Index()
            # axtest = self._figure.add_axes([0.1, 0.9, 0.1, 0.075])
            # self._btest = Button(axtest, 'Test')
            # self._btest.on_clicked(self.next)

        # if select:
        #     self.setregion(sa)

        self._sa = s.spectral_axis
        lw = this_plot_kwargs["linewidth"]
        xunit = this_plot_kwargs["xaxis_unit"]
        yunit = this_plot_kwargs["yaxis_unit"]
        if xunit is None:
            xunit = str(sa.unit)  # noqa: F821
        if "vel_frame" not in this_plot_kwargs:
            if u.Unit(xunit).is_equivalent("km/s") and "VELDEF" in s.meta:
                # If the user specified velocity units, default to
                # the velframe the data were taken in.  This we can
                # get from VELDEF keyword.  See issue #303
                this_plot_kwargs["vel_frame"] = decode_veldef(s.meta["VELDEF"])[1].lower()
            else:
                this_plot_kwargs["vel_frame"] = s.velocity_frame
        if "chan" in str(xunit).lower():
            self._sa = u.Quantity(np.arange(len(self._sa)))
            this_plot_kwargs["xlabel"] = "Channel"
        else:
            # convert the x axis to the requested
            # print(f"EQUIV {equiv} doppler_rest {sa.doppler_rest} [{rfq}] convention {convention}")
            # sa = s.spectral_axis.to( self._plot_kwargs["xaxis_unit"],
            #   equivalencies=equiv,doppler_rest=rfq, doppler_convention=convention)
            self._sa = s.velocity_axis_to(
                unit=xunit,
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
        self._axis.tick_params(axis="both", which="both", bottom=True, top=True, left=True, right=True, direction="in")
        if this_plot_kwargs["grid"]:
            self._axis.grid(visible=True, which="major", axis="both", lw=lw / 2, color="k", alpha=0.33)
            self._axis.grid(visible=True, which="minor", axis="both", lw=lw / 2, color="k", alpha=0.22, linestyle="--")

        self._set_labels(**this_plot_kwargs)
        # self._axis.axhline(y=0,color='red',lw=2)
        if self._title is not None:
            self._axis.set_title(self._title)

        if select:
            self._selector = InteractiveSpanSelector(self._axis)
            self._spectrum._selection = self._selector.get_selected_regions()

        if show:
            self.refresh()

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

    def _set_header(self, s):
        fsize_small = 9
        fsize_large = 14
        xyc = "figure fraction"

        hcoords = np.array([0.05, 0.21, 0.41, 0.59, 0.77])
        vcoords = np.array([0.84, 0.8, 0.76])

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
        velo = s.meta["VELOCITY"] * 1e-3  # * u.km / u.s # GBTIDL doesn't say km/s so neither will I (saves space)
        self._axis.annotate(
            f"V   : {velo} {s.meta['VELDEF']}", (hcoords[1], vcoords[0]), xycoords=xyc, size=fsize_small
        )
        self._axis.annotate(
            f"Int : {time_formatter(s.meta['DURATION'])}", (hcoords[1], vcoords[1]), xycoords=xyc, size=fsize_small
        )
        self._axis.annotate(
            f"LST : {time_formatter(s.meta['LST'])}", (hcoords[1], vcoords[2]), xycoords=xyc, size=fsize_small
        )

        # col 3
        # TODO: need to understand frequencies to assign correct title
        # instead of just forcing to GHz with 5 decimal points
        f0 = np.around(s.meta["RESTFREQ"] * 1e-9, 5) * u.GHz
        self._axis.annotate(f"F0   :  {f0}", (hcoords[2], vcoords[0]), xycoords=xyc, size=fsize_small)
        fsky = np.around(s.meta["OBSFREQ"] * 1e-9, 5) * u.GHz  # or CRVAL1?
        self._axis.annotate(f"Fsky :  {fsky}", (hcoords[2], vcoords[1]), xycoords=xyc, size=fsize_small)
        bw = np.around(s.meta["BANDWID"] * 1e-6, 4) * u.MHz
        self._axis.annotate(f"BW   :  {bw}", (hcoords[2], vcoords[2]), xycoords=xyc, size=fsize_small)

        # col 4
        self._axis.annotate(
            f"Pol  :   {crval4_to_pol[s.meta['CRVAL4']]}", (hcoords[3], vcoords[0]), xycoords=xyc, size=fsize_small
        )
        self._axis.annotate(f"IF   :    {s.meta['IFNUM']}", (hcoords[3], vcoords[1]), xycoords=xyc, size=fsize_small)
        self._axis.annotate(f"{s.meta['PROJID']}", (hcoords[3], vcoords[2]), xycoords=xyc, size=fsize_small)

        # col 5
        self._axis.annotate(
            f"Tsys   :  {np.around(s.meta['MEANTSYS'], 2)}", (hcoords[4], vcoords[0]), xycoords=xyc, size=fsize_small
        )
        self._axis.annotate(
            f"Tcal   :  {np.around(s.meta['TCAL'], 2)}", (hcoords[4], vcoords[1]), xycoords=xyc, size=fsize_small
        )
        self._axis.annotate(f"{s.meta['PROC']}", (hcoords[4], vcoords[2]), xycoords=xyc, size=fsize_small)

        # bottom row
        ra, dec = coord_formatter(s)
        self._axis.annotate(f"{ra}  {dec}", (hcoords[0], 0.71), xycoords=xyc, size=fsize_small)
        self._axis.annotate(
            f"{s.meta['OBJECT']}", (0.5, 0.71), xycoords=xyc, size=fsize_large, horizontalalignment="center"
        )
        az = np.around(s.meta["AZIMUTH"], 1)
        el = np.around(s.meta["ELEVATIO"], 1)
        ha = ra2ha(s.meta["LST"], s.meta["CRVAL2"])
        self._axis.annotate(
            f"Az: {az}  El: {el}  HA: {ha}", (0.95, 0.71), xycoords=xyc, size=fsize_small, horizontalalignment="right"
        )

        # bottom row
        ra, dec = coord_formatter(s)
        self._axis.annotate(f"{ra}  {dec}", (hcoords[0], 0.71), xycoords=xyc, size=fsize_small)
        self._axis.annotate(
            f"{s.meta['OBJECT']}", (0.5, 0.71), xycoords=xyc, size=fsize_large, horizontalalignment="center"
        )
        az = np.around(s.meta["AZIMUTH"], 1)
        el = np.around(s.meta["ELEVATIO"], 1)
        ha = ra2ha(s.meta["LST"], s.meta["CRVAL2"])
        self._axis.annotate(
            f"Az: {az}  El: {el}  HA: {ha}", (0.95, 0.71), xycoords=xyc, size=fsize_small, horizontalalignment="right"
        )

        # last corner
        ts = str(dt.datetime.now())[:19]
        self._axis.annotate(f"{ts}", (0.85, 0.01), xycoords=xyc, size=fsize_small, horizontalalignment="right")

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
        self.figure.savefig(file, *kwargs)

    def get_selected_regions(self):
        """ """
        regions = self._selector.get_selected_regions()
        return [tuple(np.sort([np.argmin(abs(p - self._sa.value)) for p in r])) for r in regions]


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
        self.button_ax = self.canvas.figure.add_axes([0.1, 0.025, 0.12, 0.04])
        self.clear_button = Button(self.button_ax, "Clear Regions")
        self.clear_button.on_clicked(self.clear_regions)

        # Button to clear a single region.
        self.button2_ax = self.canvas.figure.add_axes([0.24, 0.025, 0.12, 0.04])
        self.del_button = Button(self.button2_ax, "Delete Region")
        self.del_button.on_clicked(self.clear_region)

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
