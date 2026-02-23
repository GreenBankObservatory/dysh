"""
Plot a spectrum using matplotlib
"""

import datetime as dt
import os
import sys

import astropy.units as u
import matplotlib as mpl
import numpy as np
from astropy.coordinates import SkyCoord

from dysh.log import logger

from ..coordinates import (
    Observatory,
    crval4_to_pol,
    ra2ha,
)
from ..util.core import abbreviate_to, in_notebook

# Only import the appropriate GUI to avoid errors.
on_rtd = os.environ.get("READTHEDOCS") == "True"
if in_notebook() and on_rtd:
    from IPython import get_ipython

    get_ipython().run_line_magic("matplotlib", "inline")
    from .labgui import StaticLabGUI as GUI
elif in_notebook():
    from .labgui import LabGUI as GUI
elif os.environ.get("DISPLAY", "") == "":
    from .staticgui import StaticGUI as GUI
else:
    from IPython import get_ipython

    ipython = get_ipython()
    try:
        if ipython.active_eventloop != "tk":
            msg = (
                "tk event loop not started. Plotting may not work as expected.\nTo start the tk event loop use: %gui tk"
            )
            if logger._configured:
                logger.warning(msg)
            else:
                print(msg)
    except AttributeError:
        # Do not show this message if this line is encountered while loading
        # modules as part of the dysh shell startup.
        if "bin/dysh" not in sys.argv[0]:
            msg = "Not running on IPython and trying to use the ShellGUI may result in unexpected behavior."
            if logger._configured:
                logger.warning(msg)
            else:
                print(msg)
    from .shellgui import ShellGUI as GUI


_KMS = u.km / u.s

mpl.rcParams["font.family"] = "monospace"


class PlotBase:
    """This class describes describes the common interface to all Plot classes."""

    def __init__(self):
        self._init_plot()
        self._scan_numbers = None

    def _set_frontend(self):
        self._frontend = GUI(self)
        self._connect()
        self._frontend.connect_buttons(self)

    def _connect(self):
        if self.figure.canvas is not None:
            self.figure.canvas.mpl_connect("close_event", self._close)
        else:
            return

    def _close(self, event=None):
        if self.has_selector():
            self._selector.close()
            self._selector = None
        self.figure = None
        self.axes = None

    def _init_plot(self):
        if not self.has_figure():
            self.figure = mpl.figure.Figure(figsize=(10, 6), dpi=100)
        else:
            self.figure.clear()
        self.axes = self.figure.subplots(nrows=1, ncols=1)

    def has_axes(self):
        return hasattr(self, "axes") and self.axes is not None

    def has_figure(self):
        return hasattr(self, "figure") and self.figure is not None

    def has_selector(self):
        return hasattr(self, "_selector") and self._selector is not None

    def _plot_type(self):
        """The plot object"""
        return self.__class__.__name__

    def show(self):
        self._frontend.show()

    @property
    def axis(self):
        """The underlying :class:`~matplotlib.Axes` object."""
        return self.axes

    def fmt_scans(self, abbreviate=True):
        """Format the scan numbers."""
        csl = ",".join(map(str, map(int, self._scan_numbers)))
        if abbreviate:
            csl = abbreviate_to(10, csl)
        return csl

    def _set_header(self, s):
        fsize_small = 9
        fsize_large = 14
        xyc = "figure fraction"

        hcoords = np.array([0.05, 0.21, 0.41, 0.59, 0.77])
        vcoords = np.array([0.94, 0.9, 0.86])

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
                obstime=s._obstime,
                location=Observatory.get_earth_location(s.meta["SITELONG"], s.meta["SITELAT"], s.meta["SITEELEV"]),
            )
            out_str = sc.transform_to("fk5").to_string("hmsdms", sep=" ", precision=2)[:-1]
            out_ra = out_str[:11]
            out_dec = out_str[12:]
            return out_ra, out_dec

        # col 1
        self.axis.annotate(f"Scan(s) {self.fmt_scans():>10}", (hcoords[0], vcoords[0]), xycoords=xyc, size=fsize_small)
        self.axis.annotate(f"{s.meta['DATE-OBS'][:10]}", (hcoords[0], vcoords[1]), xycoords=xyc, size=fsize_small)
        self.axis.annotate(f"{s.meta['OBSERVER']}", (hcoords[0], vcoords[2]), xycoords=xyc, size=fsize_small)

        # col 2
        ljust = 4  # Left justify the column names by this amount.
        velo = s.meta["VELOCITY"] * 1e-3  # * u.km / u.s # GBTIDL doesn't say km/s so neither will I (saves space)
        self.axis.annotate(
            f"{'V':<{ljust}}: {velo} {s.meta['VELDEF']}", (hcoords[1], vcoords[0]), xycoords=xyc, size=fsize_small
        )
        self.axis.annotate(
            f"{'Int':<{ljust}}: {time_formatter(s.meta['EXPOSURE'])}",
            (hcoords[1], vcoords[1]),
            xycoords=xyc,
            size=fsize_small,
        )
        self.axis.annotate(
            f"{'LST':<{ljust}}: {time_formatter(s.meta['LST'])}",
            (hcoords[1], vcoords[2]),
            xycoords=xyc,
            size=fsize_small,
        )

        # col 3
        # TODO: need to understand frequencies to assign correct title
        # instead of just forcing to GHz with 5 decimal points
        ljust = 5
        f0 = np.around(s.rest_value.to(u.GHz), 5)
        self.axis.annotate(f"{'F0':<{ljust}}:  {f0}", (hcoords[2], vcoords[0]), xycoords=xyc, size=fsize_small)
        fsky = np.around(s.meta["OBSFREQ"] * 1e-9, 5) * u.GHz  # or CRVAL1?
        self.axis.annotate(f"{'Fsky':<{ljust}}:  {fsky}", (hcoords[2], vcoords[1]), xycoords=xyc, size=fsize_small)
        bw = np.around(s.meta["BANDWID"] * 1e-6, 4) * u.MHz
        self.axis.annotate(f"{'BW':<{ljust}}:  {bw}", (hcoords[2], vcoords[2]), xycoords=xyc, size=fsize_small)

        # col 4
        ljust = 4
        self.axis.annotate(
            f"{'Pol':<{ljust}}:   {crval4_to_pol[s.meta['CRVAL4']]}",
            (hcoords[3], vcoords[0]),
            xycoords=xyc,
            size=fsize_small,
        )
        self.axis.annotate(
            f"{'IF':<{ljust}}:    {s.meta['IFNUM']}", (hcoords[3], vcoords[1]), xycoords=xyc, size=fsize_small
        )
        self.axis.annotate(f"{s.meta['PROJID']}", (hcoords[3], vcoords[2]), xycoords=xyc, size=fsize_small)

        # col 5
        _tsys = np.around(s.meta["TSYS"], 2)
        self.axis.annotate(f"Tsys   :  {_tsys}", (hcoords[4], vcoords[0]), xycoords=xyc, size=fsize_small)
        self.axis.annotate(
            f"Tcal   :  {np.around(s.meta['TCAL'], 2)}", (hcoords[4], vcoords[1]), xycoords=xyc, size=fsize_small
        )
        self.axis.annotate(f"{s.meta['PROC']}", (hcoords[4], vcoords[2]), xycoords=xyc, size=fsize_small)

        # bottom row
        vcoord_bot = 0.82
        hcoord_bot = 0.95
        ra, dec = coord_formatter(s)
        self.axis.annotate(f"{ra}  {dec}", (hcoords[0], vcoord_bot), xycoords=xyc, size=fsize_small)
        if self.axis.get_title() == "":
            self.axis.annotate(
                f"{s.meta['OBJECT']}", (0.5, vcoord_bot), xycoords=xyc, size=fsize_large, horizontalalignment="center"
            )

        az = np.around(s.meta["AZIMUTH"], 1)
        el = np.around(s.meta["ELEVATIO"], 1)
        ha = ra2ha(s.meta["LST"], s.meta["CRVAL2"])
        self.axis.annotate(
            f"Az: {az}  El: {el}  HA: {ha}",
            (hcoord_bot, vcoord_bot),
            xycoords=xyc,
            size=fsize_small,
            horizontalalignment="right",
        )

        # last corner -- current date time.
        ts = str(dt.datetime.now())[:19]
        self.axis.annotate(
            f"{ts}", (hcoord_bot - 0.1, 0.01), xycoords=xyc, size=fsize_small, horizontalalignment="right"
        )

    def refresh(self):
        """Refresh the plot"""
        if self.axes is not None:
            self.axes.figure.canvas.draw_idle()

    def savefig(self, file, hidebuttons=True, **kwargs):
        r"""Save the plot

        Parameters
        ----------
        file : str
            The output file name.
        hidebuttons : bool
            Hide the buttons in the output.
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
        if hidebuttons:
            for button in self.figure._localaxes[1:]:
                button.set_visible(False)
            self.figure.savefig(file, *kwargs)
            for button in self.figure._localaxes[1:]:
                button.set_visible(True)
        else:
            self.figure.savefig(file, *kwargs)
