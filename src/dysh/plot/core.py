"""
Plot a spectrum using matplotlib
"""

import datetime as dt

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord

from dysh.log import logger

from ..coordinates import (
    Observatory,
    crval4_to_pol,
    ra2ha,
)
from ..util.core import abbreviate_to

_KMS = u.km / u.s


class PlotBase:
    """This class describes describes the common interface to all Plot classes."""

    def __init__(self, **kwargs):
        self.reset()
        self._figure = None
        self._plt = plt
        self._plt.rcParams["font.family"] = "monospace"
        self._scan_numbers = None

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

    def fmt_scans(self, abbreviate=True):
        """Format the scan numbers."""
        csl = ",".join(map(str, map(int, self._scan_numbers)))
        if abbreviate:
            csl = abbreviate_to(10, csl)
        return csl

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
                obstime=s._obstime,
                location=Observatory.get_earth_location(s.meta["SITELONG"], s.meta["SITELAT"], s.meta["SITEELEV"]),
            )
            out_str = sc.transform_to("fk5").to_string("hmsdms", sep=" ", precision=2)[:-1]
            out_ra = out_str[:11]
            out_dec = out_str[12:]
            return out_ra, out_dec

        # col 1
        self._axis.annotate(f"Scan(s) {self.fmt_scans():>10}", (hcoords[0], vcoords[0]), xycoords=xyc, size=fsize_small)
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
        f0 = np.around(s.rest_value.to(u.GHz), 5)
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


def check_kwargs(known_kwargs, kwargs):
    """Check if `kwargs` are in `known_kwargs`"""
    diff = set(kwargs) - set(known_kwargs)
    if len(diff) > 0:
        logger.warning(f"Unknown kwargs: {', '.join(diff)}")
