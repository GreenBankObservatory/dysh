"""
Plot spectrograms from a ScanBlock using matplotlib
"""

import datetime as dt
from copy import deepcopy

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.time import Time

from ..coordinates import (
    Observatory,
    crval4_to_pol,
    ra2ha,
)

_KMS = u.km / u.s


class ScanPlot:
    """hello"""

    def __init__(self, scanblock_or_scan, **kwargs):
        self.reset()
        self._scanblock_or_scan = scanblock_or_scan
        self._plot_kwargs.update(kwargs)
        self._plt = plt
        self._figure = None
        self._axis = None
        self._axis2 = None
        self._title = self._plot_kwargs["title"]
        acceptable_types = ["PSScan", "TPScan", "NodScan", "FSScan", "SubBeamNodScan"]

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
            self._spectral_axis = self._scanblock[0].timeaverage().spectral_axis
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
        self._s = self._scanblock_or_scan.timeaverage()

    def plot(self, **kwargs):
        """hi hello"""

        self.__init__(self._scanblock_or_scan, **kwargs)
        plt.ion()

        # self._set_xaxis_info()
        this_plot_kwargs = deepcopy(self._plot_kwargs)
        this_plot_kwargs.update(kwargs)

        cmap = kwargs.get("cmap", "inferno")
        interpolation = kwargs.get("interpolation", "nearest")

        if True:
            self._figure, self._axis = self._plt.subplots(figsize=(10, 6))
            self._axis2 = self._axis.twinx()

        self._figure.subplots_adjust(top=0.79, left=0.1, right=1.05)
        self._set_header(self._s)

        # self._axis.tick_params(axis='both',direction='inout',length=8,top=False,right=False,pad=2)
        # self._axis.yaxis.set_label_position('left')
        # self._axis.yaxis.set_ticks_position('left')

        self.im = self._axis.imshow(self.spectrogram, aspect="auto", cmap=cmap, interpolation=interpolation)

        # second "plot" to get different scales on x2, y2 axes
        # self._axis2.tick_params(axis='both',direction='inout',
        #     length=8,bottom=False,left=False,top=True,right=True,pad=2)
        # self._axis2.yaxis.set_label_position('right')
        # self._axis2.yaxis.set_ticks_position('right')
        sa = self._s.spectral_axis
        stop = self.spectrogram.shape[1]
        step = self.spectrogram.shape[1] / self.spectrogram.shape[0]
        # print(stop,step,np.arange(0,stop,step).shape)
        im2 = self._axis2.plot(np.arange(0, stop, step), sa, linewidth=0)  # noqa: F841

        self._axis.set_xlim(0, stop - 0.5)
        z_label = self._set_labels(self._s)
        self._figure.colorbar(self.im, label=z_label, pad=0.1)

    def _set_labels(self, s):
        # x1: bottom
        # x2: top
        # y1: left
        # y2: right
        # z: colorbar

        x1_label = "Integration"
        self._axis.set_xlabel(x1_label)

        y1_label = "Channel"
        self._axis.set_ylabel(y1_label)

        y2_unit = s.spectral_axis.unit
        if y2_unit.is_equivalent(u.Hz):
            nu = r"$\nu$"
            y2_label = f"{nu} ({y2_unit})"
        self._axis2.set_ylabel(y2_label)

        z_unit = s.unit
        if z_unit.is_equivalent(u.K):
            z_label = f"$T_A$ ({z_unit})"
        elif z_unit.is_equivalent(u.Jy):
            snu = r"$S_{\nu}$"
            z_label = f"{snu} ({z_unit})"
        return z_label

    def reset(self):
        self._plot_kwargs = {
            "title": None,
            "cmap": "inferno",
            "interpolation": "nearest",
        }

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
            f"Int : {time_formatter(s.meta['EXPOSURE'])}", (hcoords[1], vcoords[1]), xycoords=xyc, size=fsize_small
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
        _tsys = np.around(s.meta["TSYS"], 2)
        self._axis.annotate(f"Tsys   :  {_tsys}", (hcoords[4], vcoords[0]), xycoords=xyc, size=fsize_small)
        self._axis.annotate(
            f"Tcal   :  {np.around(s.meta['TCAL'], 2)}", (hcoords[4], vcoords[1]), xycoords=xyc, size=fsize_small
        )
        self._axis.annotate(f"{s.meta['PROC']}", (hcoords[4], vcoords[2]), xycoords=xyc, size=fsize_small)

        # bottom row
        vcoord_bot = 0.82
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
