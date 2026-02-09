"""
Plot a spectrum using matplotlib
"""

import warnings
from copy import deepcopy

import astropy.units as u
import matplotlib as mpl
import numpy as np
from astropy.io import fits
from matplotlib.text import OffsetFrom
from matplotlib.ticker import AutoLocator, MaxNLocator

from dysh.log import logger

from .plotbase import PlotBase

_KMS = u.km / u.s


class ScanPlot(PlotBase):
    r"""
    The ScanPlot class is for simple plotting of a `~dysh.spectra.scan.Scan` or `~dysh.spectra.scan.ScanBlock`
    using matplotlib functions. Plot attributes are modified using keywords
    (\*\*kwargs) described below. SpectrumPlot will attempt to make smart default
    choices for the plot if no additional keywords are given.

    Parameters
    ----------
    scanblock_or_scan : `~dysh.spectra.scan.Scan` or `~dysh.spectra.scan.ScanBlock`
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
        self.reset()
        self._scanblock_or_scan = scanblock_or_scan
        self._plot_kwargs.update(kwargs)
        self._axis2 = None
        acceptable_types = ["PSScan", "TPScan", "NodScan", "FSScan", "SubBeamNodScan"]
        self._spectrum = self._scanblock_or_scan.timeaverage()
        self._sa = self._spectrum.spectral_axis

        # determine if input is a ScanBlock or a ScanBase (raise exception if neither)
        self._type = scanblock_or_scan.__class__.__name__
        if self._type == "ScanBlock":
            self._scanblock = scanblock_or_scan
            self._num_scans = len(self._scanblock)
            self._scan_numbers = np.empty(self._num_scans, dtype=int)
        elif self._type in acceptable_types:
            self._scan = scanblock_or_scan
            self._scan_numbers = np.empty(1, dtype=int)
        else:
            raise Exception(f"Plotter input {self._type} does not appear to be a valid input object type")

        # handle scanblocks
        if self._type == "ScanBlock":
            self._nint_nos = []  # number of integrations in each scan
            self._timestamps = []  # 0-indexed timestamps in sec for every integration
            xtick_labels = []  # intnum labels for multiple-scan scanblocks
            for i, scan in enumerate(self._scanblock):
                if i == 0:
                    self.spectrogram = scan._calibrated
                else:
                    self.spectrogram = np.append(self.spectrogram, scan._calibrated, axis=0)
                self._scan_numbers[i] = scan.scan
                self._nint_nos.append(scan.nint)  # not sure if I need this
                xtick_labels.append(np.r_[0 : scan.nint])
                # TODO: figure out how to deal with generating a "time" axis
                # agnostic of scan proctype (pos sw, etc will have gaps between scans due to OFF)
                # self._timestamps.append(scan.)
            xtick_labels = np.concatenate(xtick_labels, axis=0)

        # handle scans
        elif self._type in acceptable_types:
            self._scan_numbers[0] = self._scan.scan
            self.spectrogram = self._scan._calibrated
            self._nint_nos = [self._scan.nint]
            xtick_labels = np.r_[0 : self._scan.nint]

        self._xtick_labels = xtick_labels
        self.spectrogram = self.spectrogram.T

    def reset(self):
        """Reset the plot keyword arguments to their defaults."""
        self._plot_kwargs = {
            "title": None,
            "cmap": "inferno",
            "interpolation": "nearest",
        }

    def write(self, filename, avechan=1, chan=None, overwrite=False):
        r"""
        Write the current spectrogram as a FITS image. No WCS is maintained yet.
        The integrations are in the first axis and the channels in the second.

        Parameters
        ----------
        filename : str
           Filename of FITS file to be saved. No default
        avechan : int
           Averaging number of channels. Only supports averaging by a multiple
           of the number of channels. Default: 1
        chan : [int,int]
           If given, it will select this channel range for output.
           Default: all channels will be selected
        overwrite : boolean
           Overwrite existing file. Default: False
        """
        # skip mask for now
        # add the dysh history?   not available here
        # @todo future multi-beams like argus could do 3D fits cube
        data = self.spectrogram.data
        if len(data.shape) != 2:
            logger.error(f"spectrogram is {len(data.shape)}D, can only write 2D")
            return
        if chan is not None:
            data = data[chan[0] : chan[1], :]
        if avechan == 1:
            hdu = fits.PrimaryHDU(data)
        else:
            (ny, nx) = data.shape
            if ny % avechan != 0:
                logger.info(f"{ny} not a multiple of {avechan}")
                return
            data = data.T.reshape((nx, -1, avechan)).mean(axis=2).T
            hdu = fits.PrimaryHDU(data)
        hdu.writeto(filename, overwrite=overwrite)
        logger.info(f"Wrote {filename} with size {data.shape[1]} x {data.shape[0]}")

    def plot(self, show_header=True, spectral_unit=None, **kwargs):
        r"""
        Plot the scan.

        Parameters
        ----------
        spectral_unit : `~astropy.units.Unit`
            The units to use on the frequency axis. Default: MHz if below 1 GHz, GHz if above. Otherwise, can be any valid frequency unit.
        **kwargs : various
            keyword=value arguments drawn from `~matplotlib.axes.Axes.imshow` kwargs.
            Currently implemented kwargs include `cmap`, `interpolation`, `vmin`, `vmax`, and `norm`.
        """

        this_plot_kwargs = deepcopy(self._plot_kwargs)
        this_plot_kwargs.update(kwargs)

        cmap = kwargs.get("cmap", "inferno")
        interpolation = kwargs.get("interpolation", "nearest")
        vmin = kwargs.get("vmin", None)
        vmax = kwargs.get("vmax", None)
        norm = kwargs.get("norm", None)

        if True:
            self.figure = mpl.figure.Figure(figsize=(10, 6), dpi=100)
            self.axes = self.figure.subplots(nrows=1, ncols=1)
            self._axis2 = self.axes.twinx()
            self._axis3 = self.axes.twiny()

        if show_header:
            # This has to be set before drawing the image and adding the colorbar.
            self.figure.subplots_adjust(top=0.79, left=0.1, right=1.05)
            self._set_header(self._spectrum)

        self.im = self.axes.imshow(
            self.spectrogram, aspect="auto", cmap=cmap, interpolation=interpolation, vmin=vmin, vmax=vmax, norm=norm
        )

        # address intnum labelling for len(scanblock) > 1
        self.axes.set_xticks(np.arange(self.spectrogram.shape[1]), self._xtick_labels)
        if len(self._xtick_labels) == 1:
            locator = MaxNLocator(nbins=len(self._xtick_labels), integer=True, min_n_ticks=1)
        else:
            locator = AutoLocator()
        self.axes.xaxis.set_major_locator(locator)

        # second "plot" to get different scales on x2, y2 axes
        if spectral_unit is not None:
            self._sa = self._sa.to(spectral_unit)
        elif self._sa[0] / (u.GHz) < 1:
            self._sa = self._sa.to(u.MHz)
        else:
            self._sa = self._sa.to(u.GHz)
        stop = self.spectrogram.shape[1]
        step = self.spectrogram.shape[1] / self.spectrogram.shape[0]
        im2 = self._axis2.plot(np.arange(0, stop, step), self._sa, linewidth=0)  # noqa: F841
        self._axis2.set_ylim((np.min(self._sa).value, np.max(self._sa).value))

        # third axis to plot the scan numbers
        im3 = self._axis3.plot(np.arange(0, stop, step), self._sa, linewidth=0)  # noqa: F841
        # determine tick locations and labels
        tick_locs = []
        acc = 0
        for numints in self._nint_nos:
            tick_locs.append(acc)
            acc += numints
        self._axis3.set_xticks(tick_locs)
        self._axis3.set_xticklabels(self._scan_numbers)
        fsize = 15
        x1_alt_padding = mpl.rcParams["axes.labelpad"] + fsize
        self._axis3.tick_params(
            axis="x",
            width=0,
            pad=x1_alt_padding + 9,
            # labelsize=fsize,
            bottom=True,
            top=False,
            labelbottom=True,
            labeltop=False,
        )

        self.axes.set_xlim(0, stop - 0.5)
        self._set_labels()

    def _set_labels(self):
        # x1: bottom
        # x2: top
        # y1: left
        # y2: right
        # z: colorbar

        x1_label = "Integration"
        self.axes.set_xlabel(x1_label)

        x1_alt_label = "Scan"
        off = OffsetFrom(self._axis3.get_xticklabels()[0], (0.0, 0.0))
        self._axis3.annotate(x1_alt_label, xy=(0.5, 0.5), xytext=(-10, 0.0), textcoords=off, va="bottom", ha="right")

        y1_label = "Channel"
        self.axes.set_ylabel(y1_label)

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
        self._colorbar = self.figure.colorbar(self.im, label=z_label, pad=0.1)
        # matplotlib won't set this before the Figure is drawn.
        self.figure.draw_without_rendering()
        # If there's an offset, add it to the label and make the offset invisible.
        if self._colorbar.ax.yaxis.offsetText.get_text() != "":
            off = self._colorbar.ax.yaxis.offsetText.get_text()
            e = off.split("e")[1]
            self._colorbar.set_label(z_label + rf"($\times10^{{{e}}}$)")
            self._colorbar.ax.yaxis.offsetText.set_visible(False)

    def set_clim(self, vmin=None, vmax=None):
        """
        Set the vmin and vmax parameters of the image.

        Parameters
        ----------
        vmin : float
            The minimum value of the color scale. Default None; to autoscale.
        vmax : float
            The maximum value of the color scale. Default None; to autoscale.
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

    def set_norm(self, norm=None):
        """
        Set the normalization of the image colormap.
        Can be any value supported by the `~matplotlib.axes.Axes.imshow` `norm` keyword,
        e.g. 'linear', 'log' etc or a Normalize object.

        Parameters
        ----------
        norm : str | Normalize | None
            norm used for the color scale. Default: "None", for linear scaling.
        """
        self.im.set_norm(norm)
