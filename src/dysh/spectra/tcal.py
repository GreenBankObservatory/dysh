"""
TCal class
"""

import numpy as np

from ..util.docstring_manip import insert_docstr_section
from .spectrum import Spectrum


class TCal(Spectrum):
    """
    Noise diode temperature.
    This behaves as a `~dysh.spectra.spectrum.Spectrum` with the addition of
    the `name` and `snu` attributes, which hold the name of the source used to derive
    the noise diode temperature and its flux density.
    """

    def __init__(self, name, snu, *args, **kwargs):
        Spectrum.__init__(self, *args, **kwargs)
        self._name = name
        self._snu = snu

    @property
    def name(self):
        """Calibration source name."""
        return self._name

    @property
    def snu(self):
        """Calibration source flux density."""
        return self._snu

    @classmethod
    def from_spectrum(cls, spectrum, name, snu, data=None):
        """
        Returns a `~dysh.spectra.tcal.TCal` object
        given a `~dysh.spectra.spectrum.Spectrum` object.

        Parameters
        ----------
        spectrum : `~dysh.spectra.spectrum.Spectrum`
            Spectrum object. Its attributes will be copied into the
            `~dysh.spectra.tcal.TCal` object, except for the `data` and `flux` attributes
            which can be defined by the `data` parameter.
        name : str
            Name of the calibrator used to derive the noise diode temperature.
        snu : `~numpy.ndarray`
            The flux values of the calibrator source.
        data : None or `~numpy.ndarray`
            Noise diode temperature values. If None, then it will use the `flux` attribute of `spectrum`.
            If set, the values will define the `flux` attribute of the `~dysh.spectra.tcal.TCal` object.

        Returns
        -------
        `~dysh.spectra.tcal.TCal`
            The `~dysh.spectra.tcal.TCal` object.
        """
        if data is None:
            data = spectrum.flux

        return cls(
            name,
            snu,
            flux=data,
            wcs=spectrum.wcs,
            meta=spectrum.meta,
            velocity_convention=spectrum.velocity_convention,
            radial_velocity=spectrum.radial_velocity,
            rest_value=spectrum.rest_value,
            observer=spectrum.observer,
            target=spectrum.target,
            mask=spectrum.mask,
        )

    @insert_docstr_section(Spectrum.smooth.__doc__, section="Parameters")
    def smooth(self, *args, **kwargs):
        """
        Smooth or convolve `TCal`, optionally decimating it.

        Parameters
        ----------
        {0}

        Returns
        -------
        `~dysh.spectra.tcal.TCal`
            Smoothed noise diode temperature.
        """
        kwargs.setdefault("decimate", -1)
        tcal_smo = super().smooth(*args, **kwargs)
        return TCal.from_spectrum(tcal_smo, self.name, self.snu)

    def get_tcal(self, frac=0.1, method=np.nanmean, dtype=np.float32):
        """
        Reduce the temperature of the noise diode to a single value.

        Parameters
        ----------
        frac : float
            The fraction of the channels to discard on both ends of the spectrum.
            Default: 0.1
        method : callable
            The function used to reduce the values.
            Default : `numpy.nanmean`
        dtype : str or dtype
            Typecode or data-type to which the output value, `tcal` is cast.

        Returns
        -------
        tcal : float
            The noise diode temperature reduced using `method` over the channel range defined by `frac`.
        """

        ch0 = int(self.nchan * frac)
        chf = int(self.nchan * (1 - frac))
        s = slice(ch0, chf)
        self.tcal = method(self.data[s])

        return self.tcal.astype(dtype)

    @property
    def nchan(self):
        return len(self.data)
