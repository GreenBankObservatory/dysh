"""
VaneSpectrum
"""

import numpy as np
from astropy import units as u
from astropy.time import Time

from dysh.log import logger

from ..util.docstring_manip import docstring_parameter
from ..util.gaincorrection import GBTGainCorrection
from ..util.weatherforecast import GBTWeatherForecast
from .core import mean_data
from .spectrum import Spectrum

tcal_param_docs = """tcal : None or float
    The calibration temperature in K.
    If `None`, it will be estimated using :meth:`get_tcal`.
"""

twarm_param_docs = """twarm : None or float
    The vane temperature in K.
    If `None`, it will be inferred from the TWARM column.
"""

zenith_opacity_param_docs = """zenith_opacity : float or None
    The zenith opacity in nepers.
    If `None`, it will be retrieved from the GBO weather forecast scripts (only available at GBO).
    If `None` and not at GBO, `tcal` will equal `twarm` (accurate to approximately 10%).
"""

tatm_param_docs = """tatm : float or None
    The atmospheric temperature towards the zenith in K.
    If `None`, it will be retrieved from the GBO weather forecast scripts (only available at GBO).
    If `None` and not at GBO, `tcal` will equal `twarm` (accurate to approximately 10%).
"""

tbkg_param_docs = """tbkg : float or None
    The background temperature in K.
    Default is the CMB temperature at 3 mm (2.725 K).
"""

param_docs = {
    "tcal": tcal_param_docs,
    "twarm": twarm_param_docs,
    "zenith_opacity": zenith_opacity_param_docs,
    "tatm": tatm_param_docs,
    "tbkg": tbkg_param_docs,
}

vane_kwargs = "".join(param_docs.values())


@docstring_parameter(vane_kwargs)
class VaneSpectrum(Spectrum):
    r"""
    Vane calibration spectrum.

    Parameters
    ----------
    vane : `~numpy.ndarray`
        Vane power values.
        The values will define the `data` attribute of the `~dysh.spectra.vane.VaneSpectrum` object.
    scan : int
        Scan number.
    fdnum : int
        The feed number.
    ifnum : int
        The intermediate frequency (IF) number.
    plnum : int
        The polarization number.
    {0}
    """

    def __init__(
        self,
        vane,
        scan,
        fdnum,
        ifnum,
        plnum,
        tcal=None,
        twarm=None,
        zenith_opacity=None,
        tatm=None,
        tbkg=None,
        *args,
        **kwargs,
    ):
        Spectrum.__init__(self, *args, **kwargs)
        self._vane = vane
        self._scan = scan
        self._fdnum = fdnum
        self._ifnum = ifnum
        self._plnum = plnum
        self._tcal = tcal
        self._twarm = twarm
        self._zenith_opacity = zenith_opacity
        self._tatm = tatm
        self._gbtgc = GBTGainCorrection()
        try:
            self._gbtwf = GBTWeatherForecast()
        except ValueError:
            logger.info("Weather forecast not available.")
            self._gbtwf = None
        self._validate_inputs()

    def _validate_inputs(self):
        """
        Check that the minimum information is available to compute
        the equivalent temperature of the vane.
        """
        if self._gbtwf is None:
            if self._tcal is None:
                if self.twarm is None:
                    raise ValueError("twarm not known. Set with twarm argument.")

    @property
    def scan(self):
        """The scan number."""
        return self._scan

    @property
    def ifnum(self):
        """The intermediate frequency (IF) number."""
        return self._ifnum

    @property
    def fdnum(self):
        """The feed number."""
        return self._fdnum

    @property
    def plnum(self):
        """The polarization number."""
        return self._plnum

    @classmethod
    @docstring_parameter(vane_kwargs)
    def from_spectrum(
        cls,
        spectrum,
        scan,
        fdnum,
        ifnum,
        plnum,
        vane=None,
        tcal=None,
        twarm=None,
        zenith_opacity=None,
        tatm=None,
        tbkg=None,
    ):
        """
        Returns a `~dysh.spectra.vane.VaneSpectrum` object
        given a `~dysh.spectra.spectrum.Spectrum` object.

        Parameters
        ----------
        spectrum : `~dysh.spectra.spectrum.Spectrum`
            Spectrum object. Its attributes will be copied into the
            `~dysh.spectra.vane.VaneSpectrum` object, except for the `data` and `flux` attributes
            which can be defined by the `vane` parameter.
        scan : int
            Scan number.
        fdnum : int
            The feed number.
        ifnum : int
            The intermediate frequency (IF) number.
        plnum : int
            The polarization number.
        vane : None or `~numpy.ndarray`
            Vane power values. If None, then it will use the `flux` attribute of `spectrum`.
            If set, the values will define the `data` attribute of the `~dysh.spectra.vane.VaneSpectrum` object.
        {0}

        Returns
        -------
        `~dysh.spectra.vane.VaneSpectrum`
            A `~dysh.spectra.vane.VaneSpectrum` object.
        """
        if vane is None:
            vane = spectrum.flux

        return cls(
            vane,
            scan,
            ifnum,
            plnum,
            fdnum,
            twarm=twarm,
            zenith_opacity=zenith_opacity,
            tatm=tatm,
            tbkg=tbkg,
            tcal=tcal,
            flux=vane,
            wcs=spectrum.wcs,
            meta=spectrum.meta,
            velocity_convention=spectrum.velocity_convention,
            radial_velocity=spectrum.radial_velocity,
            rest_value=spectrum.rest_value,
            observer=spectrum.observer,
            target=spectrum.target,
            mask=spectrum.mask,
        )

    @property
    def twarm(self):
        """
        Vane temperature in K.
        """
        # If twarm has not been explicitly set, then use the value in the header.
        if self._twarm is not None:
            return self._twarm
        else:
            return self.meta["TWARM"] + 273.15

    @docstring_parameter(
        "".join(f"{v}" for k, v in param_docs.items() if k in ["twarm", "tatm", "tbkg", "zenith_opacity"])
    )
    def get_tcal(
        self,
        ref,
        mjd: float | None = None,
        freq: u.Quantity | None = None,
        elev: u.Quantity | None = None,
        zenith_opacity: float | None = None,
        tatm: float | None = None,
        twarm: float | None = None,
        tbkg: float = 2.725,
    ):
        r"""
        Calibration temperature.

        Parameters
        ----------
        ref : `~dysh.spectra.spectrum.Spectrum`
            Reference spectrum used to derive :math:`T_{{\mathrm{{cal}}}}`.
        mjd : float or None
            Modified Julian date. If None, will use DATE-OBS in `ref.meta`.
        freq : `~astropy.units.Quantity` or None
            Frequency at which to compute the calibration temperature. If None, will use OBSFREQ in `ref.meta`.
        elev : `~astropy.units.Quantity` or None
            Elevation at which to evaluate the airmass. If None, will use ELEVATIO in `ref.meta`.
        {0}

        Returns
        -------
        tcal : float
            Calibration temperature for the vane in K.
        """

        # If the user set this value upstream. Use it.
        if self._tcal is not None:
            return self._tcal

        if mjd is None:
            mjd = Time(ref.meta["DATE-OBS"]).mjd
        if freq is None:
            freq = ref.meta["OBSFREQ"] * u.Hz
        if elev is None:
            elev = ref.meta["ELEVATIO"] * u.deg
        return self._get_tcal(freq, mjd, elev, zenith_opacity, tatm, twarm, tbkg)

    def _get_tcal(
        self,
        freq: u.Quantity,
        mjd: float,
        elevation: u.Quantity,
        zenith_opacity: float | None = None,
        tatm: float | None = None,
        twarm: float | None = None,
        tbkg: float = 2.725,
    ):
        # If the user set this value upstream. Use it.
        if self._tcal is not None:
            return self._tcal

        if twarm is None:
            twarm = self.twarm
            logger.info(f"Vane temperature (twarm): {twarm:.2f} K")

        if zenith_opacity is None:
            zenith_opacity = self._zenith_opacity
        if zenith_opacity is None:
            if self._gbtwf is not None:
                result = self._gbtwf.fetch(
                    vartype="Opacity",
                    specval=freq,
                    mjd=mjd,
                )
                zenith_opacity = result[:, -1]
                logger.info(f"Mean zenith opacity: {np.mean(zenith_opacity):.2f} nepers")
            else:
                logger.debug("Could not get forecasted zenith opacity ")
        logger.debug(f"Zenith opacity: {zenith_opacity} nepers")

        if tatm is None:
            tatm = self._tatm
        if tatm is None:
            if self._gbtwf is not None:
                result = self._gbtwf.fetch(vartype="Tatm", specval=freq, mjd=mjd)
                tatm = result[:, -1]
                logger.info(f"Mean atmospheric temperature: {np.mean(tatm):.2f} K")
            else:
                logger.debug("Could not get forecasted atmospheric temperature ")
        logger.debug(f"Atmospheric temperature: {tatm} K")

        if tatm is not None and zenith_opacity is not None:
            airmass = self._gbtgc.airmass(elevation, zd=False)
            tcal = (tatm - tbkg) + (twarm - tatm) * np.exp(zenith_opacity * airmass)
        else:
            tcal = twarm
            missing = []
            if zenith_opacity is None:
                missing.extend(["zenith opacity"])
            if tatm is None:
                missing.extend(["atmospheric temperature"])
            missing = " nor ".join(missing)
            logger.info(
                f"No {missing} available. Will approximate the calibration temperature to the vane temperature {tcal:.2f} K"
            )

        logger.info(f"Mean calibration temperature (tcal): {np.mean(tcal):.2f} K")
        logger.debug(f"Calibration temperature (tcal): {tcal} K")

        return tcal

    def _get_tsys(self, ref, tcal):
        """ """
        mean_off = mean_data(ref)
        mean_dif = mean_data(self.data - ref)
        tsys = tcal * mean_off / mean_dif
        return tsys

    @docstring_parameter(param_docs["tcal"])
    def get_tsys(self, ref, tcal=None):
        """
        Compute the system temperature.

        Parameters
        ----------
        ref : `~dysh.spectra.spectrum.Spectrum`
            The reference spectrum.
        {0}

        Returns
        -------
        tsys : float
            The system temperature in K.
        """

        if tcal is None:
            tcal = self.get_tcal(ref)

        return self._get_tsys(ref.data, tcal)
