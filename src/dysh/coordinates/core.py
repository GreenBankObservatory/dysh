"""
Core functions/classes for spatial and velocity coordinates and reference frames
"""

import astropy.coordinates as coord
import astropy.units as u
import numpy as np
from astropy.coordinates.spectral_coordinate import (
    DEFAULT_DISTANCE as _DEFAULT_DISTANCE,
)
from astropy.time import Time

from dysh.log import logger

_PMZERO = 0.0 * u.mas / u.yr
_PMZERORAD = 0.0 * u.rad / u.s
_VELZERO = 0.0 * u.km / u.s
_MPS = u.m / u.s
KMS = u.km / u.s

# Velocity frame conventions, partially stolen from pyspeckit.
# See also Section6.2.5.2 of the GBT observer's guide https://www.gb.nrao.edu/scienceDocs/GBTog.pdf

# string to astropy coordinate frame class
astropy_frame_dict = {
    "VLSR": coord.LSRK,
    "VRAD": coord.LSRK,
    "VELO": coord.LSRK,
    "VOPT": coord.LSRK,
    "LSRD": coord.LSRD,
    "lsrd": coord.LSRD,
    "LSRK": coord.LSRK,
    "lsrk": coord.LSRK,
    "-LSR": coord.LSRK,
    "-HEL": coord.HCRS,
    "-BAR": coord.ICRS,
    "BAR": coord.ICRS,
    "BARY": coord.ICRS,
    "icrs": coord.ICRS,
    "ICRS": coord.ICRS,
    "bary": coord.ICRS,
    "barycentric": coord.ICRS,
    "VHEL": coord.HCRS,
    "heliocentric": coord.HCRS,
    "helio": coord.HCRS,
    "hcrs": coord.HCRS,
    "HCRS": coord.HCRS,
    "VGEO": coord.GCRS,
    "geocentric": coord.GCRS,
    "gcrs": coord.GCRS,
    "GCRS": coord.GCRS,
    "-GAL": coord.Galactocentric,
    "galactocentric": coord.Galactocentric,
    "topocentric": coord.ITRS,  # but need to add observatory position
    "topo": coord.ITRS,  # but need to add observatory position
    "itrs": coord.ITRS,  # but need to add observatory position
    "fk5": coord.FK5,
    "gal": coord.Galactic,
}

astropy_convenience_frame_names = {
    "bary": "icrs",
    "barycentric": "icrs",
    "heliocentric": "hcrs",
    "helio": "hcrs",
    "geo": "gcrs",
    "geocentric": "gcrs",
    "topocentric": "itrs",
    "topo": "itrs",
    "vlsr": "lsrk",
}


# astropy-sanctioned coordinate frame string to label
frame_to_label = {
    "itrs": "Topocentric",
    "topocentric": "Topocentric",
    "topo": "Topocentric",
    "hcrs": "Heliocentric",
    "gcrs": "Geocentric",
    "icrs": "Barycentric",
    "fk5": "Barycentric",
    "lsrk": "LSRK",
    "lsrd": "LSRD",
    "galactocentric": "Galactocentric",
}

# velframe string to label
frame_dict = {
    "VLSR": "LSRK",
    "VRAD": "LSRK",
    "VELO": "LSRK",
    "VOPT": "LSRK",
    "LSRD": "LSRD",
    "LSRK": "LSRK",
    "-LSR": "LSRK",
    "-HEL": "heliocentric",
    "-BAR": "barycentric",
    "BAR": "barycentric",
    "BARY": "barycentric",
    "-OBS": "topocentric",
    "VHEL": "heliocentric",
    "VGEO": "geocentric",
    "TOPO": "topocentric",
    "TRUE": "topocentric",
    "VREST": "rest",
    "Z": "rest",
    "FREQ": "rest",
    "WAV": "rest",
    "WAVE": "rest",
    "CMB": "cmb",
    "GALAC": "galactic",
    "GALA": "galactic",
    "ALAC": "galactic",  # in case of 'VELGALAC', last 4 chars.
}


# label to velframe string
reverse_frame_dict = {
    "bary": "-BAR",
    "barycentric": "-BAR",
    "icrs": "-BAR",
    "galactocentric": "-GAL",
    "gcrs": "-GEO",
    "geocentric": "-GEO",
    "heliocentric": "-HEL",
    "hcrs": "-HEL",
    "helio": "-HEL",
    "heli": "-HEL",
    "lsr": "-LSR",
    "lsrk": "-LSR",
    "lsrd": "LSRD",
    "itrs": "-OBS",
    "topo": "-OBS",
    "topocentric": "-OBS",
    "fk5": "-BAR",
    "fk4": "-BAR",
    "cirs": "-GEO",
    "tete": "-GEO",
}
# Dictionary to convert from FITS velocity convention to specutils string.
# At GBT, VELO was written by sdfits filler for some unknown amount of
# time instead of RELA, so allow for it here
velocity_convention_dict = {
    "OPTI": "optical",
    "RADI": "radio",
    "RELA": "relativistic",
    "VELO": "relativistic",
}

# reverse of above
reverse_velocity_convention_dict = {"optical": "OPTI", "radio": "RADI", "relativistic": "VELO"}

# dictionary to convert CRVAL4 to a polarization ID for plotting purposes
# "zero" entry for more graceful handling downstream
crval4_to_pol = {
    -1: "RR",
    -2: "LL",
    -3: "RL",
    -4: "LR",
    -5: "XX",
    -6: "YY",
    -7: "YX",
    -8: "XY",
    0: "UNKNOWN",
    1: "I",
    2: "Q",
    3: "U",
    4: "V",
}


def replace_convention(veldef, doppler_convention):
    return reverse_velocity_convention_dict[doppler_convention] + veldef[4:]


# This gives the wrong answer for GBT which always writes data as topocentric
# regardless of what VELDEF says.  Therefore the kludge for GBT: Instead of
# VELDEF, pass in CTYPE1  which is always 'FREQ-OBS'.
def is_topocentric(veldef):
    # fmt: off
    if veldef[0:4] == "FREQ":  # Workaround for GBT peculiarity described above.
                               # We don't expect other observatories will
                               # use GBT's convention.
        # fmt: on
        nvd = "VELO" + veldef[4:]
        return decode_veldef(nvd)[1] == "topocentric"
    else:
        return decode_veldef(veldef)[1] == "topocentric"


def decode_veldef(veldef):
    """
    Parse the SDFITS VELDEF value into its two components, the velocity
    definition and velocity reference frame.  This value must contain
    no more than 8 characters where the first 4 characters describe the velocity
    definition and the last 4 characters describe the reference frame.

    Parameters
    ----------
    veldef : str
        The definition string, consisting of a velocity convention and a velocity frame, e.g.,  'OPTI-LSR'

    Returns
    -------
    A str tuple of velocity convention and velocity frame type, e.g., ('radio', 'LSRK')
    """
    if len(veldef) > 8:
        # in future, possibly relax this requirement
        # if string not coming from FITS
        raise ValueError(f"VELDEF string {veldef} must be no more than 8 characters.")
    vconv = veldef[:4]
    try:
        velocity_convention = velocity_convention_dict[vconv]
    except KeyError:
        raise KeyError(f"Velocity convention {vconv} not recognized.")  # noqa: B904

    frame = veldef[4:]
    try:
        frame_type = frame_dict[frame]
    except KeyError:
        raise KeyError(f"Velocity frame {frame} not recognized.")  # noqa: B904

    return velocity_convention, frame_type


def veldef_to_convention(veldef):
    """given a VELDEF, return the velocity convention expected by Spectrum(1D)

    Parameters
    ----------
    veldef : str
        Velocity definition from FITS header, e.g., 'OPTI-HELO', 'VELO-LSR'

    Returns
    -------
    convention : str
        Velocity convention string, one of {'radio', 'optical', 'relativistic'}  or None if `velframe` can't be parsed
    """

    prefix = veldef[0:4].lower()
    if prefix == "opti":
        return "optical"
    if prefix == "velo" or prefix == "radi":
        return "radio"
    if prefix == "rela":
        return "relativistic"
    return None


def sanitize_skycoord(target):
    """Method to enforce certain attributes of input SkyCoordinate in
    order to workaround astropy bug that distance and proper motions
    need to be explicitly set for certain coordinate conversions, even
    if they are zero.  See `explanation here
    <https://community.openastronomy.org/t/exception-raised-when-converting-from-lsrk-to-other-frames/841/2>_` and
    and `Astropy GitHub issue 12371 <https://github.com/astropy/astropy/issues/12731>`_.


    Parameters
    ----------
        target : `~astropy.coordinates.SkyCoordinate`
        The input Sky Coordinate

    Returns
    -------
        sanitized_target : `~astropy.coordinates.SkyCoordinate`
        Target with distance, radial velocity, and proper motion set.

    """
    if not isinstance(target, coord.SkyCoord):
        raise TypeError("Target must be instance of astropy.coordinates.SkyCoord")
    if hasattr(target, "sanitized"):  # don't do it twice.
        if target.sanitized:
            return target
    try:
        # This will actually raise an exception
        # if radial velocity wasn't specified rather than
        # just return False!
        # (but only for Galactic!)
        if hasattr(target, "radial_velocity"):
            _rv = target.radial_velocity
        else:
            _rv = _VELZERO
    except Exception:
        _rv = _VELZERO

    # print(f"{type(target.distance)}, [{target.distance.unit}], {target.distance.unit == u.dimensionless_unscaled}")
    if target.distance.unit == u.dimensionless_unscaled and round(target.distance.value) == 1:
        # distance was unset and astropy set it to 1 with a dimensionless composite unit
        newdistance = _DEFAULT_DISTANCE
    else:
        newdistance = target.distance

    if hasattr(target, "ra"):  # RADEC based
        lon = target.ra
        lat = target.dec
        pm_lon = target.pm_ra_cosdec
        pm_lat = target.pm_dec
        # could probably be clever with kwargs
        # and avoid doing this twice
        _target = coord.SkyCoord(
            lon,
            lat,
            frame=target.frame,
            distance=newdistance,
            pm_ra_cosdec=pm_lon,
            pm_dec=pm_lat,
            radial_velocity=_rv,
        )
    # ====== GALACTIC COORDS HAVE NOT BEEN FULLY TESTED. USE WITH CAUTION ====
    elif hasattr(target, "l"):  # Galactic
        lon = target.l
        lat = target.b
        try:
            # This will actually raise an exception
            # if proper motions weren't specified rather than
            # just return False!
            # (but only for Galactic!)
            hasattr(target, "pm_l_cosb")
            # WTF. If distance was give but no PM, then units come out at 'km rad /s'
            if "km" in str(target.pm_l_cosb):
                pm_lon = _PMZERORAD
                pm_lat = _PMZERORAD
            else:
                pm_lon = target.pm_l_cosb
                pm_lat = target.pm_b
        except:  # noqa: E722
            pm_lon = _PMZERORAD
            pm_lat = _PMZERORAD
        # print(
        #    f"DEBUG\n: _target = SkyCoord( {lon}, {lat}, frame={target.frame}, distance={newdistance},"
        #    f" pm_l_cosb={pm_lon}, pm_b={pm_lat}, radial_velocity={_rv})"
        # )
        _target = coord.SkyCoord(
            lon,
            lat,
            frame=target.frame,
            distance=newdistance,
            pm_l_cosb=pm_lon,
            pm_b=pm_lat,
            radial_velocity=_rv,
        )
    elif hasattr(target, "az"):  # AzEl or AltAz
        lon = target.az
        lat = target.alt
        pm_lon = target.pm_az_cosalt
        pm_lat = target.pm_alt
        _target = coord.SkyCoord(
            lon,
            lat,
            frame=target.frame,
            distance=newdistance,
            pm_az_cosalt=pm_lon,
            pm_alt=pm_lat,
            radial_velocity=_rv,
        )
    elif hasattr(target, "ha"):  # HADec
        lon = target.ha
        lat = target.dec
        pm_lon = target.pm_ha_cosdec
        pm_lat = target.pm_dec
        _target = coord.SkyCoord(
            lon,
            lat,
            frame=target.frame,
            distance=newdistance,
            pm_ha_cosdec=pm_lon,
            pm_dec=pm_lat,
            radial_velocity=_rv,
        )
    else:
        logger.warning(f"Can't sanitize {target}")
        return target

    _target.sanitized = True
    return _target


# @todo version that takes a SpectralCoord
def topocentric_velocity_to_frame(target, toframe, observer, obstime):
    """Compute the difference in topocentric velocity and the velocity in the input frame.

    Parameters
    ----------
    target: `~astropy.coordinates.SkyCoord`
        The sky coordinates of the object including proper motion and distance. Must be in ICRS
    target: `~astropy.coordinates.SkyCoord`
        The sky coordinates of the object including proper motion and distance. Must be in ICRS

    toframe: str
        The frame into which `coord` should be transformed, e.g.,  'icrs', 'lsrk', 'hcrs'.
        The string 'topo' is interpreted as 'itrs'.
        See astropy-supported reference frames (link)

    observer: `~astropy.coordinates.EarthLocation`
        The location of the observer

    obstime: `~astropy.time.Time`
        The time of the observation

    Returns
    -------
    radial_velocity : `~astropy.units.quantity.Quantity`
        The radial velocity of the source in `toframe`

    """
    if not isinstance(target.frame, coord.ICRS):
        _target = sanitize_skycoord(target.icrs)
    else:
        _target = sanitize_skycoord(target)
    # raise Exception("input frame must be ICRS")
    topocoord = observer.get_itrs(obstime=obstime)
    sc = coord.SpectralCoord(1 * u.Hz, observer=topocoord, target=_target)
    sc2 = sc.with_observer_stationary_relative_to(toframe)  # noqa: F841
    return sc.with_observer_stationary_relative_to(toframe).radial_velocity


def get_velocity_in_frame(target, toframe, observer=None, obstime=None):
    """Compute the radial velocity of a source in a new velocity frame.

    Parameters
    ----------
        target: `~astropy.coordinates.SkyCoord`
            The sky coordinates of the object including proper motion and distance.
            Note: In order to get around a bug in astropy (link), if the `target` frame or `toframe` is 'lsrk' (`~astropy.coordinates.LSRK`),

            done:

            * If proper motions attributes of `target` are not set, they will be set to zero.
            * Similarly, if distance attribute of `target` is not set, it will be set to a very large number.
            * This is done on a copy of the coordinate so as not to change the input object.

        toframe: str
            The frame into which `coord` should be transformed, e.g.,  'icrs', 'lsrk', 'hcrs'.
            The string 'topo' is interpreted as 'itrs'.
            See `astropy-supported reference frames
            <https://docs.astropy.org/en/stable/coordinates/index.html#module-astropy.coordinates.builtin_frames>`_.

        observer: `~astropy.coordinates.EarthLocation`
            The location of the observer required for certain transformations (e.g. to/from GCRS or ITRS)

        obstime: `~astropy.time.Time`
            The time of the observation, required for for certain transformations (e.g. to/from GCRS or ITRS)

    Returns
    -------
        radial_velocity : `~astropy.units.quantity.Quantity`
            The radial velocity of the source in `toframe`

    """
    _target = sanitize_skycoord(target)

    if toframe == "topo":
        vshift = -topocentric_velocity_to_frame(_target, _target.frame, observer, obstime)
        return vshift
    elif isinstance(_target.frame, coord.ITRS):
        vshift = topocentric_velocity_to_frame(_target, _target.frame, observer, obstime)
        return vshift
    # really gcrs should raise an exception because astropy does it wrong.
    if toframe.lower() == "gcrs" and obstime is not None:
        _toframe = coord.GCRS(obstime=obstime)
    else:
        _toframe = toframe
    return _target.transform_to(_toframe).radial_velocity


# @todo have a SpectralCoord version of this?
def veltofreq(velocity, restfreq, veldef):
    """Convert velocity to frequency using the given rest frequency and velocity definition.

    Parameters
    ----------

    velocity: `~astropy.units.quantity.Quantity`
        The velocity values
    restfreq: `~astropy.units.quantity.Quantity`
        The rest frequency
    veldef : str
        Velocity definition from FITS header, e.g., 'OPTI-HELO', 'VELO-LSR'.

    Returns
    -------
    frequency: `~astropy.units.quantity.Quantity`
        The velocity values converted to frequency using `restfreq` and `veldef`.

    """

    vdef = veldef_to_convention(veldef)
    if vdef is None:
        raise ValueError(f"Unrecognized VELDEF: {veldef}")
    if vdef == "radio":
        radio = u.doppler_radio(restfreq)
        frequency = velocity.to(u.Hz, equivalencies=radio)
    elif vdef == "optical":
        optical = u.doppler_optical(restfreq)
        frequency = velocity.to(u.Hz, equivalencies=optical)
    elif vdef == "relativistic":
        relativistic = u.doppler_relativistic(restfreq)
        frequency = velocity.to(u.Hz, equivalencies=relativistic)
    else:
        # should never get here.
        raise ValueError(f"Unrecognized velocity convention: {vdef}")
    return frequency


def change_veldef(ctype, toframe):
    unknown = "-UNK"
    prefix = ctype[0:4]
    newpostfix = reverse_frame_dict.get(toframe, unknown)
    newctype = prefix + newpostfix

    if newpostfix == unknown:
        logger.warning(f"Could not determine VELDEF for given frame {toframe}. Setting to {newctype}.")
    return newctype


def make_target(header):
    """
    Create a SkyCoord object from a SDFITS header dictionary CRVAL2,
    CRVAL3 are assumed to be the latitude-like and longitude-like
    coordinates. VELOCITY is taken to be the radial velocity.  Coordinate
    frame is determined from RADESYS.

    Parameters
    ----------
    header : dict
        The SDFITS header keyword dictionary

    Returns
    -------
    target : `~astropy.coordinates.SkyCoord`
        A SkyCoord object based on the input coordinates and radial velocity
    """

    # should we also require DATE-OBS or MJD-OBS?
    _required = set(["CRVAL2", "CRVAL3", "CUNIT2", "CUNIT3", "VELOCITY", "EQUINOX", "RADESYS", "DATE-OBS"])

    # for k in _required:
    #    print(f"{k} {k in header}")
    if not _required <= header.keys():
        raise ValueError(f"Header is missing one or more required keywords: {_required}")

    t1 = coord.SkyCoord(
        header["CRVAL2"],
        header["CRVAL3"],
        unit=(header["CUNIT2"], header["CUNIT3"]),
        frame=header["RADESYS"].lower(),
        radial_velocity=header["VELOCITY"] * _MPS,
        distance=_DEFAULT_DISTANCE,  # need this or PMs get units m rad /s !
        obstime=Time(header["DATE-OBS"]),
        location=gbt_location(),
    )
    # print(f"{t1},{t1.pm_ra_cosdec},{t1.pm_dec},{t1.distance},{t1.radial_velocity}")
    target = sanitize_skycoord(t1)
    return target


def gb20m_location():
    """
    Create an astropy EarthLocation for the 20 Meter using the same established by GBO.

    Returns
    -------
    gb20m : `~astropy.coordinates.EarthLocation`
        astropy EarthLocation for the GBT.
    """
    gb20m_lat = 38.43685 * u.deg
    gb20m_lon = -79.82552 * u.deg
    gb20m_height = 835.0 * u.m
    gb20m = coord.EarthLocation.from_geodetic(lon=gb20m_lon, lat=gb20m_lat, height=gb20m_height)
    return gb20m


class GB20M:
    """Singleton Green Bank 20 Meter Telescope EarthLocation object, using the Green Bank Observatory's official coordinates for the site."""

    _instance = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = gb20m_location()
        return cls._instance


def gbt_location():
    """
    Create an astropy `~astropy.coordinates.EarthLocation` for the GBT using the same established by GBO.
    See page 3 of `The GBT Tracking Local Oscillator FITS Keyword Definitions <https://www.gb.nrao.edu/GBT/MC/doc/dataproc/gbtLOFits/gbtLOFits.pdf>`_

    * latitude    = 38d 25m 59.265s N

    * longitude   = 79d 50m 23.419s W

    * height      = 854.83 m

    Note these differ from astropy's "GBT" `~astropy.coordinates.EarthLocation` by several meters.

    Returns
    -------
    gbt : `~astropy.coordinates.EarthLocation`
        astropy EarthLocation for the GBT
    """
    gbt_lat = 38.4331291667 * u.deg
    gbt_lon = -79.839838611 * u.deg
    gbt_height = 854.83 * u.m
    gbt = coord.EarthLocation.from_geodetic(lon=gbt_lon, lat=gbt_lat, height=gbt_height)
    return gbt


class GBT:
    """Singleton Green Bank Telescope EarthLocation object, using the Green Bank Observatory's official coordinates for the site."""

    _instance = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = gbt_location()
        return cls._instance


class Observatory:
    """Class that returns (latitude, longitude, altitude) of known
    observatories, using :class:`astropy.coordinates.EarthLocation`.
    This can be used for instance in transforming velocities between
    different reference frames.

    Examples
    --------
    .. code-block::

        obs = Observatory()
        print(obs['GBT'])
        print(obs['ALMA'])

    Alternatively, you can treat Observatory like a dict:

    .. code-block::

        gbt = Observatory["GBT"]

    """

    def __init__(self):
        # might be confusing API to have everyting as obs[string]
        # and just GBT as an attribute. Leave this unadvertised for now
        # in case I remove it.
        self.GBT = GBT()
        self.GB20M = GB20M()

    def __getitem__(self, key):
        # For GBT we want to use the GBO official location
        if key == "GBT" or key == "NRAO_GBT":
            return self.GBT
        elif key == "GB20M":
            return self.GB20M
        return coord.EarthLocation.of_site(key)

    def __class_getitem__(self, key):
        # For GBT we want to use the GBO official location
        if key == "GBT" or key == "NRAO_GBT":
            return GBT()
        elif key == "GB20M":
            return GB20M()
        return coord.EarthLocation.of_site(key)

    def get_earth_location(longitude, latitude, height=0.0, ellipsoid=None):
        """
        Location on Earth, initialized from geodetic coordinates.

        Parameters
        ----------
        lon : `~astropy.coordinates.Longitude` or float
            Earth East longitude.  Can be anything that initialises an
            `~astropy.coordinates.Angle` object (if float, in degrees).
        lat : `~astropy.coordinates.Latitude` or float
            Earth latitude.  Can be anything that initialises an
            `~astropy.coordinates.Latitude` object (if float, in degrees).
        height : `~astropy.units.quantity.Quantity` ['length'] or float, optional
            Height above reference ellipsoid (if float, in meters; default: 0).
        ellipsoid : str, optional
            Name of the reference ellipsoid to use (default: 'WGS84').
            Available ellipsoids are:  'WGS84', 'GRS80', 'WGS72'.

        Raises
        ------
        `~astropy.units.UnitsError`
            If the units on ``lon`` and ``lat`` are inconsistent with angular
            ones, or that on ``height`` with a length.
        ValueError
            If ``lon``, ``lat``, and ``height`` do not have the same shape, or
            if ``ellipsoid`` is not recognized as among the ones implemented.

        Notes
        -----
        See :meth:`~astropy.coordinates.EarthLocation.from_geodetic`
        """
        return coord.EarthLocation.from_geodetic(longitude, latitude, height, ellipsoid)


def eq2hor(lon, lat, frame, date_obs, unit="deg", location=GBT()):  # noqa: B008
    """
    Equatorial to horizontal coordinate conversion.

    Parameters
    ----------
    lon : float
        Longitude coordinate. E.g., RA or Galactic longitude.
    lat : float
        Latitude coordinate. E.g., Dec or Galactic latitude.
    frame : str
        Input coordinate frame. Must be recognized by `~astropy.coordinates`
    date_obs : str
        Date of observations. Must be a format compatible with `~astropy.time.Time`.
    unit : str
        Units of `lon` and `lat`.
    location : `~astropy.coordinates.EarthLocation`
        Observer location.

    Returns
    -------
    altaz : `~astropy.coordinates.AltAz`
        Horizontal coordinates.

    """

    lonlat = coord.SkyCoord(lon, lat, unit=unit, frame=frame, obstime=Time(date_obs))
    return lonlat.transform_to(coord.AltAz(location=location))


def hor2eq(az, alt, frame, date_obs, unit="deg", location=GBT()):  # noqa: B008
    """
    Horizontal to Equatorial coordinate conversion.

    Parameters
    ----------
    az : float
        Azimuth coordinate.
    alt : float
        Altitude or elevation coordinate.
    frame : str
        Output coordinate frame. Must be recognized by `~astropy.coordinates`
    date_obs : str
        Date of observations. Must be a format compatible with `~astropy.time.Time`.
    unit : str
        Units of `lon` and `lat`.
    location : `~astropy.coordinates.EarthLocation`
        Observer location.

    Returns
    -------
    eq : `~astropy.coordinates.SkyCoord`
        Celestial coordinates in `frame`.

    """

    altaz = coord.SkyCoord(az=az, alt=alt, unit=unit, frame="altaz", obstime=Time(date_obs), location=location)
    return altaz.transform_to(astropy_frame_dict[frame])


def ra2ha(lst, ra):
    """
    Take LST (sec) and RA (deg) and output wrapped HA (hr).
    Follows GBTIDL implementation (with the hour conversion included)
    """
    ha = np.around(15 * (lst / 3600) - ra, 2)
    if ha > 180:
        ha -= 360
    elif ha < -180:
        ha += 360
    return np.around(ha / 15, 2)
