"""
Core functions for FITS/SDFITS
"""

import astropy.coordinates as coord

# Velocity frame conventions, partially stolen from pyspeckit.
# See also Section6.2.5.2 of the GBT observer's guide https://www.gb.nrao.edu/scienceDocs/GBTog.pdf

# Todo:
# Deal with cmb which will require custom frame
# See https://docs.astropy.org/en/stable/coordinates/spectralcoord.html#defining-custom-velocity-frames
# Also deal with 'rest' frame
astropy_frame_dict = {
    "VLSR": coord.LSRK,
    "VRAD": coord.LSRK,
    "VELO": coord.LSRK,
    "VOPT": coord.LSRK,
    "LSRD": coord.LSRD,
    "LSRK": coord.LSRK,
    "-LSR": coord.LSRK,
    "-HEL": coord.HCRS,
    "-BAR": coord.ICRS,
    "BAR": coord.ICRS,
    "BARY": coord.ICRS,
    "barycentric": coord.ICRS,
    "VHEL": coord.HCRS,
    "heliocentric": coord.HCRS,
    "VGEO": coord.GCRS,
    "GALAC": coord.Galactic,
    "-GALA": coord.Galactic,
    "topocentric": coord.ITRS,  # but need to add observatory position
}

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
    "-OBS": "obs",
    "VHEL": "heliocentric",
    "VGEO": "topocentric",
    "TOPO": "topocentric",
    "VREST": "rest",
    "Z": "rest",
    "FREQ": "rest",
    "WAV": "rest",
    "WAVE": "rest",
    "CMB": "cmb",
    "GALAC": "galactic",
    "GALA": "galactic",
    "ALAC": "galactic",
}

# Dictionary to convert from FITS velocity convention to specutils string.
# At GBT, VELO was written by sdfits filler for some unknown amount of
# time instead of RELA, so allow for it here
vconv_dict = {
    "OPTI": "optical",
    "RADI": "radio",
    "RELA": "relativistic",
    "VELO": "relativistic",
}


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
        velocity_convention = vconv_dict[vconv]
    except KeyError:
        raise KeyError(f"Velocity convention {vconv} not recognized.")

    frame = veldef[4:]
    try:
        frame_type = frame_dict[frame]
    except KeyError:
        raise KeyError(f"Velocity frame {frame} not recognized.")

    return velocity_convention, frame_type
