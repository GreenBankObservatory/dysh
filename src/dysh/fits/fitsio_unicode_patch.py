"""
Monkey patch for fitsio to handle non-ASCII characters in FITS string columns.

The fitsio library uses ASCII decoding for byte strings in FITS files, which fails
when encountering non-ASCII bytes (e.g., byte 0x96 = en-dash in Windows-1252).
This module patches fitsio to fall back to latin-1 encoding, which can decode
any byte value 0-255.

Usage:
    import dysh.fits.fitsio_unicode_patch  # Just import to apply the patch

Or explicitly:
    from dysh.fits.fitsio_unicode_patch import patch_fitsio_unicode
    patch_fitsio_unicode()

This patch is automatically applied when dysh.fits.sdfitsload is imported.
"""

import sys

import numpy as np

from dysh.log import logger

_patch_applied = False


def patch_fitsio_unicode():
    """
    Monkey-patch fitsio to handle non-ASCII characters in FITS string columns.

    Replaces the ASCII decoding with latin-1 fallback, which allows reading
    FITS files that contain non-ASCII bytes in string columns (e.g., GB20M data).

    This is safe to call multiple times (only applies once).
    """
    global _patch_applied

    if _patch_applied:
        return

    try:
        from fitsio.hdu.table import TableHDU
    except ImportError:
        logger.warning("fitsio not installed, Unicode patch not applied")
        return

    def patched_decode(array, copy_if_needed=True):
        """
        Decode FITS byte strings to Unicode using latin-1 encoding.

        This is a patched version of TableHDU._maybe_decode_fits_ascii_strings_to_unicode_py3
        that uses latin-1 encoding instead of ASCII. latin-1 is a superset of ASCII that can
        decode any byte value (0-255), preventing UnicodeDecodeError on non-ASCII characters.

        Note: latin-1 gives identical results to ASCII for bytes 0-127, so this works
        correctly for all standard FITS files while also handling non-standard files.
        """
        IS_PY3 = sys.version_info[0] >= 3

        if not IS_PY3:
            return array

        # Check if any string columns need conversion
        do_conversion = False
        new_dt = []

        for _dt in array.dtype.descr:
            if "S" in _dt[1]:  # Byte string column (dtype like '|S32')
                do_conversion = True
                # Convert dtype from bytes (S) to unicode (U)
                if len(_dt) == 3:  # Has shape tuple
                    new_dt.append(
                        (
                            _dt[0],
                            _dt[1].replace("S", "U").replace("|", ""),
                            _dt[2],
                        )
                    )
                else:
                    new_dt.append((_dt[0], _dt[1].replace("S", "U").replace("|", "")))
            else:
                new_dt.append(_dt)

        if do_conversion:
            # Always use latin-1 encoding (superset of ASCII, handles bytes 0-255)
            # Manually convert each field
            new_array = np.empty(len(array), dtype=new_dt)

            for field_name in array.dtype.names:
                if array.dtype[field_name].kind == "S":  # Byte string field
                    # Decode byte strings with latin-1 (handles any byte 0-255)
                    new_array[field_name] = np.array(
                        [
                            s.decode("latin-1", errors="replace") if isinstance(s, bytes) else s
                            for s in array[field_name]
                        ]
                    )
                else:
                    # Copy non-string fields directly
                    new_array[field_name] = array[field_name]

            array = new_array

        return array

    # Apply the monkey patch
    TableHDU._maybe_decode_fits_ascii_strings_to_unicode_py3 = staticmethod(patched_decode)

    logger.debug("Applied fitsio Unicode patch (latin-1 fallback)")
    _patch_applied = True


# Automatically apply patch on import
patch_fitsio_unicode()
