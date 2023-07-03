
import os
import pytest
import numpy as np

from astropy.io import fits
from astropy.utils.data import (
                    get_pkg_data_filename,
                    get_pkg_data_filenames,
                    )

from dysh.fits import gbtfitsload


class TestGBTFITSLoad():
    """
    """

    def setup_method(self):
        self._file_list = list(get_pkg_data_filenames("data/", pattern="*.fits"))

    def test_load(self):

        expected = {"TGBT21A_501_11.raw.vegas.fits": 4,
                    "TGBT21A_501_11_getps_scan_152_intnum_0_ifnum_0_plnum_0.fits": 1,
                    }

        for fnm in self._file_list:
            print(fnm)

            filename = os.path.basename(fnm)
            sdf = gbtfitsload.GBTFITSLoad(fnm)
            assert len(sdf._ptable[0]) == expected[filename]


    def test_getps(self):
        """
        """

        gbtidl_file = get_pkg_data_filename("data/TGBT21A_501_11_getps_scan_152_intnum_0_ifnum_0_plnum_0.fits")
        # We should probably use dysh to open the file...
        hdu = fits.open(gbtidl_file)
        gbtidl_getps = hdu[1].data["DATA"][0]

        sdf_file = get_pkg_data_filename("data/TGBT21A_501_11.raw.vegas.fits")
        sdf = gbtfitsload.GBTFITSLoad(sdf_file)
        psscan = sdf.getps(152)
        psscan.calibrate()
        dysh_getps = psscan.calibrated(0).flux.to("K").value

        diff = gbtidl_getps - dysh_getps
        assert np.nanmedian(diff) == 0.0
        assert np.all(abs(diff[~np.isnan(diff)]) < 5e-7)
        assert np.isnan(diff[3072])
