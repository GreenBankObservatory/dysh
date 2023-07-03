
import os
import pytest

from astropy.utils.data import (
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

