import os
import glob
import pytest
import pathlib

import numpy as np

from astropy.io import fits

import dysh
from dysh import util
from dysh.fits import gbtfitsload

dysh_root = pathlib.Path(dysh.__file__).parent.resolve()


class TestGBTFITSLoad:
    """ """

    def setup_method(self):
        self.root_dir = util.get_project_root()
        self.data_dir = f"{self.root_dir}/testdata"
        self._file_list = glob.glob(f"{self.data_dir}/TGBT21A_501_11/*.fits")

    def test_load(self):
        expected = {
            "TGBT21A_501_11.raw.vegas.fits": 4,
            "TGBT21A_501_11_getps_scan_152_intnum_0_ifnum_0_plnum_0.fits": 1,
            "TGBT21A_501_11_gettp_scan_152_intnum_0_ifnum_0_plnum_0_cal_state_0.fits": 1,
            "TGBT21A_501_11_gettp_scan_152_intnum_0_ifnum_0_plnum_0_cal_state_1.fits": 1,
            "TGBT21A_501_11_ifnum_0_int_0-2.fits": 24,
            "TGBT21A_501_11_ifnum_0_int_0-2_getps_152_plnum_0.fits": 1,
            "TGBT21A_501_11_ifnum_0_int_0-2_getps_152_plnum_1.fits": 1,
            "TGBT21A_501_11_gettp_scan_152_ifnum_0_plnum_0.fits": 1,
        }

        for fnm in self._file_list:
            print(fnm)

            filename = os.path.basename(fnm)
            sdf = gbtfitsload.GBTFITSLoad(fnm)
            assert len(sdf.index(bintable=0)) == expected[filename]

    def test_getps_single_int(self):
        """ """

        gbtidl_file = f"{self.data_dir}/TGBT21A_501_11/TGBT21A_501_11_getps_scan_152_intnum_0_ifnum_0_plnum_0.fits"
        # We should probably use dysh to open the file...
        hdu = fits.open(gbtidl_file)
        gbtidl_getps = hdu[1].data["DATA"][0]

        sdf_file = f"{self.data_dir}/TGBT21A_501_11/TGBT21A_501_11.raw.vegas.fits"
        sdf = gbtfitsload.GBTFITSLoad(sdf_file)
        # psscan is a ScanList
        psscan = sdf.getps(152)
        assert len(psscan) == 1
        psscan.calibrate()
        dysh_getps = psscan[0].calibrated(0).flux.to("K").value

        diff = gbtidl_getps - dysh_getps
        assert np.nanmedian(diff) == 0.0
        assert np.all(abs(diff[~np.isnan(diff)]) < 5e-7)
        assert np.isnan(diff[3072])

    def test_gettp_single_int(self):
        """ """

        # Get the answer from GBTIDL.
        gbtidl_file = f"{self.data_dir}/TGBT21A_501_11/TGBT21A_501_11_gettp_scan_152_intnum_0_ifnum_0_plnum_0_cal_state_1.fits"
        hdu = fits.open(gbtidl_file)
        gbtidl_gettp = hdu[1].data["DATA"][0]

        # Get the answer from dysh.
        sdf_file = f"{self.data_dir}/TGBT21A_501_11/TGBT21A_501_11.raw.vegas.fits"
        sdf = gbtfitsload.GBTFITSLoad(sdf_file)
        tps_on = sdf.gettp(152, sig=True, cal=True, calibrate=False)
        assert len(tps_on) == 1

        # Compare.
        diff = tps_on[0].total_power(0).flux.value - gbtidl_gettp
        assert np.nanmean(diff) == 0.0

        # Now with the noise diode Off.
        tps_off = sdf.gettp(152, sig=True, cal=False, calibrate=False)
        assert len(tps_off) == 1
        gbtidl_file = f"{self.data_dir}/TGBT21A_501_11/TGBT21A_501_11_gettp_scan_152_intnum_0_ifnum_0_plnum_0_cal_state_0.fits"
        hdu = fits.open(gbtidl_file)
        gbtidl_gettp = hdu[1].data["DATA"][0]
        diff = tps_off[0].total_power(0).flux.value - gbtidl_gettp
        assert np.nanmean(diff) == 0.0

        # Now, both on and off.
        tps = sdf.gettp(152, sig=True, cal=True)
        assert len(tps) == 1
        tps_tavg = tps.timeaverage()
        assert len(tps_tavg) == 1
        gbtidl_file = f"{self.data_dir}/TGBT21A_501_11/TGBT21A_501_11_gettp_scan_152_ifnum_0_plnum_0.fits"
        hdu = fits.open(gbtidl_file)
        table = hdu[1].data
        spec = table["DATA"][0]
        diff = tps[0].total_power(0).flux.value - spec
        assert np.nanmean(diff) == 0.0
        # what about tps_tavg

    def test_load_multifits(self):
        """
        Loading multiple SDFITS files under a directory.
        """

        fits_path = f"{self.data_dir}/AGBT18B_354_03.raw.vegas"

        sdf = gbtfitsload.GBTFITSLoad(fits_path)
