from gzip import GzipFile

import astropy.units as u
import numpy as np

from dysh.fits import gbtfitsload


class TestDopplerConvention:
    def read_ascii(self, filename, unzip=True):
        if unzip:
            g = GzipFile(filename)
        else:
            g = filename
        values = np.loadtxt(g, skiprows=3, unpack=True)
        velokms = values[:][0] * u.km / u.s
        flux = values[:][1] * u.ct
        return velokms, flux

    def compare_gbtidl(self, spectrum, filename, doppler_convention, maxdiff):
        velokms, _flux = self.read_ascii(filename)
        s = spectrum.with_velocity_convention(doppler_convention)
        x = np.mean(velokms - s.axis_velocity())
        assert np.abs(x.value) < maxdiff

    def test_compare_with_GBTIDL(self, data_dir):
        """ """
        # get filenames
        sdf_file = f"{data_dir}/TGBT21A_501_11/TGBT21A_501_11.raw.vegas.fits"

        sdf = gbtfitsload.GBTFITSLoad(sdf_file)
        sp = sdf.getps(scan=152, fdnum=0, ifnum=0, plnum=0)[0].getspec(0)
        sp.set_frame("hcrs")  # data were tracked in heliocentric. change in place.

        conventions = {"OPTI-HEL": "optical", "RADI-HEL": "radio", "TRUE-HEL": "relativistic"}
        maxdiff = 1.0
        for k, v in conventions.items():
            gbtidl_file = f"{data_dir}/gbtidl_spectra/onoff-L_getps_152_{k}.ascii.gz"
            self.compare_gbtidl(filename=gbtidl_file, spectrum=sp, doppler_convention=v, maxdiff=maxdiff)
