import astropy.units as u
import numpy as np

from dysh.fits import gbtfitsload


class TestVelocityFrames:
    def read_ascii(self, filename):
        values = np.loadtxt(filename, skiprows=3, unpack=True)
        # print(np.shape(values))
        freqGHz = values[:][0] * u.GHz
        flux = values[:][1] * u.ct
        return freqGHz, flux

    def compare_gbtidl(self, spectrum, filename, frame, maxdiff):
        freqGHz, flux = self.read_ascii(filename)
        # use version of method that makes a copy, so it gets some test coverage.
        s = spectrum.with_frame(frame)
        x = np.mean(freqGHz - s.spectral_axis)
        y = np.nanmean(flux - s.flux)
        # print(x, x.to("Hz"))
        # print(y)
        assert np.abs(y.value) < 1
        assert np.abs(x.to("Hz").value) < maxdiff

    def test_compare_with_GBTIDL(self, data_dir):
        """ """
        # get filenames
        sdf_file = f"{data_dir}/TGBT21A_501_11/TGBT21A_501_11.raw.156.fits"

        sdf = gbtfitsload.GBTFITSLoad(sdf_file)
        sp = sdf.gettp(scan=156, fdnum=0, plnum=0, ifnum=0)[0].total_power(0)

        # dict of veldef frame type to astropy frame type string
        framedict = {
            "BAR": "icrs",
            "GAL": "galactocentric",
            "GEO": "gcrs",
            "HEL": "hcrs",
            "LSD": "lsrd",
            "LSR": "lsrk",
            "TOPO": sp.observer,  # must be an actual coordinate for topographical
        }
        maxHzdiff = {"BAR": 2.0, "GAL": 13000.0, "GEO": 2.0, "HEL": 2.0, "LSD": 30.0, "LSR": 2.0, "TOPO": 2.0}
        for k, v in framedict.items():
            gbtidl_file = f"{data_dir}/gbtidl_spectra/onoff-L_gettp_156_intnum_0_{k}.ascii"
            self.compare_gbtidl(filename=gbtidl_file, spectrum=sp, frame=v, maxdiff=maxHzdiff[k])
