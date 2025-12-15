from astropy.coordinates import EarthLocation, UnknownSiteException

import dysh.coordinates as coords
from dysh.coordinates import GBT, Observatory, astropy_convenience_frame_names, decode_veldef
from dysh.spectra.spectrum import Spectrum


class TestCore:
    """Test dysh.coordinates core functions"""

    def test_observatory(self):
        """Test the Observator class"""
        obs = Observatory()
        for k in ["alma", "Hat Creek", "La Silla Observatory", "Mars Hill", "Whipple"]:
            # Try to get the observatory locations using the existing cache.
            try:
                eloc = EarthLocation.of_site(k)
            # Force cache refresh if this fails.
            except UnknownSiteException:
                eloc = EarthLocation.of_site(k, refresh_cache=True)
            assert obs[k] == eloc
        assert obs["GBT"] is GBT()  # instance method
        assert Observatory["GBT"] is GBT()  # static method
        try:
            Observatory["FOOBAR"]
        except Exception as e:
            assert isinstance(e, UnknownSiteException)

    def test_veldef(self):
        # first make sure we get correct answers for normal inputs
        inputs = [
            "RADILSRK",
            "RADI-LSR",
            "RADILSRD",
            "OPTICMB",
            "VELO-BAR",
            "OPTIBARY",
            "RELATOPO",
            "RADIGALA",
            "OPTI-HEL",
        ]
        outputs = [
            ("radio", "LSRK"),
            ("radio", "LSRK"),
            ("radio", "LSRD"),
            ("optical", "cmb"),
            ("relativistic", "barycentric"),
            ("optical", "barycentric"),
            ("relativistic", "topocentric"),
            ("radio", "galactic"),
            ("optical", "heliocentric"),
        ]
        for i, j in zip(inputs, outputs, strict=False):
            assert decode_veldef(i) == j

        # Now test that bad input raises an exception
        try:
            decode_veldef("This is more than 8 chars")
        except ValueError:
            assert True
        try:
            # frame fails
            decode_veldef("OPTI-LRS")
        except KeyError:
            assert True
        try:
            # convention fails
            decode_veldef("MAXILSRK")
        except KeyError:
            assert True

    def test_veldef_to_convention(self):
        pairs = {"OPTI-HELO": "optical", "VELO-LSR": "radio", "RADI-LSR": "radio", "RELA-LSR": "relativistic"}

        for k, v in pairs.items():
            assert coords.veldef_to_convention(k) == v

    def test_frame_switch(self):
        f = Spectrum.fake_spectrum()
        assert f.velocity_frame == "itrs"
        for k, v in astropy_convenience_frame_names.items():
            fnew = f.with_frame(k)
            assert fnew.velocity_frame == v
