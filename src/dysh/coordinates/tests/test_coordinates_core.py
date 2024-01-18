from astropy.coordinates import EarthLocation, UnknownSiteException

from dysh.coordinates import GBT, Observatory, decode_veldef


class TestCore:
    """Test dysh.coordinates core functions"""

    def test_observatory(self):
        """Test the Observator class"""
        obs = Observatory()
        for k in ["alma", "Hat Creek", "La Silla Observatory", "Mars Hill", "Whipple"]:
            assert obs[k] == EarthLocation.of_site(k)
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
        for i, j in zip(inputs, outputs):
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
