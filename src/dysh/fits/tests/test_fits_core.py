from dysh.fits import decode_veldef


class TestCore:
    """Test dysh.fits core functions"""

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
            ("doppler_radio", "LSRK"),
            ("doppler_radio", "LSRK"),
            ("doppler_radio", "LSRD"),
            ("doppler_optical", "cmb"),
            ("doppler_relativistic", "barycentric"),
            ("doppler_optical", "barycentric"),
            ("doppler_relativistic", "topocentric"),
            ("doppler_radio", "galactic"),
            ("doppler_optical", "heliocentric"),
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
