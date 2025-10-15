"""
Tests for dysh.util.calibrator
"""

import pytest
from astropy import units as u

from dysh.util import calibrator


class TestCalibratorTable:
    def test_valid_scales(self):
        ct = calibrator.CalibratorTable()
        assert ct.valid_scales() == ["Perley-Butler 2017", "Ott 1994"]


class TestCalibrator:
    def test_from_name(self):
        # Loop over valid names.
        # Scale aware.
        ct = calibrator.CalibratorTable()
        for s in ct.valid_scales():
            for name in ct.valid_names(s):
                c = calibrator.Calibrator.from_name(name, scale=s)
        # Any scale.
        for name in ct.valid_names():
            c = calibrator.Calibrator.from_name(name)

        # Test some names.
        c = calibrator.Calibrator.from_name("3C123")
        assert c.name == "3C123"
        assert c.las == 1.2e-2
        assert c.cal is True
        c = calibrator.Calibrator.from_name("3C286", scale="Ott 1994")
        assert c.name == "3C286"
        assert c.las == 8.3e-4
        assert c.cal is True
        # Lowercase scale.
        cl = calibrator.Calibrator.from_name("3C286", scale="ott 1994")
        assert cl.name == "3C286"
        assert cl.las == c.las
        assert cl.cal is True
        assert cl.coefs == c.coefs
        # Typo in scale name.
        cl = calibrator.Calibrator.from_name("3C286", scale="OTt 1994")
        assert cl.name == "3C286"
        assert cl.las == c.las
        assert cl.cal is True
        assert cl.coefs == c.coefs
        c = calibrator.Calibrator.from_name("Hydra A")
        assert c.name == "Hydra A"
        assert c.las == 0.13
        assert c.cal is False
        # Lowercase name.
        cl = calibrator.Calibrator.from_name("hydra a")
        assert cl.name == "hydra a"
        assert cl.las == c.las
        assert cl.cal is False
        assert cl.coefs == c.coefs
        # Short scale name.
        cs = calibrator.Calibrator.from_name("hydra A", scale="Perley")
        assert cs.name == "hydra A"
        assert cs.las == c.las
        assert cs.cal is False
        assert cs.coefs == c.coefs

        # Invalid calibrator name.
        with pytest.raises(ValueError):
            c = calibrator.Calibrator.from_name("Lemming A")
        # Invalid scale.
        with pytest.raises(ValueError):
            c = calibrator.Calibrator.from_name("3C123", scale="Salas-Pound 2025")

    def test_get_snu(self):
        c = calibrator.Calibrator.from_name("3C123")
        snu = c.compute_sed(1 * u.GHz)
        assert snu == 63.34320007697665 * u.Jy

        # Warnings if outside frequency range.
        with pytest.warns(UserWarning):
            c.compute_sed(100 * u.GHz)
        with pytest.warns(UserWarning):
            c.compute_sed(100 * u.Hz)

        c = calibrator.Calibrator.from_name("3C286", scale="Ott 1994")
        with pytest.warns(UserWarning):
            snu = c.compute_sed(1 * u.GHz)
        assert snu == 16.91998602503004 * u.Jy

    def test_Calibrator(self):
        c = calibrator.Calibrator("3", "l", [1, 1, 1], calibrator.poly_ott, nu_max=1 * u.m)
        assert c.name == "3"
        assert c.scale == "l"
        assert c.coefs == [1, 1, 1]
        assert c.nu_min is None
        assert c.nu_max == 1 * u.m
        with pytest.warns(UserWarning):
            snu = c.compute_sed(100 * u.GHz)
        assert snu.value == pytest.approx(1e31)
        assert snu.unit == u.Jy

        with pytest.raises(TypeError):
            calibrator.Calibrator("3", "l", [1, 1, 1], calibrator.poly_ott, nu_max=1)
