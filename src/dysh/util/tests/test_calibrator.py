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
        # Test some names.
        c = calibrator.Calibrator.from_name("3C123")
        assert c.name == "3C123"
        assert c.las == 1.2e-2
        assert c.cal is True
        c = calibrator.Calibrator.from_name("3C286", scale="Ott 1994")
        assert c.name == "3C286"
        assert c.las == 8.3e-4
        assert c.cal is True
        c = calibrator.Calibrator.from_name("Hydra A")
        assert c.name == "Hydra A"
        assert c.las == 0.13
        assert c.cal is False

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

        c = calibrator.Calibrator.from_name("3C286", scale="Ott 1994")
        snu = c.compute_sed(1 * u.GHz)
        assert snu == 16.91998602503004 * u.Jy
