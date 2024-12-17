import pathlib

import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import Angle
from astropy.table import QTable
from astropy.time import Time

import dysh.util as util
from dysh.util.gaincorrection import GBTGainCorrection


class TestGainCorrection:
    """Test the GBT Gain Correction functions"""

    def setup_method(self):
        self.gbtgc = GBTGainCorrection()
        self.angles = Angle(np.arange(91), unit=u.degree)
        self.dates = Time(
            [
                "2001-02-03",
                "2002-02-28",
                "2004-07-07",
                "2008-08-04",
                "2009-11-01",
                "2014-04-30",
                "2014-05-01",
                "2018-06-14",
            ]
        )

    def test_surface_error(self):
        answer = np.array([490.0, 490.0, 440.0, 418.0, 230.0, 230.0, 230.0, 230.0]) * u.micron
        for i in range(len(self.dates)):
            se = self.gbtgc.surface_error(self.dates[i])
            assert se == answer[i]

    def test_gaincorrection_factor(self):
        angles = Angle([0.0, 25.0, 45.0, 60.0, 77.0, 90.0], unit=u.degree)
        for i in range(len(self.dates)):
            gc = self.gbtgc.gain_correction(angles, self.dates[i], zd=True)
            gc = self.gbtgc.gain_correction(angles, self.dates[i], zd=False)

    def test_aperture_efficiency(self):
        pass

    def test_airmass(self):
        angles = Angle([0.0, 25.0, 40.0, 45.0, 50, 65.0, 90.0], unit=u.degree)
        answer = np.array([37.5541407, 2.35964333, 1.5501963, 1.40793864, 1.29840534, 1.09473653, 0.99060048])
        am = self.gbtgc.airmass(angles, zd=False)
        assert np.mean(am - answer) < 1e-10
        assert all(self.gbtgc.airmass(angles, zd=True) == am[::-1])
