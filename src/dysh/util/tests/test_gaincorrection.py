import pathlib

import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import Angle
from astropy.table import QTable
from astropy.time import Time
from scipy.optimize import minimize_scalar

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
        answer = np.array([490.0, 490.0, 440.0, 418.0, 240.0, 240.0, 230.0, 230.0]) * u.micron
        for i in range(len(self.dates)):
            se = self.gbtgc.surface_error(self.dates[i])
            # from print(f"doing {self.dates[i]}")
            assert se == answer[i]

    def test_gaincorrection_factor(self):
        angles = Angle([0.0, 5.0, 25.0, 40.0, 45.0, 50, 65.0, 85.0, 90.0], unit=u.degree)
        answer = np.array([])
        for i in range(len(self.dates)):
            gc1 = self.gbtgc.gain_correction(angles, self.dates[i], zd=True)
            gc2 = self.gbtgc.gain_correction(angles, self.dates[i], zd=False)
            assert all(gc1 == gc2[::-1])

    def test_aperture_efficiency(self):
        """This attempts to reproduce 2009b and 2014 rows of Table 2 in GBT Memo 301"""
        freqs = np.array([10.0, 30.0, 43.0, 80.0, 110.0]) * u.GHz
        answer = np.array(
            [
                [0.70, 0.65, 0.59, 0.37, 0.21],
                [0.70, 0.65, 0.60, 0.39, 0.23],
            ]
        )
        # test only for 2009b and 2014
        f = self.gbtgc.gain_correction
        i = 0
        for d in self.dates[[5, 7]]:
            # first find the elevation angle where the gain curve reaches a maximum
            maxpoint = minimize_scalar(
                lambda x: -f(angle=x * u.degree, date=d, zd=False), bounds=[0, 90], method="bounded"
            )
            # Evaluate the aperture efficiency at the given requencies and the elevation of the gain maximum
            a = maxpoint.x * u.degree
            ap = self.gbtgc.aperture_efficiency(freqs, a, d, zd=False)
            apr = np.round(ap, decimals=2)  # this is all the precision the table has
            # if they differ by 1% efficiency, that's acceptable for this test.
            assert np.max(np.abs(apr - answer[i])) < 0.011
            i += 1

    def test_airmass(self):
        angles = Angle([0.0, 25.0, 40.0, 45.0, 50, 65.0, 90.0], unit=u.degree)
        answer = np.array([37.5541407, 2.35964333, 1.5501963, 1.40793864, 1.29840534, 1.09473653, 0.99060048])
        am = self.gbtgc.airmass(angles, zd=False)
        assert np.mean(am - answer) < 1e-10
        assert all(self.gbtgc.airmass(angles, zd=True) == am[::-1])
