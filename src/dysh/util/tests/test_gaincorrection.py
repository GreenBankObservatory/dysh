import pathlib

import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import Angle
from astropy.table import QTable
from astropy.time import Time
from scipy.optimize import minimize_scalar

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
            assert se == answer[i]

    def test_gaincorrection_factor(self):
        # first test use symmetric angles around 45
        angles = Angle([0.0, 25.0, 40.0, 45.0, 50, 65.0, 90.0], unit=u.degree)
        # answer = np.array([])
        for i in range(len(self.dates)):
            gc1 = self.gbtgc.gain_correction(angles, self.dates[i], zd=True)
            gc2 = self.gbtgc.gain_correction(angles, self.dates[i], zd=False)
            assert all(gc1 == gc2[::-1])
        # now try to reproduce portion of figure 1 in GBT memo 301
        angles = Angle(np.arange(10, 90, 5), unit=u.degree)
        scale = 0.983  # 2009b gmax from Dave Frayer email 12/18/2024
        gc = scale * self.gbtgc.gain_correction(angles, self.dates[4], zd=False)
        # measured using graphreader.com on a screenshot of Figure 1 from GBT memo 301 for "2009b" curve
        answer = np.array(
            [
                0.777,
                0.819,
                0.858,
                0.893,
                0.922,
                0.944,
                0.961,
                0.974,
                0.981,
                0.985,
                0.981,
                0.971,
                0.957,
                0.939,
                0.912,
                0.882,
            ]
        )
        assert np.abs(np.mean(gc - answer)) < 0.001
        # measured using graphreader.com on a screenshot of Figure 1 from GBT memo 301 for "2009a" curve
        gc *= 0.728  # 2009a gmax from Dave Frayer email 12/18/2024.  Yes it is a double scale 0.983*0.728
        answer = np.array(
            [
                0.565,
                0.597,
                0.625,
                0.65,
                0.67,
                0.688,
                0.701,
                0.71,
                0.715,
                0.716,
                0.714,
                0.706,
                0.696,
                0.683,
                0.665,
                0.642,
            ]
        )
        assert np.abs(np.mean(gc - answer)) < 0.001
        # 2014 curve.
        gc = self.gbtgc.gain_correction(angles, self.dates[-1], zd=False)
        # measured using graphreader.com on a screenshot of Figure 1 from GBT memo 301 for "2014" curve
        # harder to do because of graph crowding.
        answer = np.array([0.987, 0.991, 0.994, 0.997, 0.999, 1, 1, 1, 1, 1, 0.998, 0.995, 0.991, 0.987, 0.983, 0.978])
        assert np.abs(np.mean(gc - answer)) < 0.001
        #
        # 2003 curve was not generated using G(ZD) so we do not test it here.
        # according to Frayer it is (IDL code):
        # ;;No Gravity corrections Condon 2003 PTCS/PN/31.1
        # ;;max eta_a=43% (390um rms) scaled to max of 60% fac=0.43/0.60=0.717
        # nograv=0.717*exp(-1.0*((xel-52.0)/54.)^2.)

    def test_aperture_efficiency(self):
        """This attempts to reproduce Table 2 in GBT Memo 301"""
        freqs = np.array([10.0, 30.0, 43.0, 80.0, 110.0]) * u.GHz
        answer = np.array(
            [
                [0.69, 0.56, 0.43, 0.13, 0.03],  # 2003
                [0.69, 0.56, 0.43, 0.13, 0.03],  # 2009a
                [0.70, 0.65, 0.59, 0.37, 0.21],  # 2009b
                [0.70, 0.65, 0.60, 0.39, 0.23],  # 2014
            ]
        )
        f = self.gbtgc.gain_correction
        i = 0
        for d in self.dates[[5, 7]]:
            # first find the elevation angle where the gain curve reaches a maximum
            maxpoint = minimize_scalar(
                lambda x: -f(angle=x * u.degree, date=d, zd=False), bounds=[0, 90], method="bounded"
            )
            # Evaluate the aperture efficiency at the given requencies and the elevation of the gain maximum
            a = maxpoint.x * u.degree
            if i < 2:
                # For 2003 and 2009a, a smaller surface error is used in memo 301 than the standard in our gain table.
                # So pass in the memo value in these cases.
                ap = self.gbtgc.aperture_efficiency(freqs, a, d, zd=False, eps0=390 * u.micron)
            else:
                ap = self.gbtgc.aperture_efficiency(freqs, a, d, zd=False)
            apr = np.round(ap, decimals=2)  # this is all the precision the table has
            # print(np.abs(apr - answer[i]))
            assert all(apr == answer[i])
            i += 1

    def test_airmass(self):
        angles = Angle([0.0, 25.0, 40.0, 45.0, 50, 65.0, 90.0], unit=u.degree)
        answer = np.array([37.5541407, 2.35964333, 1.5501963, 1.40793864, 1.29840534, 1.09473653, 0.99060048])
        am = self.gbtgc.airmass(angles, zd=False)
        assert np.mean(am - answer) < 1e-10
        assert all(self.gbtgc.airmass(angles, zd=True) == am[::-1])
