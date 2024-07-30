import numpy as np

import dysh.util as du


class TestUtil:
    """Test dysh.fits core functions"""

    def test_sq_weighted_avg(self):
        """Test square weighted average function"""
        a = np.random.rand(1024)
        w = np.random.rand(1024)
        u = 0
        ws = 0
        for i in range(len(a)):
            u += a[i] * a[i] * w[i]
            ws += w[i]
        result = np.sqrt(u / ws)
        diff = result - du.sq_weighted_avg(a, 0, w)
        assert np.abs(diff) < 2e-15

    def test_match(self):
        """Test minimum_string_match function"""
        s = ["alpha", "beta", "gamma", "gemma"]
        assert du.minimum_string_match("a", s) == "alpha"
        assert du.minimum_string_match("A", s) == None
        assert du.minimum_string_match("g", s) == None
        assert du.minimum_string_match("ga", s) == "gamma"
        assert du.minimum_string_match("am", s) == None

    def test_powerof2(self):
        """Test powerof2 function"""
        inout = {2**0: 0, 2**15: 15, 2**15.49: 15, 2**15.5: 16}

        for k, v in inout.items():
            assert du.powerof2(k) == v
