import numpy as np
import pytest

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
        s = ["alpha", "beta", "gamma", "gemma", "beta2"]
        assert du.minimum_string_match("a", s) == "alpha"
        assert du.minimum_string_match("A", s) == None  # noqa: E711
        assert du.minimum_string_match("g", s) == "gamma"
        assert du.minimum_string_match("ga", s) == "gamma"
        assert du.minimum_string_match("am", s) == None  # noqa: E711
        assert du.minimum_string_match("beta", s) == "beta"
        assert du.minimum_string_match("BEt", s, casefold=True) == "beta"
        assert du.minimum_string_match("gem", s, casefold=True) == "gemma"
        assert du.minimum_string_match("gem", list(map(str.upper, s)), casefold=True) == "GEMMA"

        assert du.minimum_list_match(["b", "al"], s, casefold=False) == ["beta", "alpha"]
        assert du.minimum_list_match(["bE", "ALP"], s, casefold=True) == ["beta", "alpha"]
        # ensure single string is treated correctly
        assert du.minimum_list_match("bag", s, casefold=False) is None
        assert du.minimum_list_match("be", s, casefold=False) == ["beta"]

    def test_powerof2(self):
        """Test powerof2 function"""
        inout = {2**0: 0, 2**15: 15, 2**15.49: 15, 2**15.5: 16}
        for k, v in inout.items():
            assert du.powerof2(k) == v

    def test_merge_ranges(self):
        """Test merge_ranges function"""
        from astropy import units as u

        r = [(1 * u.GHz, 2 * u.GHz), (1.5 * u.GHz, 3 * u.GHz)]
        assert list(du.merge_ranges(r)) == [(1 * u.GHz, 3 * u.GHz)]
        assert list(du.merge_ranges([])) == []

    def test_get_valid_channel_range(self):
        a = [100, 200]
        assert du.get_valid_channel_range(a) == a
        a = [400, 100]
        with pytest.raises(ValueError):
            du.get_valid_channel_range(a)
        a = [1, [24, 99], 33]
        with pytest.raises(ValueError):
            du.get_valid_channel_range(a)
        a = [4567, 9678]
        assert du.get_valid_channel_range(np.array(a)) == a
        b = [a]
        assert du.get_valid_channel_range(b) == a
        assert du.get_valid_channel_range(np.array(b)) == a

    def test_isot_to_mjd(self):
        isot = [
            "2024-10-11T06:49:29.00",
            "2024-10-11T06:49:31.01",
            "2024-10-11T06:49:34.01",
            "2024-10-11T06:49:36.01",
            "2024-10-11T06:49:54.00",
            "2024-10-11T06:49:56.01",
            "2024-10-11T06:49:59.01",
            "2024-10-11T06:50:01.01",
            "2024-10-11T07:31:52.51",
            "2024-10-11T07:31:53.51",
            "2024-10-11T07:31:54.51",
            "2024-10-11T07:31:55.51",
            "2024-10-11T07:31:57.51",
            "2024-10-11T07:31:58.51",
            "2024-10-11T07:31:59.51",
            "2024-10-11T07:32:00.51",
            "2024-10-11T07:32:18.51",
            "2024-10-11T07:32:19.51",
            "2024-10-11T07:32:20.51",
            "2024-10-11T07:32:21.51",
        ]
        expected = np.array(
            [
                60594.28436343,
                60594.28438669,
                60594.28442141,
                60594.28444456,
                60594.28465278,
                60594.28467604,
                60594.28471076,
                60594.28473391,
                60594.3138022,
                60594.31381377,
                60594.31382535,
                60594.31383692,
                60594.31386007,
                60594.31387164,
                60594.31388322,
                60594.31389479,
                60594.31410312,
                60594.3141147,
                60594.31412627,
                60594.31413785,
            ]
        )
        result = du.isot_to_mjd(isot)
        assert np.all(abs(result - expected) * 24 * 3600 < 1e-3)  # Less than 1 ms difference.
