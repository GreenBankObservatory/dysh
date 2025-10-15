import time

from dysh.util import timers


class TestUtil:
    """Test dysh.util files functions"""

    def test_dysh_timers(self):
        """Test dysh timers"""
        slop = 1.2  # one of windows OS complained 49.8838 > 50
        n_ms = 50
        dt = timers.DTime()
        dt.tag("test1")
        time.sleep(n_ms / 1000)  # 50 ms
        dt.tag("test2")
        dt.report()
        dt.close()
        dt_total = dt.total()
        assert dt_total * slop > n_ms
