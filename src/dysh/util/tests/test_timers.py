
import time
from dysh.util import timers 


class TestUtil:
    """Test dysh.util files functions"""

    def test_dysh_timers(self):
        """Test dysh timers"""
        n_ms = 50
        dt = timers.DTime()
        dt.tag("test1")
        time.sleep(n_ms/1000)   # 50 ms
        dt.tag("test2")
        dt.report()
        dt.close()
        dt_total = dt.total()
        assert(dt_total > n_ms)
