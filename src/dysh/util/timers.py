#! /usr/bin/env python
#
import os
import sys
import time

from pathlib import Path
import numpy as np
from astropy.time import Time

from astropy.table import Table,vstack
import astropy.units as u



__all__ = ["DTime"]


class DTime(object):
    r"""This class encapsulated some popular timing/performance tools.

    dt = DTime()
    dt.tag("test1")
    dt.tag("test2")
    dt.tag("test3")
    dt.done()

    By default it simply builds a delta-time of the time it took between the different tags, as 
    labeled by their tag name. If DTime() is supplied a number of data items for extra columns,
    these will be reported, or stored in a table, if out= is supplied.
    

    """

    def __init__(self, benchmark="generic", units="ms",
                 out = None,  mode="overwrite", #   append, "w" "w!", "a"
                 data_cols = None, data_units = None, data_types = None):
        #   @todo  add table entries for future data!=None
        #   @todo  add keyword so this class can be de-activated
        self.active = 1
        self.state = 0
        if data_cols is not None and data_units is not None and data_types is not None:
            # @todo check if all lengths are the same
            self.table = Table(names=data_cols, meta={"name": f"Dysh Benchmark {benchname}"}, units=data_units, dtype=data_types)
        else:
            print("Warning: skipping table saving")
            self.table = None
        
        self.stats = []
        #   also start the first timer (the extra TBD)
        self.stats.append(["start", time.perf_counter_ns(), [0, 0, 0, 0, ""]])

    def tag(self, name, data=None):
        self.stats.append([name, time.perf_counter_ns(), data])

    def close(self, data=None):
        self.state = 1
        self.stats.append(["end", time.perf_counter_ns(), data])

    def report(self, debug=False):
        assert(self.state == 1)
        if debug:
            print(self.stats)
        n = len(self.stats)
        #  for now just print the timings
        for i in range(1,n):
            dt = (self.stats[i][1] - self.stats[i-1][1])/1e6   # now in ms, check units
            print(i,self.stats[i][0],dt)
        if self.table is not None:
            print("@todo process and write the table")
        
if __name__ == "__main__":
    #  Also compare the output of this with /usr/bin/time on the executable
    #  to find any overhead of the timer we didn't account. Optionally run
    #  a benchmark once and twice and use the difference.
    dt = DTime()
    dt.tag("nothing    ")
    dt.tag("test0      ")
    dt.tag("test1      ")
    dt.tag("test2      ")
    dt.tag("test3      ")
    dt.tag("test4      ")
    a = np.arange(1e3)
    dt.tag("arange(1e3)")
    a = np.arange(1e4)
    dt.tag("arange(1e4)")
    a = np.arange(1e5)
    dt.tag("arange(1e5)")
    a = np.arange(1e6)
    dt.tag("arange(1e6)")
    a = np.arange(1e7)
    dt.tag("arange(1e7)")
    a = np.arange(1e8)
    dt.tag("arange(1e8)")
    a = np.arange(1e9)
    dt.tag("arange(1e9)")
    dt.close()
    dt.report()
