#! /usr/bin/env python
#
import os
import sys
import time
import pstats
from pstats import SortKey
import cProfile

from pathlib import Path
import numpy as np
from astropy.time import Time

from astropy.table import Table,vstack
import astropy.units as u



__all__ = ["DTime"]

# see also ADMIT's util.utils.Dtime()

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

    def __init__(self,
                 benchname="generic", units="ms",
                 out = None, overwrite=False, append=False, profile=False, statslines=25,
                 data_cols = None, data_units = None, data_types = None):
        self.benchname = benchname
        self.active = 1                    # @todo
        self.state = 0
        self.out = out                     # args
        self.append = append               # args
        self.overwrite = overwrite         # args
        self.profile = profile             # args
        self.statslines = int(statslines)  # args
        if self.profile:
            self.pr = cProfile.Profile()
            self.pr.enable()
            
        if data_cols is not None and data_units is not None and data_types is not None:
            # @todo check if all lengths are the same
            my_cols =  ["name", "time"]
            my_unit =  ["", "ms"]
            my_type =  [str, float]
            self.table = Table(meta={"name": f"Dysh Benchmark {benchname}"},
                               names=my_cols+data_cols, 
                               units=my_unit+data_units,
                               dtype=my_type+data_types)
            # prepare dummy data[] for the first row in the table
            # or define a data=None in the constructor
            my_data = []
            for t in data_types:
                if type(t) == type(str):
                    my_data.append("")
                elif type(t) == type(float):
                    my_data.append(0.0)
                else:
                    my_data.append(None)
        else:
            print("Warning: skipping table saving")
            self.table = None
            my_data = None
        
        self.stats = []
        self.stats.append(["start", time.perf_counter_ns(), my_data])

    def tag(self, name, data=None):
        self.stats.append([name, time.perf_counter_ns(), data])

    def close(self):
        self.state = 1

    def report(self, debug=False):
        #assert(self.state == 1)
        print(f"# Dysh Benchmark: {self.benchname}")
        n = len(self.stats)
        if debug:
            print(f"Found {n} entries")
            for i in range(n):
                print(self.stats[i])

        for i in range(1,n):
            dt = (self.stats[i][1] - self.stats[i-1][1])/1e6   # in ms, @todo check units
            if True:
                print(self.stats[i][0],dt)
            if self.table is not None:
                self.table.add_row([self.stats[i][0], dt] + self.stats[i][2])
        if self.table is not None:
            self.table["time"].info.format = "0.1f"
            if self.out is not None:
                if os.path.exists(self.out):
                    if self.append:
                        oldtable = Table.read(self.out, format="ipac")
                        table2 = vstack([oldtable, self.table])
                    elif self.overwrite:
                        table2 = self.table
                    else:
                        raise Exception(f"{self.out} exists. Use -w to overwrite.")
                else:
                    table2 = self.table
                print(f"Overwriting {self.out}")
                table2.write(self.out, format="ipac", overwrite=True)
            else:
                self.table.pprint_all()

        if self.profile:
            self.pr.disable()
            ps = pstats.Stats(self.pr).sort_stats(SortKey.CUMULATIVE)
            #ps = pstats.Stats(pr).sort_stats(SortKey.CUMULATIVE, SortKey.TIME)
            #ps.print_stats(int(args.statslines))
            ps.print_stats(self.statslines)

    def total(self):
        """report total time so far"""
        n = len(self.stats)
        dt =  (self.stats[n-1][1] - self.stats[0][1])/1e6
        return dt
        
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
    #   stop here, too much memory
    #a = np.arange(1e9)
    #dt.tag("arange(1e9)")
    dt.close()
    dt.report(debug=False)
