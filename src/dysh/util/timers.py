#! /usr/bin/env python
#
import cProfile
import os
import platform
import pstats
import time
from pstats import SortKey

import numpy as np
from astropy.table import Table, vstack

__all__ = ["DTime"]

__ostype__ = None

# see also ADMIT's util.utils.Dtime()


class DTime:
    r"""This class encapsulated some popular timing/performance tools.

    Parameters
    ----------

    benchname : str
         Identifying name of the benchmark stored in the metadata of the table

    units : str
         Units. Allowed are "ms" (the default), others not implemented yet, if ever.

    data_cols : list
         List of names of the extra columns (in addition to the default name and time) written
         to an Astropy at the report stage of this class.

    data_units : list
         List of units names of the extra columns.

    data_types : list
         List of data types of the extra columns.

    args: dict
         This dictionary controls a number of common variables used in dysh benchmarking.

         out        : output filename (astropy Table). Default is no file is written.
         append     : append to previous output file (astropy Table).
         overwrite  : overwrite a previous output file (astropy Table).
         profile    : run the profiler: Default False
         statslines : number of profiler statistics lines to print. Default 25
         sortkey    : how to sort the profiler statistics, "cumulative" or "time". Default "cumulative" (SortKey.CUMULATIVE).


    Examples
    --------

    >>> dt = DTime()
    >>> dt.tag("test1")
    >>> dt.tag("test2")
    >>> dt.tag("test3")
    >>> dt.report()

    By default it simply builds a delta-time of the time it took between the different tags, as
    labeled by their tag name. If DTime() is supplied a number of data items for extra columns,
    these will be reported, or stored in a table, if out= is supplied.


    """

    def __init__(
        self,
        benchname="generic",
        units="ms",
        active=True,
        data_cols=None,
        data_units=None,
        data_types=None,
        args=None,
    ):
        # out = None, overwrite=False, append=False, profile=False, statslines=25,
        self.active = active
        if not self.active:
            return
        self.benchname = benchname
        self.state = 0
        self.ndata = 0
        self._sortkeys = {"cumulative": SortKey.CUMULATIVE, "time": SortKey.TIME}
        if args is not None:
            self.out = args["out"]  # @todo check the dictionary
            self.append = args["append"]
            self.overwrite = args["overwrite"]
            self.profile = args["profile"]
            self.statslines = int(args["statslines"])
            self.sortkey = self._sortkeys[args["sortkey"]]
        else:
            self.out = None
            self.append = False
            self.overwrite = False
            self.profile = False
            self.statslines = 0
            self.sortkey = self._sortkeys["cumulative"]

        if self.profile:
            self.pr = cProfile.Profile()
            self.pr.enable()

        # standard columns if no data[] are given
        my_cols = ["name", "time", "VmSize", "VmRSS"]
        my_unit = ["", "ms", "MB", "MB"]
        my_type = [str, float, float, float]
        my_data = []

        if data_cols is not None and data_units is not None and data_types is not None:
            self.ndata = len(data_cols)
            # print("PJT ndata",self.ndata)
            if len(data_units) != self.ndata:
                raise ValueError(f"data_units length ({len(data_units)}) must match data_cols length ({self.ndata})")
            if len(data_types) != self.ndata:
                raise ValueError(f"data_types length ({len(data_types)}) must match data_cols length ({self.ndata})")
            my_cols = my_cols + data_cols
            my_unit = my_unit + data_units
            my_type = my_type + data_types

        self.table = Table(
            meta={"name": f"Dysh Benchmark {benchname}"},
            names=my_cols,
            units=my_unit,
            dtype=my_type,
        )
        self.stats = []

        # prepare for the first row in the table
        my_data = []
        for i in range(self.ndata):
            t = data_types[i]
            if type(t) is type(str):
                my_data.append("")
            elif type(t) is type(float):
                my_data.append(0.0)
            else:
                my_data.append(None)

        self.stats.append(["start", time.perf_counter_ns(), 0.0, 0.0, my_data])

    def active(self):
        return self.active

    def tag(self, name, data=None):
        """ """
        if not self.active:
            return
        if data is not None and len(data) != self.ndata:
            raise ValueError(f"data length ({len(data)}) must match expected length ({self.ndata})")
        mem1, mem2 = self._mem()
        self.stats.append([name, time.perf_counter_ns(), mem1, mem2, data])

    def close(self):
        """ """
        # print("PJT close")
        self.state = 1
        if not self.active:
            return

    def enable(self):
        """Enable the profiler"""
        self.pr.enable()

    def disable(self):
        """Disable the profiler"""
        self.pr.disable()

    def _mem(self):
        """Read memory usage info from /proc/pid/status
        Return Virtual and Resident memory size in MBytes.

        @todo   add implementation for Mac (see ADMIT)
        """
        global __ostype__

        if __ostype__ is None:
            # __ostype__ = os.uname()[0].lower()
            __ostype__ = platform.uname()[0].lower()
            print("Found ostype=", __ostype__)

        scale = {"MB": 1024.0}
        lines = []

        try:
            if __ostype__ == "linux":
                proc_status = f"/proc/{os.getpid()}/status"  # linux only
                # open pseudo file  /proc/<pid>/status
                t = open(proc_status)
                # get value from line e.g. 'VmRSS:  9999  kB\n'
                for it in t.readlines():
                    if "VmSize" in it or "VmRSS" in it:
                        lines.append(it)
                t.close()
            else:
                print("no get_mem yet")
                return np.array([0, 0])
        except Exception:
            print("error get_mem")
            return np.array([])

        mem = {}
        if __ostype__ != "darwin":
            for line in lines:
                words = line.strip().split()
                key = words[0][:-1]
                scaled = float(words[1]) / scale["MB"]
                mem[key] = scaled

        return np.array([mem["VmSize"], mem["VmRSS"]])

    def report(self, debug=False):
        """ """
        if not self.active:
            return
            #        assert(self.state == 1)
        # print("PJT state",self.state)

        n = len(self.stats)
        if debug:
            print(f"Found {n} entries")
            for i in range(n):
                print(self.stats[i])

        for i in range(1, n):
            dt = (self.stats[i][1] - self.stats[i - 1][1]) / 1e6  # in ms, @todo check units
            if debug:
                print(self.stats[i][0], dt)
            if self.table is not None:
                mem1 = self.stats[i][2]
                mem2 = self.stats[i][3]
                if self.ndata == 0:
                    self.table.add_row([self.stats[i][0], dt, mem1, mem2])
                else:
                    self.table.add_row([self.stats[i][0], dt, mem1, mem2] + self.stats[i][4])
        if self.table is not None:
            self.table["time"].info.format = "0.1f"
            self.table["VmSize"].info.format = "0.1f"
            self.table["VmRSS"].info.format = "0.1f"
            if self.out is not None:
                if os.path.exists(self.out):
                    if self.append:
                        oldtable = Table.read(self.out, format="ascii.ecsv")
                        table2 = vstack([oldtable, self.table])
                    elif self.overwrite:
                        table2 = self.table
                    else:
                        raise Exception(f"{self.out} exists. Use -w to overwrite.")
                else:
                    table2 = self.table
                print(f"Overwriting {self.out}")
                table2.write(self.out, format="ascii.ecsv", overwrite=True)
            else:
                print(f"# Dysh Benchmark: {self.benchname}")
                self.table.pprint_all()

        if self.profile:
            self.pr.disable()
            ps = pstats.Stats(self.pr).sort_stats(self.sortkey)
            # ps = pstats.Stats(pr).sort_stats(SortKey.CUMULATIVE, SortKey.TIME)
            # ps.print_stats(int(args.statslines))
            ps.print_stats(self.statslines)

    def total(self):
        """
        report total CPU time so far
        """
        if not self.active:
            return 0
        n = len(self.stats)
        dt = (self.stats[n - 1][1] - self.stats[0][1]) / 1e6
        return dt


if __name__ == "__main__":
    #  Also compare the output of this with /usr/bin/time on the executable
    #  to find any overhead of the timer we didn't account. Optionally run
    #  a benchmark once and twice and use the difference.
    dt = DTime(active=True)
    dt.tag("nothing    ")
    dt.tag("test0      ")
    dt.tag("test1      ")
    dt.tag("test2      ")
    dt.tag("test3      ")
    dt.tag("test4      ")
    for k in range(3, 9):
        a = np.arange(10**k)
        dt.tag(f"arange(1e{k})")
    dt.close()
    dt.report(debug=False)
    print("Final total:", dt.total())
