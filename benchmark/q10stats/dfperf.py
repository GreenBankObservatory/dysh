#!/usr/bin/env python3
"""
Created on Thu Jul 17 11:26:52 2025

@author: mpound
"""

import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal

from dysh.fits import GBTFITSLoad
from dysh.util import get_project_testdata
from dysh.util.files import dysh_data, valid_dysh_accept
from dysh.util.timers import DTime


def bymerge(df1, df2):
    merged_df = pd.merge(df1, df2, how="outer", indicator=True)
    final_df = merged_df[merged_df["_merge"] != "both"]
    return final_df.empty


def bycompare(df1, df2):
    return (df1.compare(df2)).empty


def byassert(df1, df2):
    # print("BYASSERT")
    try:
        # assert_frame_equal(pd.DataFrame(df1), pd.DataFrame(df2))  # check_names=True)
        assert_frame_equal(df1, df2)
    except AssertionError:
        print("error")
        return False
    # print("EQUAL")
    return True


def bynumpy(df1, df2):
    return np.all(df1.to_numpy(na_value=0) == df2.to_numpy(na_value=0))


data_cols = ["#files", "nchan", "nrow", "EQUAL?"]
data_units = ["", "", "", ""]
data_types = [int, int, int, bool]
sdf_file = get_project_testdata() / "AGBT05B_047_01/AGBT05B_047_01.raw.acs/AGBT05B_047_01.raw.acs.fits"
accept = list(valid_dysh_accept.keys())[9:]
print(f"{accept=}")
num = 1
dt = DTime("df", data_cols=data_cols, data_units=data_units, data_types=data_types)
for k in accept:
    sdf_file = dysh_data(dysh_data="/bigdisk/data/gbt/dysh_data", accept=accept[0])
    sdf = GBTFITSLoad(sdf_file)
    s = sdf.stats()
    df = sdf._selection.copy()
    dt.tag("load " + k, data=[s["nfiles"], s["nchan"], s["nrows"], False])

    for _i in range(num):
        eq = df.equals(sdf._selection)
    dt.tag("equals ", data=[s["nfiles"], s["nchan"], s["nrows"], eq])
    for _i in range(num):
        eq = all(df.eq(sdf._selection, axis=1))

    dt.tag("eq ", data=[s["nfiles"], s["nchan"], s["nrows"], eq])
    for _i in range(num):
        eq = bymerge(df, sdf._selection)
    dt.tag("merge ", data=[s["nfiles"], s["nchan"], s["nrows"], eq])
    for _i in range(num):
        eq = bycompare(df, sdf._selection.copy())
    dt.tag("compare ", data=[s["nfiles"], s["nchan"], s["nrows"], eq])
    for _i in range(num):
        eq = byassert(df, sdf._selection.copy())
    dt.tag("assert ", data=[s["nfiles"], s["nchan"], s["nrows"], eq])
    for _i in range(num):
        eq = bynumpy(df, sdf._selection)
    dt.tag("numpy ", data=[s["nfiles"], s["nchan"], s["nrows"], eq])
dt.report()
