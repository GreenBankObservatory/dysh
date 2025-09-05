"""
Core functions for FITS/SDFITS
"""

__all__ = [
    "default_sdfits_columns",
    "summary_column_definitions",
]


from collections import namedtuple


def default_sdfits_columns():
    """The default column names for GBT SDFITS.

    Returns
    -------
        colnames : list
        A list of the GBT SDFITS column names
    """
    colnames = [
        "OBJECT",
        "BANDWID",
        "DATE-OBS",
        "DURATION",
        "EXPOSURE",
        "TSYS",
        "TDIM7",
        "TUNIT7",
        "CTYPE1",
        "CRVAL1",
        "CRPIX1",
        "CDELT1",
        "CTYPE2",
        "CRVAL2",
        "CTYPE3",
        "CRVAL3",
        "CRVAL4",
        "OBSERVER",
        "OBSID",
        "SCAN",
        "OBSMODE",
        "FRONTEND",
        "TCAL",
        "VELDEF",
        "VFRAME",
        "RVSYS",
        "OBSFREQ",
        "LST",
        "AZIMUTH",
        "ELEVATIO",
        "TAMBIENT",
        "PRESSURE",
        "HUMIDITY",
        "RESTFREQ",
        "DOPFREQ",
        "FREQRES",
        "EQUINOX",
        "RADESYS",
        "TRGTLONG",
        "TRGTLAT",
        "SAMPLER",
        "FEED",
        "SRFEED",
        "FEEDXOFF",
        "FEEDEOFF",
        "SUBREF_STATE",
        "SIDEBAND",
        "PROCSEQN",
        "PROCSIZE",
        "PROCSCAN",
        "PROCTYPE",
        "LASTON",
        "LASTOFF",
        "TIMESTAMP",
        "QD_XEL",
        "QD_EL",
        "QD_BAD",
        "QD_METHOD",
        "VELOCITY",
        "FOFFREF1",
        "ZEROCHAN",
        "ADCSAMPF",
        "VSPDELT",
        "VSPRVAL",
        "VSPRPIX",
        "SIG",
        "CAL",
        "CALTYPE",
        "TWARM",
        "TCOLD",
        "CALPOSITION",
        "BACKEND",
        "PROJID",
        "TELESCOP",
        "SITELONG",
        "SITELAT",
        "SITEELEV",
        "IFNUM",
        "PLNUM",
        "FDNUM",
        "INT",
        "INTNUM",  # not all SDFITS files have INT, so we always create INTNUM
        "NSAVE",
        # The following are added by the GBTFITSLoad constructor.
        # Arguable whether they should be included or not.
        "HDU",
        "BINTABLE",
        "ROW",
        "PROC",
        "OBSTYPE",
        "SUBOBSMODE",
    ]
    return colnames


def summary_column_definitions():
    """
    Column definitions for the summary function of `~dysh.fits.GBTFITSLoad`.

    Returns
    -------
    col_defs : dict
        Dictionary of column definitions.
        For every column this lists the operation to be performed on the column
        when aggregating it by scan number and project id. The data operation is
        stored in the operation attribute. For a description of the operations see
        `https://pandas.pydata.org/pandas-docs/stable/user_guide/groupby.html#built-in-aggregation-methods`_.
        The type attribute is the data type for the column. The name attribute is
        the name to be shown if `verbose=False`.
        The scale attribute is the factor by which to scale the column.
    """

    col_def = namedtuple(
        "col_def",
        ["operation", "type", "name", "scale"],
        defaults=(None, 1),
    )

    str_col = col_def("first", object)

    col_defs = {
        "OBJECT": col_def("first", object),
        "VELOCITY": col_def("mean", float, scale=1e-3),
        "PROC": col_def("first", object),
        "PROCSEQN": col_def("mean", int),
        "PROCSIZE": col_def("mean", int),
        "RESTFREQ": col_def("mean", float, scale=1e-9),
        "DOPFREQ": col_def("mean", float, scale=1e-9),
        "IFNUM": col_def("nunique", int, name="# IF"),
        "FEED": col_def("nunique", int),
        "AZIMUTH": col_def("mean", float),
        "ELEVATIO": col_def("mean", float, name="ELEVATION"),
        "FDNUM": col_def("nunique", int, name="# FEED"),
        "INTNUM": col_def("nunique", int, name="# INT"),
        "PLNUM": col_def("nunique", int, name="# POL"),
        "SIG": col_def("nunique", object, name="# SIG"),
        "CAL": col_def("nunique", object, name="# CAL"),
        "DATE-OBS": col_def("first", object),
        "FITSINDEX": col_def("mean", int, name="# FITS"),
        "SCAN": col_def("mean", int),
        "BANDWID": col_def("mean", float, name="BW", scale=1e-6),
        "DURATION": col_def("mean", float),
        "EXPOSURE": col_def("mean", float),
        "TSYS": col_def("mean", float),
        "CTYPE1": col_def("first", object),
        "CTYPE2": col_def("first", object),
        "CTYPE3": col_def("first", object),
        "CRVAL1": col_def("mean", int),
        "CRVAL2": col_def("mean", int),
        "CRVAL3": col_def("mean", int),
        "CRVAL4": col_def("mean", int),
        "OBSERVER": str_col,
        "OBSID": str_col,
        "OBSMODE": str_col,
        "FRONTEND": str_col,
        "TCAL": col_def("mean", float),
        "VELDEF": str_col,
        "VFRAME": col_def("mean", float, scale=1e-3),
        "RVSYS": col_def("mean", float, scale=1e-3),
        "OBSFREQ": col_def("mean", float, scale=1e-9),
        "LST": col_def("mean", float),
        "FREQRES": col_def("mean", float, scale=1e-9),
        "EQUINOX": str_col,
        "CALTYPE": str_col,
        "TWARM": col_def("mean", float),
        "TCOLD": col_def("mean", float),
        "TAMBIENT": col_def("mean", float),
        "OBSTYPE": str_col,
        "SUBOBSMODE": str_col,
        "PROJID": str_col,
        "SUBREF_STATE": col_def("nunique", int, name="# SUBREF"),
        "BINTABLE": col_def("mean", int),
    }

    return col_defs
