"""
Core functions for FITS/SDFITS
"""

__all__ = [
    "ColDef",
    "default_sdfits_columns",
    "summary_column_definitions",
]


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


class ColDef:
    """
    Class to hold column definitions.

    For every column this lists the operation to be performed on the column
    when aggregating it by scan number and project id. The data operation is
    stored in the operation attribute. For a description of the operations see
    `https://pandas.pydata.org/pandas-docs/stable/user_guide/groupby.html#built-in-aggregation-methods`_.
    The type attribute is the data type for the column. The name attribute is
    the name to be shown when using `summary(verbose=False)`.
    The scale attribute is the factor by which to scale the column.
    The post attribute is any post operations to apply to the column after being aggregated.
    """

    def __init__(self, operation, type, name=None, scale=1, post=None):
        self.operation = operation
        self.type = type
        self.name = name
        self.scale = scale
        self.post = post


def summary_column_definitions():
    """
    Column definitions for the summary function of `~dysh.fits.GBTFITSLoad`.

    Returns
    -------
    col_defs : dict
        Dictionary of column definitions (`~dysh.fits.core.ColDef`).
        See `~dysh.fits.core.ColDef` for details.
    """

    str_col = ColDef("first", object)

    # Use 4 significant figures for Az and El
    # since the servo system accuracy is roughly 1''~0.0003 deg

    col_defs = {
        "OBJECT": ColDef("first", object),
        "VELOCITY": ColDef("mean", float, scale=1e-3, post=lambda x: round(x, 3)),
        "PROC": ColDef("first", object),
        "PROCSEQN": ColDef("mean", int),
        "PROCSIZE": ColDef("mean", int),
        "RESTFREQ": ColDef("mean", float, scale=1e-9),
        "DOPFREQ": ColDef("mean", float, scale=1e-9),
        "IFNUM": ColDef("nunique", int, name="# IF"),
        "FEED": ColDef("nunique", int),
        "AZIMUTH": ColDef("mean", float, post=lambda x: round(x, 4)),
        "ELEVATIO": ColDef("mean", float, name="ELEVATION", post=lambda x: round(x, 4)),
        "FDNUM": ColDef("nunique", int, name="# FEED"),
        "INTNUM": ColDef("nunique", int, name="# INT"),
        "PLNUM": ColDef("nunique", int, name="# POL"),
        "SIG": ColDef(lambda x: set(x), object, name="SIG", post=lambda x: "/".join(x)),
        "CAL": ColDef(lambda x: set(x), object, name="CAL", post=lambda x: "/".join(x)),
        "DATE-OBS": ColDef("first", object),
        "FITSINDEX": ColDef("mean", int, name="# FITS"),
        "SCAN": ColDef("mean", int),
        "BANDWID": ColDef("mean", float, name="BW", scale=1e-6),
        "DURATION": ColDef("mean", float),
        "EXPOSURE": ColDef("mean", float),
        "TSYS": ColDef("mean", float),
        "CTYPE1": ColDef("first", object),
        "CTYPE2": ColDef("first", object),
        "CTYPE3": ColDef("first", object),
        "CRVAL1": ColDef("mean", int),
        "CRVAL2": ColDef("mean", int),
        "CRVAL3": ColDef("mean", int),
        "CRVAL4": ColDef("mean", int),
        "OBSERVER": str_col,
        "OBSID": str_col,
        "OBSMODE": str_col,
        "FRONTEND": str_col,
        "TCAL": ColDef("mean", float),
        "VELDEF": str_col,
        "VFRAME": ColDef("mean", float, scale=1e-3),
        "RVSYS": ColDef("mean", float, scale=1e-3),
        "OBSFREQ": ColDef("mean", float, scale=1e-9),
        "LST": ColDef("mean", float),
        "FREQRES": ColDef("mean", float, scale=1e-9),
        "EQUINOX": str_col,
        "CALTYPE": str_col,
        "TWARM": ColDef("mean", float),
        "TCOLD": ColDef("mean", float),
        "TAMBIENT": ColDef("mean", float),
        "OBSTYPE": str_col,
        "SUBOBSMODE": str_col,
        "PROJID": str_col,
        "SUBREF_STATE": ColDef("nunique", int, name="# SUBREF"),
        "BINTABLE": ColDef("mean", int),
    }

    return col_defs
