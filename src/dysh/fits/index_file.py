"""
SDFITS .index file support

This module provides functions to read and write SDFITS .index files
which contain metadata from SDFITS files in ASCII format for fast access.
The format is compatible with GBTIDL's index format.

Optional fast parsing is available via the rsdfits package (Rust implementation).
Install rsdfits from ~/repos/rsdfits for ~100x faster index file parsing.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass, fields
from datetime import datetime
from pathlib import Path

import pandas as pd

from dysh.log import logger

# Optional fast Rust-based parser
try:
    from rsdfits import parse_sdfits_index_file as _rsdfits_parse

    _RSDFITS_AVAILABLE = True
except ImportError:
    _rsdfits_parse = None
    _RSDFITS_AVAILABLE = False


def is_rsdfits_available() -> bool:
    """Check if the fast Rust-based index parser (rsdfits) is available.

    Returns
    -------
    bool
        True if rsdfits is installed and importable, False otherwise.
    """
    return _RSDFITS_AVAILABLE


# SDFITS index file section markers and constants
HEADER_SECTION_ID = "[header]"
ROWS_SECTION_ID = "[rows]"
TABLE_INDEX_COLUMN = "#INDEX#"


@dataclass
class IndexMetadata:
    """
    Metadata header for SDFITS .index files.

    Attributes
    ----------
    created : str
        Creation timestamp
    last_modified : str
        Last modification timestamp
    version : str
        Index file format version (typically "1.7")
    observer : str
        Observer name
    backend : str
        Backend name
    tcal_rx_table : str
        Tcal receiver table name
    created_by : str
        Program that created the index ("gbtidl" or "dysh")
    sprotect : int
        Protection flag
    """

    created: str
    last_modified: str
    version: str = "1.7"
    observer: str = "Unknown"
    backend: str = "Unknown"
    tcal_rx_table: str = "unknown"
    created_by: str = "dysh"
    sprotect: int = 1


# Column name mapping between SDFITS index and dysh
# The SDFITS index format uses certain legacy column names that differ from
# the FITS binary table column names
SDFITS_INDEX_TO_DYSH_MAP = {
    "INT": "INTNUM",
    "PROJECT": "PROJID",
    "DATEOBS": "DATE-OBS",
    "EXT": "HDU",  # EXT in .index = HDU in dysh (both are HDU number, 1-based)
    "SOURCE": "OBJECT",  # .index calls it SOURCE, dysh calls it OBJECT
    "PROCEDURE": "PROC",  # .index has full name (too long for FITS convention), dysh uses abbreviation
    "PROCS": "PROCSEQN",  # .index has abbreviation, dysh uses full name
    "ELEVATION": "ELEVATIO",  # .index name is too long for FITS
    "LONGITUDE": "CRVAL2",  # .index LONGITUDE = FITS CRVAL2 (lon-like coordinate)
    "LATITUDE": "CRVAL3",  # .index LATITUDE = FITS CRVAL3 (lat-like coordinate)
    "CENTFREQ": "CRVAL1",  # .index CENTFREQ = FITS CRVAL1 (center frequency)
}

DYSH_TO_SDFITS_INDEX_MAP = {v: k for k, v in SDFITS_INDEX_TO_DYSH_MAP.items()}


# --- Writer format specification matching sparrow3/GBTIDL IndexWriter exactly ---

_HEADER_LINE_LEN = 256  # Header lines padded to 256 chars (sparrow3's headerLineLen)

# Writer column spec: (canonical_name, col_header_fmt, row_value_fmt, skip_spacing)
# Matches sparrow3 IndexWriter.cells exactly.
# - col_header_fmt: formats the column name for the header row (may truncate, e.g. EXTENSION → EXT)
# - row_value_fmt: formats the data value for each row
# - skip_spacing: if True, no trailing space in header; in data rows, space added only when int val < 1e6
_WRITER_SPEC = [
    ("#INDEX#", "%7.7s", "%6d", True),
    ("PROJECT", "%16.16s", "%16.16s", False),
    ("FILE", "%64.64s", "%64.64s", False),
    ("EXTENSION", "%3.3s", "%3d", False),
    ("ROW", "%7.7s", "%6d", True),
    ("SOURCE", "%32.32s", "%32.32s", False),
    ("PROCEDURE", "%9.9s", "%9.9s", False),
    ("OBSID", "%32.32s", "%32.32s", False),
    ("E2ESCAN", "%5.5s", "%5d", False),
    ("PROCSEQN", "%5.5s", "%5d", False),
    ("SCAN", "%10.10s", "%10d", False),
    ("POLARIZATION", "%3.3s", "%3.3s", False),
    ("PLNUM", "%5.5s", "%5d", False),
    ("IFNUM", "%5.5s", "%5d", False),
    ("FEED", "%5.5s", "%5d", False),
    ("FDNUM", "%5.5s", "%5d", False),
    ("INT", "%10.10s", "%10d", False),
    ("NUMCHN", "%6.6s", "%6d", False),
    ("SIG", "%3.3s", "%3.3s", False),
    ("CAL", "%3.3s", "%3.3s", False),
    ("SAMPLER", "%12.12s", "%12.12s", False),
    ("AZIMUTH", "%16.16s", "%16.9e", False),
    ("ELEVATION", "%16.16s", "%16.9e", False),
    ("LONGITUDE", "%16.16s", "%16.9e", False),
    ("LATITUDE", "%16.16s", "%16.9e", False),
    ("TRGTLONG", "%16.16s", "%16.9e", False),
    ("TRGTLAT", "%16.16s", "%16.9e", False),
    ("SUBREF", "%3.3s", "%3d", False),
    ("LST", "%16.16s", "%16.9e", False),
    ("CENTFREQ", "%16.16s", "%16.9e", False),
    ("RESTFREQ", "%16.16s", "%16.9e", False),
    ("VELOCITY", "%16.16s", "%16.9e", False),
    ("FREQINT", "%16.16s", "%16.9e", False),
    ("FREQRES", "%16.16s", "%16.9e", False),
    ("DATEOBS", "%22.22s", "%22.22s", False),
    ("TIMESTAMP", "%22.22s", "%22.22s", False),
    ("BANDWIDTH", "%16.16s", "%16.9e", False),
    ("EXPOSURE", "%16.16s", "%16.9e", False),
    ("TSYS", "%16.16s", "%16.9e", False),
    ("NSAVE", "%10.10s", "%10d", False),
    ("PROCSCAN", "%16.16s", "%16.16s", False),
    ("PROCTYPE", "%16.16s", "%16.16s", False),
    ("WCALPOS", "%16.16s", "%16.16s", False),
]

# Column names from the writer spec (#INDEX# → INDEX for DataFrame compatibility)
_WRITER_COLUMNS = [name if name != "#INDEX#" else "INDEX" for name, _, _, _ in _WRITER_SPEC]

# Stokes/polarization code mapping (FITS CRVAL4 integer → string)
_POLARIZATION_MAP = {
    1: "I",
    2: "Q",
    3: "U",
    4: "V",
    -1: "RR",
    -2: "LL",
    -3: "RL",
    -4: "LR",
    -5: "XX",
    -6: "YY",
    -7: "XY",
    -8: "YX",
}


def get_index_path(fits_path: str | Path) -> Path:
    """
    Generate .index filename from FITS path.

    Parameters
    ----------
    fits_path : str or Path
        Path to FITS file

    Returns
    -------
    Path
        Path to corresponding .index file
    """
    fits_path = Path(fits_path)
    # .index file replaces .fits extension with .index
    # e.g., "file.fits" -> "file.index", "dir.vegas.A.fits" -> "dir.vegas.A.index"
    return fits_path.parent / f"{fits_path.stem}.index"


def parse_table_header(line: str) -> dict[str, list[str] | list[int]]:
    """Parse the header row from the [rows] table in SDFITS index file

    The header row is a space-delimited list of column names. Column widths are fixed per-file, but not
    across all files, so we have to parse it each time to determine column widths.

    The header row will start with TABLE_INDEX_COLUMN
    """

    cols: list[str] = []
    starts: list[int] = []
    ends: list[int] = []
    prev_end = None
    start = 0
    for m in re.finditer(r"\S+", line):
        # Everything other than spaces is a column name
        name = m.group()

        if prev_end is not None:
            start = prev_end
        end = m.end()
        prev_end = end
        # For some strange reason, the index column is one character longer than its contents. So, the end index
        # that we use for the table data needs to be reduced by one to avoid having a single space character
        # at the right of each index value
        # This isn't strictly necessary since pandas strips whitespace anyway though
        if name == TABLE_INDEX_COLUMN:
            end -= 1
        cols.append(name)
        starts.append(start)
        ends.append(end)
    return {"cols": cols, "starts": starts, "ends": ends}


def parse_sdfits_index_file(path: Path, use_rust: bool | None = None) -> pd.DataFrame:
    """Given a path to an SDFITS index file, parse it into a DataFrame

    Returns DataFrame with column names mapped from SDFITS index convention to dysh convention
    (e.g., INT -> INTNUM, PROJECT -> PROJID, etc.)

    This function is used internally by sdfitsload.py which expects dysh column names.
    For raw SDFITS index format, use read_index() instead.

    Parameters
    ----------
    path : Path
        Path to the SDFITS .index file
    use_rust : bool or None, optional
        Whether to use the fast Rust parser (rsdfits). If None (default), uses Rust
        if available, otherwise falls back to Python. Set to False to force Python
        parser even when rsdfits is installed.

    Returns
    -------
    pd.DataFrame
        DataFrame with dysh column names and header metadata in df.attrs
    """
    # Determine whether to use Rust parser
    if use_rust is None:
        use_rust = _RSDFITS_AVAILABLE
    elif use_rust and not _RSDFITS_AVAILABLE:
        raise ImportError("rsdfits package not available. Install from ~/repos/rsdfits or set use_rust=False")

    if use_rust:
        logger.debug(f"Parsing {path} with rsdfits parser (fast)")
        return _parse_with_rsdfits(path)
    else:
        logger.debug(f"Parsing {path} with rsdfits parser (slow)")
        return _parse_with_python(path)


def _parse_with_rsdfits(path: Path) -> pd.DataFrame:
    """Parse index file using fast Rust implementation."""
    df = _rsdfits_parse(path)

    # Rename #INDEX# to INDEX if present (rsdfits may or may not do this)
    if "#INDEX#" in df.columns:
        df = df.rename(columns={"#INDEX#": "INDEX"})

    # Apply column name mapping from SDFITS index convention to dysh convention
    df = df.rename(columns=SDFITS_INDEX_TO_DYSH_MAP)

    return df


def _parse_with_python(path: Path) -> pd.DataFrame:
    """Parse index file using pure Python implementation."""
    header = {}
    with open(path) as file:
        header_section_id = next(file)
        if header_section_id != f"{HEADER_SECTION_ID}\n":
            raise ValueError("Invalid SDFITS index file")

        for line in file:
            if line != f"{ROWS_SECTION_ID}\n":
                try:
                    key, value = line.split("=", maxsplit=1)
                except ValueError:
                    logger.debug(f"index file {line} not useful")
                else:
                    header[key.strip()] = value.strip()
            else:
                break

        # This should be the row that contains the header (column names) for the index table
        table_header_raw = next(file)
        if not table_header_raw.startswith(TABLE_INDEX_COLUMN):
            raise ValueError(f"Expected next line of file to be {ROWS_SECTION_ID} table index; got {table_header_raw}")

        table_header = parse_table_header(table_header_raw)
        df = pd.read_fwf(
            file,
            colspecs=tuple(zip(table_header["starts"], table_header["ends"], strict=True)),
            names=table_header["cols"],
            engine="c",
        )
        df.attrs.update(header)

    # Rename #INDEX# to INDEX first
    if "#INDEX#" in df.columns:
        df = df.rename(columns={"#INDEX#": "INDEX"})

    # Apply column name mapping from SDFITS index convention to dysh convention
    df = df.rename(columns=SDFITS_INDEX_TO_DYSH_MAP)

    return df


def read_index(index_path: str | Path) -> tuple[IndexMetadata, pd.DataFrame]:
    """
    Read an SDFITS .index file.

    Parameters
    ----------
    index_path : str or Path
        Path to .index file

    Returns
    -------
    metadata : IndexMetadata
        Index file metadata from header section
    index_df : pd.DataFrame
        Index data as DataFrame (with SDFITS index column names, NOT converted to dysh)

    Notes
    -----
    SDFITS .index files have the following structure:
    - [header] section with key=value pairs
    - [rows] section with column headers and data
    - Column names are kept in SDFITS index format (INT, PROJECT, etc.)
    - Use convert_sdfits_index_to_dysh() to convert column names if needed
    """
    index_path = Path(index_path)

    if not index_path.exists():
        raise FileNotFoundError(f"Index file not found: {index_path}")

    header = {}
    with open(index_path) as file:
        header_section_id = next(file)
        if header_section_id.strip() != HEADER_SECTION_ID:
            raise ValueError("Invalid SDFITS index file - missing [header] section")

        for line in file:
            line_stripped = line.strip()
            if line_stripped == ROWS_SECTION_ID:
                break
            if "=" in line:
                key, value = line.split("=", maxsplit=1)
                header[key.strip()] = value.strip()

        # Create metadata object
        metadata = IndexMetadata(
            created=header.get("created", ""),
            last_modified=header.get("last_modified", ""),
            version=header.get("version", "1.7"),
            observer=header.get("observer", "Unknown"),
            backend=header.get("backend", "Unknown"),
            tcal_rx_table=header.get("tcal_rx_table", "unknown"),
            created_by=header.get("created_by", "gbtidl"),
            sprotect=int(header.get("sprotect", 1)),
        )

        # This should be the row that contains the header (column names) for the index table
        table_header_raw = next(file)
        if not table_header_raw.startswith(TABLE_INDEX_COLUMN):
            raise ValueError(f"Expected table header starting with {TABLE_INDEX_COLUMN}, got {table_header_raw}")

        table_header = parse_table_header(table_header_raw)
        df = pd.read_fwf(
            file,
            colspecs=tuple(zip(table_header["starts"], table_header["ends"], strict=True)),
            names=table_header["cols"],
            engine="c",
        )

    # Rename #INDEX# to INDEX
    if "#INDEX#" in df.columns:
        df = df.rename(columns={"#INDEX#": "INDEX"})

    # Convert data types
    _convert_index_datatypes(df)

    # Return DataFrame with SDFITS column names (NOT converted to dysh)
    return metadata, df


def _convert_index_datatypes(df: pd.DataFrame) -> None:
    """Convert string columns to appropriate data types in-place.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to convert (modified in-place)
    """
    # Handle boolean columns (T/F values)
    bool_cols = ["SIMPLE", "EXTEND", "SIG", "CAL"]
    for col in bool_cols:
        if col in df.columns:
            # Convert T/F to boolean
            df[col] = df[col].map({"T": True, "F": False, "": None})

    # Integer columns
    int_cols = [
        "INDEX",
        "EXT",
        "ROW",
        "E2ESC",
        "PROCS",
        "SCAN",
        "PLNUM",
        "IFNUM",
        "INT",
        "INTNUM",
        "SUB",
        "NSAVE",
        "FDNUM",
        "HDU",
        "BINTABLE",
        "FEED",
        "NUMCHN",
    ]
    for col in int_cols:
        if col in df.columns and df[col].dtype != bool:
            try:
                df[col] = pd.to_numeric(df[col], errors="coerce").astype("Int64")
            except (ValueError, TypeError):
                pass  # Keep as string if conversion fails

    # Float columns
    float_cols = [
        "AZIMUTH",
        "ELEVATIO",
        "ELEVATION",
        "TRGTLONG",
        "TRGTLAT",
        "LST",
        "CENTFREQ",
        "RESTFREQ",
        "VELOCITY",
        "FREQINT",
        "FREQRES",
        "BANDWIDTH",
        "EXPOSURE",
        "TSYS",
        "BANDWID",
        "DURATION",
        "CRVAL1",
        "CRPIX1",
        "CDELT1",
        "LONGITUDE",
        "LATITUDE",
    ]
    for col in float_cols:
        if col in df.columns:
            try:
                df[col] = pd.to_numeric(df[col], errors="coerce")
            except (ValueError, TypeError):
                pass  # Keep as string if conversion fails


def convert_sdfits_index_to_dysh(df: pd.DataFrame) -> pd.DataFrame:
    """Convert column names from SDFITS index format to dysh format.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with SDFITS index column names

    Returns
    -------
    pd.DataFrame
        DataFrame with dysh column names
    """
    return df.rename(columns=SDFITS_INDEX_TO_DYSH_MAP)


def convert_dysh_to_sdfits_index(df: pd.DataFrame) -> pd.DataFrame:
    """Convert column names from dysh format to SDFITS index format.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with dysh column names

    Returns
    -------
    pd.DataFrame
        DataFrame with SDFITS index column names, with dysh-only columns removed
    """
    # Only rename if target column doesn't already exist (to avoid duplicates)
    rename_dict = {}
    for dysh_col, sdfits_col in DYSH_TO_SDFITS_INDEX_MAP.items():
        if dysh_col in df.columns and sdfits_col not in df.columns:
            rename_dict[dysh_col] = sdfits_col

    df = df.rename(columns=rename_dict)

    # Remove dysh-specific columns that don't exist in SDFITS index
    dysh_only_cols = ["HDU", "BINTABLE", "FITSINDEX"]
    cols_to_drop = [col for col in dysh_only_cols if col in df.columns]
    if cols_to_drop:
        df = df.drop(columns=cols_to_drop)

    return df


def count_index_rows(index_path: Path) -> int:
    """Count the number of data rows in an index file.

    Parameters
    ----------
    index_path : Path
        Path to .index file

    Returns
    -------
    int
        Number of data rows (excluding header and table header)
    """
    with open(index_path) as f:
        in_rows_section = False
        count = 0
        for line in f:
            if line.strip() == ROWS_SECTION_ID:
                in_rows_section = True
                next(f)  # Skip the column header line
                continue
            if in_rows_section:
                count += 1
        return count


def read_index_incremental(index_path: Path, start_row: int = 0, end_row: int | None = None) -> pd.DataFrame:
    """Read a subset of rows from an index file.

    Parameters
    ----------
    index_path : Path
        Path to .index file
    start_row : int, optional
        First row to read (0-based), by default 0
    end_row : int, optional
        Last row to read (exclusive), by default None (read all)

    Returns
    -------
    pd.DataFrame
        DataFrame with the requested rows (column names NOT converted to dysh format)
    """
    with open(index_path) as f:
        # Skip to [rows] section
        for line in f:
            if line.strip() == ROWS_SECTION_ID:
                break

        # Parse table header
        table_header_raw = next(f)
        if not table_header_raw.startswith(TABLE_INDEX_COLUMN):
            raise ValueError(f"Expected table header, got: {table_header_raw}")

        table_header = parse_table_header(table_header_raw)

        # Read only the requested rows
        lines_to_read = []
        for i, line in enumerate(f):
            if i >= start_row:
                if end_row is not None and i >= end_row:
                    break
                lines_to_read.append(line)

        if not lines_to_read:
            # Return empty DataFrame with correct columns
            return pd.DataFrame(columns=table_header["cols"])

        # Parse the selected lines
        from io import StringIO

        df = pd.read_fwf(
            StringIO("".join(lines_to_read)),
            colspecs=tuple(zip(table_header["starts"], table_header["ends"], strict=True)),
            names=table_header["cols"],
            engine="c",
        )

    return df


def _generate_rows_header() -> str:
    """Generate the column header row matching sparrow3's IndexWriter format.

    Column names are formatted using the col_header_fmt from _WRITER_SPEC,
    which truncates long names (e.g., EXTENSION → EXT, POLARIZATION → POL).
    Skip-spacing columns (#INDEX#, ROW) get no trailing space in the header.
    """
    header = ""
    last_name = _WRITER_SPEC[-1][0]
    for name, col_fmt, _row_fmt, skip_spacing in _WRITER_SPEC:
        header += col_fmt % name
        if name != last_name and not skip_spacing:
            header += " "
    return header


def _get_polarization(crval4) -> str:
    """Translate CRVAL4 Stokes parameter integer to polarization string.

    Matches sparrow3's IndexWriter.getPolarization().
    """
    try:
        return _POLARIZATION_MAP.get(int(crval4), "??")
    except (ValueError, TypeError):
        return "??"


def _get_center_frequency(crval1, crpix1, cdelt1, numchn) -> float:
    """Compute center frequency from WCS parameters.

    Matches sparrow3's IndexWriter.getCenterFrequency().
    """
    center_chan = (float(numchn) / 2.0) - 0.5
    return ((center_chan - float(crpix1)) * float(cdelt1)) + float(crval1)


def _get_procedure(obsmode) -> str:
    """Extract procedure name from OBSMODE string.

    Matches sparrow3's IndexWriter.getProcedure().
    """
    if pd.isna(obsmode) or not obsmode:
        return ""
    return str(obsmode).split(":")[0]


def _translate_boolean(val) -> str:
    """Convert a boolean-like value to T/F string.

    Matches sparrow3's IndexWriter.translateBoolean().
    """
    if isinstance(val, str):
        upper = val.strip().upper()
        if upper in ("T", "TRUE"):
            return "T"
        return "F"
    if pd.isna(val):
        return "F"
    return "T" if val else "F"


def _prepare_for_writing(df: pd.DataFrame) -> pd.DataFrame:
    """Convert a dysh/FITS DataFrame to canonical SDFITS index column names for writing.

    Computes derived columns (CENTFREQ, POLARIZATION, PROCEDURE, etc.) and renames
    columns from dysh/FITS names to canonical SDFITS index names matching sparrow3.
    """
    df = df.copy()

    # --- Compute derived columns before renaming ---

    # CENTFREQ: compute from WCS parameters if all are available
    if "CRVAL1" in df.columns and "CRPIX1" in df.columns and "CDELT1" in df.columns and "NUMCHN" in df.columns:
        df["CENTFREQ"] = [
            _get_center_frequency(r["CRVAL1"], r["CRPIX1"], r["CDELT1"], r["NUMCHN"]) for _, r in df.iterrows()
        ]
    elif "CRVAL1" in df.columns and "CENTFREQ" not in df.columns:
        # Fallback: use CRVAL1 as approximate CENTFREQ
        df["CENTFREQ"] = df["CRVAL1"]

    # POLARIZATION: derive from CRVAL4
    if "CRVAL4" in df.columns and "POLARIZATION" not in df.columns and "POL" not in df.columns:
        df["POLARIZATION"] = df["CRVAL4"].apply(_get_polarization)

    # PROCEDURE: extract from OBSMODE
    if "OBSMODE" in df.columns and "PROCEDURE" not in df.columns and "PROC" not in df.columns:
        df["PROCEDURE"] = df["OBSMODE"].apply(_get_procedure)

    # FREQINT: from CDELT1
    if "CDELT1" in df.columns and "FREQINT" not in df.columns:
        df["FREQINT"] = df["CDELT1"]

    # BANDWIDTH: from BANDWID
    if "BANDWID" in df.columns and "BANDWIDTH" not in df.columns:
        df["BANDWIDTH"] = df["BANDWID"]

    # SUBREF: from SUBREF_STATE
    if "SUBREF_STATE" in df.columns and "SUBREF" not in df.columns and "SUB" not in df.columns:
        df["SUBREF"] = df["SUBREF_STATE"]

    # E2ESCAN: default to 0 (sparrow3 always sets to 0)
    if "E2ESCAN" not in df.columns and "E2ESC" not in df.columns:
        df["E2ESCAN"] = 0

    # --- Rename columns from dysh/FITS names to canonical SDFITS index names ---
    rename_map = {
        "INTNUM": "INT",
        "PROJID": "PROJECT",
        "DATE-OBS": "DATEOBS",
        "HDU": "EXTENSION",
        "OBJECT": "SOURCE",
        "PROC": "PROCEDURE",
        "ELEVATIO": "ELEVATION",
        "CRVAL2": "LONGITUDE",
        "CRVAL3": "LATITUDE",
        # Abbreviated → canonical (for DataFrames from old dysh-written index files)
        "EXT": "EXTENSION",
        "POL": "POLARIZATION",
        "E2ESC": "E2ESCAN",
        "PROCS": "PROCSEQN",
        "SUB": "SUBREF",
    }
    actual_renames = {src: dst for src, dst in rename_map.items() if src in df.columns and dst not in df.columns}
    if actual_renames:
        df = df.rename(columns=actual_renames)

    # --- Convert SIG/CAL to T/F strings ---
    for col in ("SIG", "CAL"):
        if col in df.columns:
            df[col] = df[col].apply(_translate_boolean)

    return df


def write_index(index_path: str | Path, metadata: IndexMetadata, df: pd.DataFrame):
    """
    Write an SDFITS-compatible .index file.

    Parameters
    ----------
    index_path : str or Path
        Path to output .index file
    metadata : IndexMetadata
        Index file metadata
    df : pd.DataFrame
        Index data (will be converted to SDFITS index format)

    Notes
    -----
    Creates an SDFITS-compatible ASCII format matching sparrow3/GBTIDL:
    - Header section with 256-character padded lines
    - Data section with fixed-width columns matching sparrow3's IndexWriter
    - Scientific notation for floats (%16.9e)
    - T/F for booleans
    - Column header generated from format specs (matches sparrow3 exactly)
    """
    index_path = Path(index_path)

    # Convert to canonical SDFITS index format
    df_sdfits = _prepare_for_writing(df)

    # Ensure all standard columns exist with appropriate defaults
    for col in _WRITER_COLUMNS:
        if col not in df_sdfits.columns:
            # Find the spec entry for this column
            spec_name = "#INDEX#" if col == "INDEX" else col
            spec = next((s for s in _WRITER_SPEC if s[0] == spec_name), None)
            if spec:
                row_fmt = spec[2]
                if col == "NSAVE":
                    df_sdfits[col] = -1  # sparrow3 defaults NSAVE to -1
                elif "d" in row_fmt:
                    df_sdfits[col] = 0
                elif "e" in row_fmt:
                    df_sdfits[col] = 0.0
                else:
                    df_sdfits[col] = ""

    with open(index_path, "w") as f:
        # Write header section
        f.write("[header]\n")

        # Write each metadata field with 256-char padding (matching sparrow3)
        for field in fields(IndexMetadata):
            key = field.name
            value = getattr(metadata, key)
            line = f"{key} = {value}"
            f.write(line.ljust(_HEADER_LINE_LEN) + "\n")

        # Write [rows] marker
        f.write("[rows]\n")

        # Write column header (generated dynamically to match sparrow3 format)
        f.write(_generate_rows_header() + "\n")

        # Write data rows
        for _, row in df_sdfits.iterrows():
            f.write(_format_sdfits_row(row) + "\n")


def _format_sdfits_row(row: pd.Series) -> str:
    """Format a single row matching sparrow3's IndexWriter format.

    Parameters
    ----------
    row : pd.Series
        Row with canonical SDFITS index column names

    Returns
    -------
    str
        Formatted row string
    """
    result = ""
    last_name = _WRITER_SPEC[-1][0]

    for name, _col_fmt, row_fmt, skip_spacing in _WRITER_SPEC:
        # Handle #INDEX# which is stored as INDEX in DataFrame
        lookup_name = "INDEX" if name == "#INDEX#" else name
        val = row.get(lookup_name)

        # Determine value type from format string and apply defaults
        if "d" in row_fmt:
            # Integer format
            if pd.isna(val) or val is None or val == "":
                val = 0
            else:
                val = int(val)
        elif "e" in row_fmt:
            # Float scientific notation
            if pd.isna(val) or val is None or val == "":
                val = 0.0
            else:
                val = float(val)
        elif pd.isna(val) or val is None:
            # String format, missing value
            val = ""
        else:
            # String format
            val = str(val)

        result += row_fmt % val

        # Spacing logic matching sparrow3's initRowString:
        # - skip_spacing columns (#INDEX#, ROW): add space only when int value < 1e6
        # - other columns: always add space (except the last column)
        if name != last_name:
            if not skip_spacing or (isinstance(val, int) and val < 1000000):
                result += " "

    return result


def validate_index(fits_path: str | Path, index_path: str | Path) -> bool:
    """
    Check if .index file is up-to-date with FITS file.

    Parameters
    ----------
    fits_path : str or Path
        Path to FITS file
    index_path : str or Path
        Path to .index file

    Returns
    -------
    bool
        True if index is valid and up-to-date, False otherwise
    """
    fits_path = Path(fits_path)
    index_path = Path(index_path)

    # Check if index file exists
    if not index_path.exists():
        return False

    # Check if FITS file exists
    if not fits_path.exists():
        return False

    # Compare modification times
    fits_mtime = os.path.getmtime(fits_path)
    index_mtime = os.path.getmtime(index_path)

    # Index is valid if it's not significantly older than FITS
    # Allow small time differences (e.g., filesystem precision, clock skew)
    time_diff = fits_mtime - index_mtime
    if time_diff > 1.0:  # More than 1 second older
        return False

    return True


def create_index_metadata(
    observer: str = "Unknown", backend: str = "Unknown", tcal_rx_table: str = "unknown"
) -> IndexMetadata:
    """
    Create index metadata with current timestamp.

    Parameters
    ----------
    observer : str, optional
        Observer name
    backend : str, optional
        Backend name
    tcal_rx_table : str, optional
        Tcal receiver table name

    Returns
    -------
    IndexMetadata
        Index metadata object with current timestamps
    """
    now = datetime.now().strftime("%a %b %d %H:%M:%S %Y")

    return IndexMetadata(
        created=now,
        last_modified=now,
        version="1.7",
        observer=observer,
        backend=backend,
        tcal_rx_table=tcal_rx_table,
        created_by="dysh",
        sprotect=1,
    )
