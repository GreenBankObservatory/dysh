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
    "PROCEDURE": "PROC",  # .index has full name, dysh uses abbreviation
    "ELEVATION": "ELEVATIO",  # .index has correct spelling, dysh has typo
    "LONGITUDE": "CRVAL2",  # .index LONGITUDE = FITS CRVAL2 (lon-like coordinate)
    "LATITUDE": "CRVAL3",  # .index LATITUDE = FITS CRVAL3 (lat-like coordinate)
    "CENTFREQ": "CRVAL1",  # .index CENTFREQ = FITS CRVAL1 (center frequency)
}

DYSH_TO_SDFITS_INDEX_MAP = {v: k for k, v in SDFITS_INDEX_TO_DYSH_MAP.items()}


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
                    print(f"{line} not useful")
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
    Creates an SDFITS-compatible ASCII format with:
    - Header section with 200-character padded lines
    - Data section with fixed-width columns
    - Scientific notation for floats
    - T/F for booleans
    """
    index_path = Path(index_path)

    # Convert dysh column names to SDFITS index names
    df_sdfits = convert_dysh_to_sdfits_index(df.copy())

    # Ensure all 43 standard SDFITS index columns exist (add missing ones with default values)
    standard_cols = [
        "INDEX",
        "PROJECT",
        "FILE",
        "EXT",
        "ROW",
        "SOURCE",
        "PROCEDURE",
        "OBSID",
        "E2ESC",
        "PROCS",
        "SCAN",
        "POL",
        "PLNUM",
        "IFNUM",
        "FEED",
        "FDNUM",
        "INT",
        "NUMCHN",
        "SIG",
        "CAL",
        "SAMPLER",
        "AZIMUTH",
        "ELEVATION",
        "LONGITUDE",
        "LATITUDE",
        "TRGTLONG",
        "TRGTLAT",
        "SUB",
        "LST",
        "CENTFREQ",
        "RESTFREQ",
        "VELOCITY",
        "FREQINT",
        "FREQRES",
        "DATEOBS",
        "TIMESTAMP",
        "BANDWIDTH",
        "EXPOSURE",
        "TSYS",
        "NSAVE",
        "PROCSCAN",
        "PROCTYPE",
        "WCALPOS",
    ]

    for col in standard_cols:
        if col not in df_sdfits.columns:
            # Add missing column with appropriate default value
            if col in ["SIG", "CAL"]:
                df_sdfits[col] = False
            elif col in [
                "AZIMUTH",
                "ELEVATION",
                "LONGITUDE",
                "LATITUDE",
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
            ]:
                df_sdfits[col] = 0.0
            elif col in [
                "INDEX",
                "EXT",
                "ROW",
                "E2ESC",
                "PROCS",
                "SCAN",
                "PLNUM",
                "IFNUM",
                "FEED",
                "FDNUM",
                "INT",
                "NUMCHN",
                "SUB",
                "NSAVE",
            ]:
                df_sdfits[col] = 0
            else:
                df_sdfits[col] = ""

    # Reorder columns to match standard order
    df_sdfits = df_sdfits[standard_cols]

    with open(index_path, "w") as f:
        # Write header section
        f.write("[header]\n")

        # Write each metadata field with 200-char padding
        for field in fields(IndexMetadata):
            key = field.name
            value = getattr(metadata, key)
            line = f"{key} = {value}"
            # Pad to 200 characters
            line = line.ljust(200)
            f.write(line + "\n")

        # Write [rows] marker
        f.write("[rows]\n")

        # Write column header
        f.write(
            "#INDEX#  PROJECT          FILE                                                            EXT ROW     SOURCE                           PROCEDURE OBSID                            E2ESC PROCS SCAN       POL PLNUM IFNUM FEED  FDNUM INT        NUMCHN SIG CAL SAMPLER      AZIMUTH          ELEVATION        LONGITUDE        LATITUDE         TRGTLONG         TRGTLAT          SUB LST              CENTFREQ         RESTFREQ         VELOCITY         FREQINT          FREQRES          DATEOBS               TIMESTAMP             BANDWIDTH        EXPOSURE         TSYS             NSAVE      PROCSCAN         PROCTYPE         WCALPOS\n"
        )

        # Write data rows
        for _, row in df_sdfits.iterrows():
            line = _format_sdfits_row(row)
            f.write(line + "\n")


def _format_sdfits_row(row: pd.Series) -> str:
    """Format a single row in SDFITS index format.

    Parameters
    ----------
    row : pd.Series
        Row with SDFITS index column names

    Returns
    -------
    str
        Formatted row string
    """
    # SDFITS column format specification
    # Format: (column_name, width, format_type)
    # format_type: 'i'=integer (right-aligned), 'f'=float (scientific), 's'=string (left-aligned), 'b'=boolean (T/F)
    sdfits_spec = [
        ("INDEX", 7, "i"),
        ("PROJECT", 16, "s"),
        ("FILE", 64, "s"),
        ("EXT", 3, "i"),
        ("ROW", 7, "i"),
        ("SOURCE", 32, "s"),
        ("PROCEDURE", 9, "s"),
        ("OBSID", 32, "s"),
        ("E2ESC", 5, "i"),
        ("PROCS", 5, "i"),
        ("SCAN", 10, "i"),
        ("POL", 3, "s"),
        ("PLNUM", 5, "i"),
        ("IFNUM", 5, "i"),
        ("FEED", 5, "i"),
        ("FDNUM", 5, "i"),
        ("INT", 10, "i"),
        ("NUMCHN", 6, "i"),
        ("SIG", 3, "b"),
        ("CAL", 3, "b"),
        ("SAMPLER", 12, "s"),
        ("AZIMUTH", 16, "f"),
        ("ELEVATION", 16, "f"),
        ("LONGITUDE", 16, "f"),
        ("LATITUDE", 16, "f"),
        ("TRGTLONG", 16, "f"),
        ("TRGTLAT", 16, "f"),
        ("SUB", 3, "i"),
        ("LST", 16, "f"),
        ("CENTFREQ", 16, "f"),
        ("RESTFREQ", 16, "f"),
        ("VELOCITY", 16, "f"),
        ("FREQINT", 16, "f"),
        ("FREQRES", 16, "f"),
        ("DATEOBS", 21, "s"),
        ("TIMESTAMP", 21, "s"),
        ("BANDWIDTH", 16, "f"),
        ("EXPOSURE", 16, "f"),
        ("TSYS", 16, "f"),
        ("NSAVE", 10, "i"),
        ("PROCSCAN", 16, "s"),
        ("PROCTYPE", 16, "s"),
        ("WCALPOS", 16, "s"),
    ]

    parts = []
    for col_name, width, fmt in sdfits_spec:
        val = row.get(col_name, "")

        # Format based on type
        if fmt == "i":
            # Integer: right-aligned
            if pd.isna(val) or val == "":
                formatted = str(0).rjust(width)
            else:
                formatted = str(int(val)).rjust(width)
        elif fmt == "f":
            # Float scientific notation
            if pd.isna(val) or val == "":
                formatted = "0.000000000e+00".rjust(width)
            else:
                formatted = f"{float(val):.9e}".rjust(width)
        elif fmt == "b":
            # Boolean: T/F, left-aligned
            if pd.isna(val) or val == "":
                formatted = "F".ljust(width)
            else:
                val_str = "T" if val else "F"
                formatted = val_str.ljust(width)
        elif fmt == "s":
            # String: left-aligned
            if pd.isna(val) or val == "":
                formatted = "".ljust(width)
            else:
                formatted = str(val).ljust(width)
        else:
            formatted = str(val).ljust(width)

        # Truncate if too long
        formatted = formatted[:width]
        parts.append(formatted)
        parts.append(" ")  # Space between columns

    return "".join(parts).rstrip()  # Remove trailing space


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
