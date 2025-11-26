"""
Online Mode Simulator for SDFITS Sessions

This module provides utilities to simulate SDFITS session creation in real-time,
writing data scan-by-scan to allow testing of online monitoring features.

The simulator reads from existing SDFITS sessions and recreates them incrementally,
either at a fixed rate or respecting the original observation timing.
"""

from __future__ import annotations

import logging
import os
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path
from typing import Literal

import fitsio
import numpy as np
import pandas as pd

from dysh.fits.index_file import parse_sdfits_index_file
from dysh.fits.sdfitsload import SDFITSLoad

logger = logging.getLogger(__name__)


def simulate_session(
    source_dir: str | Path,
    output_dir: str | Path,
    mode: Literal["fixed", "realtime"] = "fixed",
    interval: float = 1.0,
    speedup: float = 1.0,
    create_index: bool = True,
    parallel: bool = True,
) -> None:
    """
    Simulate SDFITS session creation by writing scans incrementally.

    Reads from an existing SDFITS session directory and recreates it scan-by-scan
    in a new location, using minimal memory by loading only one scan at a time.

    Parameters
    ----------
    source_dir : str or Path
        Path to source SDFITS session directory (e.g., "PROJECT.raw.vegas/")
    output_dir : str or Path
        Path where simulated session will be created
    mode : {'fixed', 'realtime'}, default='fixed'
        Timing mode:
        - 'fixed': Write one scan every `interval` seconds
        - 'realtime': Use actual time deltas from observation, scaled by `speedup`
    interval : float, default=1.0
        For 'fixed' mode: seconds between scan writes
    speedup : float, default=1.0
        For 'realtime' mode: factor to speed up/slow down replay (1.0 = original speed)
    create_index : bool, default=True
        If True, also create/update .index files incrementally
    parallel : bool, default=True
        If True, write multiple FITS files (A, B, C, D) in parallel

    Examples
    --------
    >>> # Simulate at 1 scan per second
    >>> simulate_session("AGBT21A.raw.vegas", "/tmp/sim_output", mode="fixed", interval=1.0)

    >>> # Simulate at 10x original observation speed
    >>> simulate_session("AGBT21A.raw.vegas", "/tmp/sim_output", mode="realtime", speedup=10.0)

    Notes
    -----
    - Only loads one scan's worth of rows at a time (memory efficient)
    - Creates sdfitsStatus.txt file that updates as scans are written
    - Output directory will be created if it doesn't exist
    - Will fail loudly on any errors (no recovery attempts)
    """
    source_dir = Path(source_dir).resolve()
    output_dir = Path(output_dir).resolve()

    if not source_dir.is_dir():
        raise ValueError(f"Source directory does not exist: {source_dir}")

    # Check if source_dir is a project directory (contains .raw.* subdirs)
    # or a session directory (contains .fits files directly)
    fits_files = sorted(source_dir.glob("*.fits"))
    if not fits_files:
        # Look for .raw.vegas or .raw.acs subdirectories
        raw_dirs = sorted(source_dir.glob("*.raw.*"))
        if raw_dirs:
            # Use the first one found (typically only one per project)
            source_dir = raw_dirs[0]
            logger.info(f"Found session directory: {source_dir.name}")
            fits_files = sorted(source_dir.glob("*.fits"))

    if not fits_files:
        raise ValueError(f"No FITS files found in {source_dir}")

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # Check disk space before starting
    _check_disk_space(fits_files, output_dir)

    logger.info(f"Simulating session from {source_dir} to {output_dir}")
    logger.info(f"Mode: {mode}, interval={interval}s, speedup={speedup}x, create_index={create_index}")

    logger.info(f"Found {len(fits_files)} FITS files: {[f.name for f in fits_files]}")

    # Extract project name and backend from directory name
    # e.g., "AGBT21A_501_11.raw.vegas" -> project="AGBT21A_501_11", backend="vegas"
    dir_name = source_dir.name
    if ".raw." in dir_name:
        project_name, backend = dir_name.split(".raw.")
    else:
        raise ValueError(f"Cannot parse project/backend from directory name: {dir_name}")

    # Prepare sdfitsStatus.txt path
    status_file = output_dir / "sdfitsStatus.txt"
    _initialize_status_file(status_file)

    # Process files in parallel or sequentially
    if parallel and len(fits_files) > 1:
        logger.info("Writing FITS files in parallel")
        with ThreadPoolExecutor(max_workers=len(fits_files)) as executor:
            futures = {
                executor.submit(
                    _simulate_single_file,
                    fits_file,
                    output_dir,
                    mode,
                    interval,
                    speedup,
                    create_index,
                    project_name,
                    backend,
                    status_file,
                ): fits_file
                for fits_file in fits_files
            }

            for future in as_completed(futures):
                fits_file = futures[future]
                try:
                    future.result()
                    logger.info(f"✓ Completed: {fits_file.name}")
                except Exception as e:
                    logger.error(f"✗ Failed: {fits_file.name}")
                    raise RuntimeError(f"Error simulating {fits_file.name}: {e}") from e
    else:
        logger.info("Writing FITS files sequentially")
        for fits_file in fits_files:
            _simulate_single_file(
                fits_file,
                output_dir,
                mode,
                interval,
                speedup,
                create_index,
                project_name,
                backend,
                status_file,
            )
            logger.info(f"✓ Completed: {fits_file.name}")

    logger.info(f"✓ Session simulation complete: {output_dir}")


def _simulate_single_file(
    source_file: Path,
    output_dir: Path,
    mode: str,
    interval: float,
    speedup: float,
    create_index: bool,
    project_name: str,
    backend: str,
    status_file: Path,
) -> None:
    """Simulate a single FITS file by writing it scan-by-scan."""
    output_file = output_dir / source_file.name
    output_index = output_dir / f"{source_file.stem}.index"

    logger.info(f"Processing {source_file.name} -> {output_file}")

    # Read scan metadata from index file (or generate it)
    scan_metadata = _read_scan_metadata(source_file)

    if scan_metadata.empty:
        logger.warning(f"No scans found in {source_file.name}, skipping")
        return

    unique_scans = sorted(scan_metadata["SCAN"].unique())
    total_rows = len(scan_metadata)
    logger.info(f"  Found {len(unique_scans)} scans, {total_rows} total rows")

    # Open source FITS file for reading
    source_fits = fitsio.FITS(str(source_file), "r")

    # Get headers from source for use when writing
    bintable_header = source_fits[1].read_header()

    # Create index file if requested
    if create_index:
        _initialize_index_file(output_index, scan_metadata)

    # Write each scan incrementally
    scan_times = scan_metadata.groupby("SCAN")["UTC"].first() if mode == "realtime" else None
    prev_time = None
    rows_written = 0

    try:
        for i, scan_num in enumerate(unique_scans):
            # Get row indices for this scan
            scan_rows = scan_metadata[scan_metadata["SCAN"] == scan_num]["ROW"].values
            rows_written += len(scan_rows)

            logger.info(
                f"  [{i+1}/{len(unique_scans)}] Scan {scan_num}: "
                f"{len(scan_rows)} rows ({rows_written}/{total_rows} total) -> {output_file.name}"
            )

            # Load only this scan's data from source
            scan_data = source_fits[1].read(rows=list(scan_rows))

            # Write to output FITS
            if i == 0:
                # First scan: create file and write data with header
                with fitsio.FITS(str(output_file), "rw", clobber=True) as output_fits:
                    output_fits.write(scan_data, header=bintable_header)
            else:
                # Subsequent scans: append to existing extension
                with fitsio.FITS(str(output_file), "rw") as output_fits:
                    output_fits[-1].append(scan_data)

            # Update index file if requested
            if create_index:
                _append_to_index_file(output_index, scan_metadata[scan_metadata["SCAN"] == scan_num])

            # Update sdfitsStatus.txt
            _update_status_file(status_file, backend, project_name, scan_num, output_file, output_index)

            # Calculate sleep time
            if mode == "realtime":
                current_time = scan_times.loc[scan_num]
                if prev_time is not None:
                    delta = (current_time - prev_time).total_seconds() / speedup
                    if delta > 0:
                        time.sleep(delta)
                prev_time = current_time
            else:  # fixed mode
                time.sleep(interval)

    except KeyboardInterrupt:
        logger.info(f"  Interrupted after {rows_written}/{total_rows} rows")
        raise
    finally:
        source_fits.close()

    logger.info(f"  ✓ Finished {source_file.name}: {rows_written} rows written")


def _read_scan_metadata(fits_file: Path) -> pd.DataFrame:
    """
    Read scan metadata from .index file, or generate it from FITS if missing.

    Returns DataFrame with columns: ROW, SCAN, UTC, and other metadata needed
    for determining scan boundaries and timing.
    """
    index_file = fits_file.parent / f"{fits_file.stem}.index"

    if index_file.exists():
        logger.debug(f"Reading existing index file: {index_file.name}")
        df = parse_sdfits_index_file(index_file)
    else:
        logger.debug(f"No index file found, generating from FITS: {fits_file.name}")
        # Load FITS and generate index
        sdf = SDFITSLoad(str(fits_file))
        df = sdf.index()

    # Ensure ROW column exists (0-indexed row number)
    if "ROW" not in df.columns:
        df["ROW"] = range(len(df))

    # Parse UTC timestamps if present
    if "DATE-OBS" in df.columns:
        df["UTC"] = pd.to_datetime(df["DATE-OBS"], errors="coerce")
    elif "UTC" not in df.columns:
        # Generate fake timestamps if needed (1 second apart)
        df["UTC"] = pd.date_range(start="2025-01-01", periods=len(df), freq="1s")

    return df


def _check_disk_space(fits_files: list[Path], output_dir: Path) -> None:
    """Check if there's enough disk space for the simulation output.

    Raises IOError if insufficient space is available.
    """
    import shutil

    # Calculate total size of source files
    total_source_size = sum(f.stat().st_size for f in fits_files)

    # Add 10% safety margin for overhead
    required_space = int(total_source_size * 1.1)

    # Get available space on output filesystem
    disk_usage = shutil.disk_usage(output_dir)
    available_space = disk_usage.free

    if available_space < required_space:
        # Format sizes for human readability
        def format_size(size_bytes: int) -> str:
            for unit in ["B", "KB", "MB", "GB", "TB"]:
                if size_bytes < 1024:
                    return f"{size_bytes:.1f} {unit}"
                size_bytes /= 1024
            return f"{size_bytes:.1f} PB"

        raise IOError(
            f"Insufficient disk space in {output_dir}. "
            f"Need {format_size(required_space)}, but only {format_size(available_space)} available. "
            f"Free up space or use a different output directory."
        )

    logger.debug(
        f"Disk space check passed: need {total_source_size / 1e9:.1f} GB, "
        f"have {available_space / 1e9:.1f} GB available"
    )


def _initialize_index_file(index_file: Path, metadata: pd.DataFrame) -> None:
    """Create initial empty .index file with proper header structure."""
    # For now, don't create the index file until we have data
    # The append function will create it if needed
    pass


def _append_to_index_file(index_file: Path, scan_metadata: pd.DataFrame) -> None:
    """Append scan rows to .index file."""
    # TODO: Implement .index file writing
    # This would require formatting the DataFrame rows in the GBTIDL .index format
    # For MVP, we can skip this and rely on FITS-only simulation
    pass


def _initialize_status_file(status_file: Path) -> None:
    """Create initial sdfitsStatus.txt file with header."""
    with open(status_file, "w") as f:
        f.write("# GBT SDFITS Status File\n")
        f.write("# backend,project,scan,timestamp,datetime,file,index\n")
    logger.debug(f"Initialized status file: {status_file}")


def _update_status_file(
    status_file: Path,
    backend: str,
    project: str,
    scan: int,
    fits_file: Path,
    index_file: Path,
) -> None:
    """Update sdfitsStatus.txt with latest scan information."""
    timestamp = datetime.now()
    iso_timestamp = timestamp.strftime("%Y-%m-%dT%H:%M:%S")
    datetime_str = timestamp.strftime("%Y-%m-%d %H:%M:%S")

    # Format: backend,project,scan,timestamp,datetime,file,index
    # Use only the filename, not full path (like GBTIDL does)
    entry = (
        f"{backend},{project},{scan},{iso_timestamp},{datetime_str},"
        f"{fits_file.name},{index_file.name}\n"
    )

    # For now, overwrite with latest (could extend to track all files A,B,C,D)
    # In reality, each FITS file would have its own line, but for simplicity
    # we'll just update with the most recent scan
    with open(status_file, "a") as f:
        f.write(entry)
