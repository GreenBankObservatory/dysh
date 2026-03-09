#!/usr/bin/env python
"""
dysh-cli: Command-line interface for dysh utilities

Provides various command-line tools for working with SDFITS data.
"""

import argparse
import logging
import os
import sys
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

from dysh.fits import GBTOnline, gbtfitsload, index_file, sdfitsload
from dysh.fits.gbtfitsload import _parse_sdfits_status_file
from dysh.log import logger
from dysh.shell import main as shell_main
from dysh.util.online_simulator import simulate_session


def format_size(size_bytes):
    """Format size in bytes to human-readable format."""
    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if size_bytes < 1024.0:
            return f"{size_bytes:.1f} {unit}"
        size_bytes /= 1024.0
    return f"{size_bytes:.1f} PB"


def get_total_size(filepath):
    """Get total size of FITS file(s) in bytes."""
    path = Path(filepath)
    total_size = 0

    if path.is_file():
        total_size = path.stat().st_size
    elif path.is_dir():
        # Sum all .fits files in directory
        for fits_file in path.glob("*.fits"):
            total_size += fits_file.stat().st_size

    return total_size


def cmd_info(args):
    """Display information about a FITS file."""
    filepath = Path(args.file)
    if not filepath.exists():
        logger.error(f"File not found: {filepath}")
        return 1

    # Try to load the file
    try:
        start_time = time.time()

        # Check for index file
        index_path = index_file.get_index_path(filepath)
        has_index = index_path.exists()

        if args.verbosity >= 1:
            if has_index:
                print(f"📋 .index file: {index_path}")
            else:
                print("📋 No .index file (will read from FITS)")
            print()

        if filepath.is_dir() or "vegas" in str(filepath).lower():
            # Multi-file dataset
            loader = gbtfitsload.GBTFITSLoad(str(filepath), index_only=True)
        else:
            # Single FITS file
            loader = sdfitsload.SDFITSLoad(str(filepath), index_only=True)

        load_time = time.time() - start_time

        print(f"File: {filepath}")
        print(f"Total rows: {loader.total_rows}")
        # Handle both SDFITSLoad (_bintable) and GBTFITSLoad (_sdf list)
        if hasattr(loader, "_bintable"):
            num_bintables = len(loader._bintable)
        elif hasattr(loader, "_sdf"):
            # GBTFITSLoad: count total bintables across all SDFITSLoad objects
            num_bintables = sum(len(sdf._bintable) for sdf in loader._sdf)
        else:
            num_bintables = 0
        print(f"Number of bintables: {num_bintables}")
        if has_index:
            print(f"Index: ⚡ Loaded from .index file ({load_time:.2f}s)")
        else:
            print(f"Index: 📖 Read from FITS ({load_time:.2f}s)")
        print()

        if args.verbosity >= 2:
            loader.summary()
        else:
            print("Columns:", ", ".join(loader.columns[:20]))
            if len(loader.columns) > 20:
                print(f"  ... and {len(loader.columns) - 20} more")

        return 0
    except Exception as e:
        logger.error(f"Error loading file: {e}")
        traceback.print_exc()
        return 1


def cmd_index(args):
    """Create or display index file for a FITS file."""
    filepath = Path(args.file)
    if not filepath.exists():
        logger.error(f"File not found: {filepath}")
        return 1

    index_path = index_file.get_index_path(filepath)

    if args.create:
        # Create new index file
        try:
            loader = sdfitsload.SDFITSLoad(str(filepath))
            loader.create_index()

            if loader._index_metadata is None:
                loader._index_metadata = index_file.create_index_metadata(
                    observer=loader._header.get("OBSERVER", "Unknown"),
                    backend=loader._header.get("BACKEND", "Unknown"),
                )

            index_file.write_index(index_path, loader._index_metadata, loader._index)
            print(f"Created index file: {index_path}")
            print(f"  {len(loader._index)} rows, {len(loader._index.columns)} columns")
            return 0
        except Exception as e:
            logger.error(f"Error creating index: {e}")
            traceback.print_exc()
            return 1

    elif args.validate:
        # Validate existing index
        if index_file.validate_index(filepath, index_path):
            print(f"Index file is valid: {index_path}")
            return 0
        else:
            print(f"Index file is invalid or out of date: {index_path}")
            return 1

    else:
        # Display index info
        if not index_path.exists():
            print(f"No index file found at: {index_path}")
            print("Use --create to create one")
            return 1

        try:
            metadata, df = index_file.read_index(index_path)
            print(f"Index file: {index_path}")
            print(f"Created: {metadata.created}")
            print(f"Observer: {metadata.observer}")
            print(f"Backend: {metadata.backend}")
            print(f"Rows: {len(df)}")
            print(f"Columns: {len(df.columns)}")

            if args.verbosity >= 1:
                print()
                print("First 10 rows:")
                print(df.head(10))

            return 0
        except Exception as e:
            logger.error(f"Error reading index: {e}")
            traceback.print_exc()
            return 1


def cmd_summary(args):
    """Display detailed summary of a FITS file."""
    filepath = Path(args.file)
    if not filepath.exists():
        logger.error(f"File not found: {filepath}")
        return 1

    try:
        start_time = time.time()

        # Check for index file
        index_path = index_file.get_index_path(filepath)
        has_index = index_path.exists()

        # Load the file
        if has_index:
            if args.verbosity >= 1:
                print(f"Found index file; loading {index_path}...")
            size_to_load = get_total_size(index_path)
        else:
            if args.verbosity >= 1:
                print(f"No index file found; loading {filepath}...")
            size_to_load = get_total_size(filepath)

        load_start = time.time()
        loader = gbtfitsload.GBTFITSLoad(str(filepath), index_only=True)
        load_time = time.time() - load_start

        if args.verbosity >= 1:
            print(f"Loaded {format_size(size_to_load)} in {load_time:.2f}s")
            print(f"  Total rows: {loader.total_rows}")
            print(f"  FITS file size: {format_size(get_total_size(filepath))}")
            # Handle both SDFITSLoad (_bintable) and GBTFITSLoad (_sdf list)
            if hasattr(loader, "_bintable"):
                num_bintables = len(loader._bintable)
            elif hasattr(loader, "_sdf"):
                # GBTFITSLoad: count total bintables across all SDFITSLoad objects
                num_bintables = sum(len(sdf._bintable) for sdf in loader._sdf)
            else:
                num_bintables = 0
            print(f"  Bintables: {num_bintables}")
            print(f"  Index columns: {len(loader._index.columns)}")
            print()

        # Generate summary
        if args.verbosity >= 1:
            print("Generating summary...")

        # Just print directly - users can pipe to less if they want paging
        if hasattr(loader, "get_summary"):
            # GBTFITSLoad has rich summary options
            if args.verbosity >= 2:
                cols_requested = args.columns.split(",") if args.columns else "default"
                print(f"   Columns: {cols_requested}")
                print("   Accessing metadata from index (no DATA column read)")
                print()

            loader.summary(
                scan=args.scan,
                verbose=args.verbosity >= 2,
                max_rows=args.max_rows,
                show_index=args.show_index,
                columns=args.columns.split(",") if args.columns else None,
                add_columns=args.add_columns.split(",") if args.add_columns else None,
            )
        else:
            # SDFITSLoad has simple summary
            if args.verbosity >= 2:
                print("   Using simple summary format")
                print()
            loader.summary()

        if args.verbosity >= 1:
            print()
            total_time = time.time() - start_time
            print(f"Summary generated in {total_time:.2f}s")

        return 0
    except Exception as e:
        logger.error(f"Error loading file: {e}")
        traceback.print_exc()
        return 1


def cmd_shell(args):
    """Launch dysh IPython shell."""
    # Reconstruct sys.argv for the shell command
    # Remove 'shell' subcommand but keep all other arguments
    original_argv = sys.argv.copy()
    sys.argv = [sys.argv[0], *args.shell_args]

    try:
        shell_main()
        return 0
    finally:
        sys.argv = original_argv


def cmd_online(args):
    """Monitor live GBT observations."""
    # Determine sdfits_root for session checking
    if "SDFITS_DATA" in os.environ:
        sdfits_root = os.environ["SDFITS_DATA"]
    elif "DYSH_DATA" in os.environ:
        sdfits_root = os.environ["DYSH_DATA"] + "/sdfits"
    else:
        sdfits_root = "/home/sdfits"
    status_file_path = os.path.join(sdfits_root, "sdfitsStatus.txt")

    # Connect to observation
    print("Connecting to GBT online data...")
    try:
        if args.project:
            sdf = GBTOnline(args.project, backend=args.backend)
            # If user specified a project explicitly, don't check for session changes during monitoring
            check_session_changes = False
            monitored_path = None

            # But DO check if they connected to a non-active session at startup
            if os.path.exists(status_file_path):
                try:
                    status_info = _parse_sdfits_status_file(status_file_path, backend=args.backend)
                    if status_info:
                        active_session = os.path.join(sdfits_root, status_info["project"], status_info["file"])
                        if sdf._online != active_session:
                            print("\n" + "=" * 80)
                            print("⚠️  WARNING: NOT MONITORING ACTIVE SESSION ⚠️")
                            print("=" * 80)
                            print(f"You are monitoring: {sdf._online}")
                            print(f"Active session is:  {active_session}")
                            print(f"Active backend: {status_info['backend'].upper()}")
                            print("\nTo monitor the active session instead, use:")
                            print("  dysh online")
                            print("=" * 80 + "\n")
                except Exception as e:
                    # Don't fail if we can't check - just skip the warning
                    logger.debug(f"Could not check active session: {e}")
        else:
            sdf = GBTOnline(backend=args.backend)
            # Auto-discovered mode: check for session changes during monitoring
            check_session_changes = True
            monitored_path = sdf._online
    except FileNotFoundError as e:
        print(f"\nError: {e}", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"\nError: Could not connect: {e}", file=sys.stderr)
        if args.verbosity >= 1:
            traceback.print_exc()
        return 1

    # Parse columns if provided
    # For online mode, we want to show DATE-OBS (scan start time) by default
    if args.columns:
        columns = [c.strip() for c in args.columns.split(",")]
        add_columns = None
    else:
        columns = None
        add_columns = ["DATE-OBS"]

    # Get initial summary DataFrame
    try:
        df = sdf.get_summary(columns=columns, add_columns=add_columns)

        # Add timestamp column showing when we detected this data
        df["RENDERED"] = datetime.now(timezone.utc).strftime("%H:%M:%S")

        last_summary_count = len(df)
        last_raw_count = len(sdf._index) if sdf._index is not None else 0

        # Track seen scan numbers to detect session restarts
        seen_scans = set(df["SCAN"].unique()) if "SCAN" in df.columns else set()

        # Track last scan number to detect scan regression
        last_scan_num = int(df["SCAN"].iloc[-1]) if "SCAN" in df.columns and len(df) > 0 else None

        # Print initial scans
        print(f"\nMonitoring: {sdf._online}")
        print(f"Polling every {args.interval} seconds. Press Ctrl+C to stop.\n")
        print(df.to_string(index=False))
        print(f"\n[{last_raw_count} raw rows, {last_summary_count} scans]")
        print()

    except Exception as e:
        print(f"\nError getting initial summary: {e}", file=sys.stderr)
        if args.verbosity >= 1:
            traceback.print_exc()
        return 1

    # Monitor loop
    try:
        while True:
            time.sleep(args.interval)

            # Check if session has changed (if in auto-discover mode)
            if check_session_changes and os.path.exists(status_file_path):
                try:
                    status_info = _parse_sdfits_status_file(status_file_path, backend=args.backend)
                    if status_info:
                        current_session = os.path.join(sdfits_root, status_info["project"], status_info["file"])
                        if current_session != monitored_path:
                            print("\n" + "=" * 80)
                            print("⚠️  SESSION CHANGED! ⚠️")
                            print("=" * 80)
                            print(f"Currently monitoring: {monitored_path}")
                            print(f"Active session is now: {current_session}")
                            print(f"Backend: {status_info['backend'].upper()}")
                            print("\nTo switch to the new session, restart with:")
                            print("  dysh online")
                            if args.backend:
                                print(f"  or: dysh online --backend {args.backend}")
                            print("=" * 80 + "\n")
                except Exception as e:
                    # Don't let session check errors break monitoring
                    logger.debug(f"Error checking for session changes: {e}")

            # Get updated summary
            try:
                # First, reload new data from disk
                if args.verbosity >= 3:
                    print(
                        f"[DEBUG] Checking for new data at {datetime.now(timezone.utc).strftime('%H:%M:%S')}",
                        file=sys.stderr,
                    )

                has_new_data = sdf._reload()

                if args.verbosity >= 3:
                    print(f"[DEBUG] _reload() returned: {has_new_data}", file=sys.stderr)

                # Check raw row count to detect new integrations within same scan
                current_raw_count = len(sdf._index) if sdf._index is not None else 0

                # Detect if file was reset (fewer rows than before)
                if current_raw_count < last_raw_count:
                    timestamp = datetime.now(timezone.utc).strftime("%H:%M:%S")
                    logger.warning(
                        f"Data reset detected: {last_raw_count} → {current_raw_count} rows. "
                        "File may have been recreated."
                    )
                    print(
                        f"[{timestamp}] Data reset: {last_raw_count} → {current_raw_count} rows (file recreated)",
                        flush=True,
                    )
                    # Reset counters to track from the new file
                    last_raw_count = current_raw_count
                    last_summary_count = len(sdf.get_summary(columns=columns, add_columns=add_columns))
                    seen_scans = set()  # Reset seen scans on data reset
                    last_scan_num = None  # Reset last scan on data reset

                if has_new_data or current_raw_count > last_raw_count:
                    df = sdf.get_summary(columns=columns, add_columns=add_columns)

                    # Add timestamp column showing when we detected this data
                    df["RENDERED"] = datetime.now(timezone.utc).strftime("%H:%M:%S")

                    current_summary_count = len(df)

                    # Detect duplicate scan numbers (session restart indicator)
                    if "SCAN" in df.columns:
                        current_scans = set(df["SCAN"].unique())
                        new_scan_nums = current_scans - seen_scans
                        if not new_scan_nums and current_summary_count > last_summary_count:
                            # We have new rows but no new scan numbers - must be duplicates
                            timestamp = datetime.now(timezone.utc).strftime("%H:%M:%S")
                            logger.warning("Duplicate scan numbers detected. Session may have restarted.")
                            print(
                                f"[{timestamp}] WARNING: Duplicate scan numbers detected (session may have restarted)",
                                flush=True,
                            )
                        seen_scans.update(current_scans)

                        # Detect scan number regression (current scan lower than previous)
                        current_scan_num = int(df["SCAN"].iloc[-1])
                        if last_scan_num is not None and current_scan_num < last_scan_num:
                            timestamp = datetime.now(timezone.utc).strftime("%H:%M:%S")
                            logger.warning(
                                f"Scan number regression: {last_scan_num} → {current_scan_num}. "
                                "Unexpected scan sequence."
                            )
                            print(
                                f"[{timestamp}] WARNING: Scan number regression: "
                                f"{last_scan_num} → {current_scan_num} (unexpected sequence)",
                                flush=True,
                            )
                        last_scan_num = current_scan_num

                    # Print new scan rows if any
                    if current_summary_count > last_summary_count:
                        new_rows = df.iloc[last_summary_count:]
                        print(new_rows.to_string(index=False, header=False))
                        last_summary_count = current_summary_count

                    # Always show status update when raw data changes
                    if current_raw_count > last_raw_count:
                        timestamp = datetime.now(timezone.utc).strftime("%H:%M:%S")
                        latest_scan = df["SCAN"].iloc[-1] if "SCAN" in df.columns else "?"
                        print(
                            f"[{timestamp}] Scan {latest_scan}: "
                            f"{last_raw_count} → {current_raw_count} rows (+{current_raw_count - last_raw_count})",
                            flush=True,
                        )
                        last_raw_count = current_raw_count
                    elif has_new_data:
                        # File changed but no new rows - still notify
                        timestamp = datetime.now(timezone.utc).strftime("%H:%M:%S")
                        print(f"[{timestamp}] Data refreshed (no new rows)", flush=True)

            except Exception as e:
                logger.warning(f"Error during update: {e}")
                if args.verbosity >= 2:
                    traceback.print_exc()
                continue

    except KeyboardInterrupt:
        print("\n\nMonitoring stopped.")
        return 0

    return 0


def cmd_simulate_online(args):
    """Simulate SDFITS session creation for testing online mode."""
    try:
        print(f"Simulating session: {args.source_dir} → {args.output_dir}")
        print(f"Mode: {args.mode}, interval: {args.interval}s, speedup: {args.speedup}x")
        print(f"Parallel: {not args.sequential}, Create index: {args.create_index}")
        print()

        simulate_session(
            source_dir=args.source_dir,
            output_dir=args.output_dir,
            mode=args.mode,
            interval=args.interval,
            speedup=args.speedup,
            create_index=args.create_index,
            parallel=not args.sequential,
        )

        print()
        print(f"✓ Simulation complete! Output: {args.output_dir}")
        return 0

    except Exception as e:
        logger.error(f"Simulation failed: {e}")
        if args.verbosity >= 1:
            raise
        return 1


def main():
    """Main CLI entry point."""
    # Special handling for 'shell' command: pass all arguments directly to shell
    # This is needed because argparse.REMAINDER doesn't work well with -- prefixed args
    # that look like options (e.g., --colors, --profile) - argparse tries to parse them
    # at the parent level before REMAINDER gets a chance to capture them.
    if len(sys.argv) > 1 and sys.argv[1] == "shell":
        # Pass everything after 'shell' directly to the shell's parser
        original_argv = sys.argv.copy()
        sys.argv = [sys.argv[0], *sys.argv[2:]]  # Remove 'shell' from argv
        try:
            shell_main()
            return 0
        finally:
            sys.argv = original_argv

    parser = argparse.ArgumentParser(
        prog="dysh",
        description="Command-line utilities for dysh SDFITS data processing",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "-v",
        action="count",
        default=0,
        dest="verbosity",
        help="Increase verbosity (-v, -vv, -vvv for levels 1, 2, 3)",
    )

    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # Shell command - Interactive IPython shell
    # Note: actual shell handling is done above, before argparse, to properly pass through args
    shell_parser = subparsers.add_parser("shell", help="Launch dysh interactive IPython shell")
    shell_parser.add_argument("shell_args", nargs=argparse.REMAINDER, help="Arguments to pass to dysh shell")
    shell_parser.set_defaults(func=cmd_shell)

    # Info command
    info_parser = subparsers.add_parser("info", help="Display information about a FITS file")
    info_parser.add_argument("file", help="FITS file or directory")
    info_parser.set_defaults(func=cmd_info)

    # Summary command
    summary_parser = subparsers.add_parser("summary", help="Display detailed summary of a FITS file")
    summary_parser.add_argument("file", help="FITS file or directory")
    summary_parser.add_argument("--scan", type=int, help="Show specific scan number")
    summary_parser.add_argument(
        "--max-rows", type=int, default=-1, help="Maximum rows to display (-1=config default, None=unlimited)"
    )
    summary_parser.add_argument("--show-index", action="store_true", help="Show DataFrame index")
    summary_parser.add_argument("--columns", type=str, help="Comma-separated list of columns to display")
    summary_parser.add_argument("--add-columns", type=str, help="Comma-separated list of columns to add to defaults")
    summary_parser.set_defaults(func=cmd_summary)

    # Index command
    index_parser = subparsers.add_parser("index", help="Create or display index file")
    index_parser.add_argument("file", help="FITS file")
    index_parser.add_argument("-c", "--create", action="store_true", help="Create new index file")
    index_parser.add_argument("--validate", action="store_true", help="Validate existing index file")
    index_parser.set_defaults(func=cmd_index)

    # Online command - Monitor live observations
    online_parser = subparsers.add_parser("online", help="Monitor live GBT observations")
    online_parser.add_argument(
        "project", nargs="?", help="Project name or directory to monitor (default: auto-discover)"
    )
    online_parser.add_argument("--backend", choices=["vegas", "acs", "sp"], help="Specific backend to monitor")
    online_parser.add_argument("--interval", type=int, default=5, help="Polling interval in seconds (default: 5)")
    online_parser.add_argument("--columns", type=str, help="Comma-separated list of columns to display")
    online_parser.set_defaults(func=cmd_online)

    # Simulate-online command - Simulate SDFITS session creation
    simulate_parser = subparsers.add_parser(
        "simulate-online", help="Simulate SDFITS session creation scan-by-scan for testing"
    )
    simulate_parser.add_argument("source_dir", help="Source SDFITS session directory (e.g., PROJECT.raw.vegas/)")
    simulate_parser.add_argument("output_dir", help="Output directory for simulated session")
    simulate_parser.add_argument(
        "--mode",
        choices=["fixed", "realtime"],
        default="fixed",
        help="Timing mode: 'fixed' interval or 'realtime' with original timing (default: fixed)",
    )
    simulate_parser.add_argument(
        "--interval", type=float, default=1.0, help="For fixed mode: seconds between scans (default: 1.0)"
    )
    simulate_parser.add_argument(
        "--speedup",
        type=float,
        default=1.0,
        help="For realtime mode: speed multiplier (e.g., 10.0 = 10x faster) (default: 1.0)",
    )
    simulate_parser.add_argument(
        "--no-index",
        action="store_false",
        dest="create_index",
        help="Don't create .index files (FITS only)",
    )
    simulate_parser.add_argument(
        "--sequential",
        action="store_true",
        help="Write FITS files sequentially instead of in parallel",
    )
    simulate_parser.set_defaults(func=cmd_simulate_online)

    args = parser.parse_args()

    # Configure logging based on verbosity level
    # Default: INFO level with timestamps
    # -v: DEBUG level
    # -vv: DEBUG level with more detail
    if args.verbosity >= 2:
        log_level = logging.DEBUG
        log_format = "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
    elif args.verbosity >= 1:
        log_level = logging.DEBUG
        log_format = "%(asctime)s [%(levelname)s] %(message)s"
    else:
        log_level = logging.INFO
        log_format = "%(asctime)s %(message)s"

    # Configure root logger for CLI
    logging.basicConfig(
        level=log_level,
        format=log_format,
        datefmt="%H:%M:%S",
        stream=sys.stderr,
    )

    if not hasattr(args, "func"):
        parser.print_help()
        return 1

    try:
        return args.func(args)
    except KeyboardInterrupt:
        print("\nInterrupted by user")
        return 130
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        if args.verbosity >= 1:
            raise
        return 1


if __name__ == "__main__":
    sys.exit(main())
