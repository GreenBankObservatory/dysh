"""Tests for dysh.log module."""

import logging
from pathlib import Path

from dysh.log import init_logging


def test_init_logging_with_unwritable_path(tmp_path):
    """Test that init_logging doesn't crash with an unwritable path."""
    # Create a read-only directory
    readonly_dir = tmp_path / "readonly"
    readonly_dir.mkdir()
    readonly_dir.chmod(0o444)

    log_file = readonly_dir / "test.log"

    # This should not raise an exception
    init_logging(verbosity=2, path=log_file)

    # Verify logging still works (to stderr)
    logger = logging.getLogger("dysh")
    logger.info("Test message")

    # Clean up
    readonly_dir.chmod(0o755)


def test_init_logging_with_nonexistent_directory():
    """Test that init_logging doesn't crash with a non-existent directory."""
    # Use a path that definitely doesn't exist
    log_file = Path("/nonexistent/directory/test.log")

    # This should not raise an exception
    init_logging(verbosity=2, path=log_file)

    # Verify logging still works (to stderr)
    logger = logging.getLogger("dysh")
    logger.info("Test message")


def test_init_logging_with_valid_path(tmp_path):
    """Test that init_logging works correctly with a valid path."""
    log_file = tmp_path / "test.log"

    # This should work without issues
    init_logging(verbosity=2, path=log_file)

    # Verify logging works
    logger = logging.getLogger("dysh")
    logger.info("Test message")

    # The log file should exist (after delay) when we actually log something
    # Note: with delay=True, the file won't exist until first write


def test_init_logging_without_path():
    """Test that init_logging works without a path."""
    # This should work without issues
    init_logging(verbosity=2)

    # Verify logging works
    logger = logging.getLogger("dysh")
    logger.info("Test message")
