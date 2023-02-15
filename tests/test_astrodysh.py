#!/usr/bin/env python

"""Tests for `astrodysh` package."""

import subprocess

from astrodysh.astrodysh import astrodysh


def test_cli():
    """Ensure that the CLI has been installed"""
    subprocess.check_call(["astrodysh", "foo", "bar"])


def test_hello_world():
    assert astrodysh("foo", "bar") is True
