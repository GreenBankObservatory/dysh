import subprocess
import sys
import uuid
from io import StringIO

from dysh.shell import init_shell


def test_shell_cli():
    """Simply prove that we can launch $ dysh from CLI"""
    subprocess.check_call(["dysh"])


def test_shell_cli_with_args():
    """Simply prove that we can launch $ dysh from CLI"""
    subprocess.check_call(["dysh", "--colors", "Linux", "--no-banner", "--profile", "foo"])


def test_init_shell(monkeypatch):
    """Prove that we can open the shell, print something, and exit without error"""

    # Generate a unique string that we will look for in stdout later
    test_str = uuid.uuid4().hex.upper()
    mock_input = StringIO(f"print('{test_str}')\n")
    monkeypatch.setattr(sys, "stdin", mock_input)

    mock_stdout = StringIO()
    monkeypatch.setattr(sys, "stdout", mock_stdout)
    init_shell("--no-banner", colors="NoColor", profile="dysh")

    mock_stdout.seek(0)
    assert test_str in mock_stdout.read()
