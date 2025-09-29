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


def test_shell_cli_with_script(tmp_path):
    """Test that we can run a script from the CLI"""

    script = tmp_path / "test_cli_script.py"
    code = """fnm = dysh_data(test='getps')
sdf = GBTFITSLoad(fnm)
sb = sdf.getps(ifnum=0, fdnum=0, plnum=0)
assert len(sb) == 4"""

    with open(script, "w") as o:
        o.write(code)
    subprocess.check_call(["dysh", script])
