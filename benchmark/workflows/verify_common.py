"""
Shared logic for the workflow ``verify.py`` scripts.

Each verifier is invoked by ``run_bench.py`` as ``python verify.py dysh_stdout reference_stdout``,
where the reference is GBTIDL's stdout or a committed ``golden.txt``. A verifier extracts
``KEY=value`` markers from both files and calls `run`.
"""

import re
import sys

ATOL_K = 0.001  # absolute tolerance, 1 mK


def extract(path, pattern, key=str):
    """Read ``KEY=value`` markers from a captured stdout file.

    Parameters
    ----------
    path : str or `~pathlib.Path`
        File to read.
    pattern : str
        Regular expression with two groups: the key, then the numeric value.
    key : callable, optional
        Applied to group 1 to make the dictionary key, e.g. `int` for feed numbers.
        Default is `str`.

    Returns
    -------
    dict
        ``{key: float(value)}``. If a key appears more than once, the last value wins.

    Raises
    ------
    FileNotFoundError
        If `path` does not exist.
    """
    with open(path) as fh:
        return {key(m.group(1)): float(m.group(2)) for m in re.finditer(pattern, fh.read())}


def compare(dysh, ref, expected, label, fmt=".4f", atol=ATOL_K):
    """Compare dysh values against reference values and print a PASS/FAIL line for each key.

    A key missing from either side fails, so empty or garbled output cannot pass vacuously.
    Keys present in the data but not in `expected` are ignored.

    Parameters
    ----------
    dysh : dict
        Values from dysh's output, as returned by `extract`.
    ref : dict
        Values from the reference (GBTIDL or golden) output.
    expected : iterable
        Keys that must be present on both sides and agree.
    label : str
        Prefix for each printed name, e.g. ``"RMS_"``.
    fmt : str, optional
        Format spec for the printed values. Default is ``".4f"``.
    atol : float, optional
        Absolute tolerance, in the same units as the values (K). A difference equal to `atol`
        fails. Default is `ATOL_K`.

    Returns
    -------
    bool
        `True` if every expected key is present on both sides and within tolerance.
    """
    ok = True
    for k in expected:
        name = f"{label}{k}"
        if k not in dysh or k not in ref:
            side = "dysh" if k not in dysh else "reference"
            print(f"FAIL {name}: missing from {side} output")
            ok = False
            continue
        d, r = dysh[k], ref[k]
        diff = abs(d - r)
        passed = diff < atol
        ok &= passed
        print(f"{'PASS' if passed else 'FAIL'} {name}: dysh={d:{fmt}}  ref={r:{fmt}}  abs_diff={diff * 1000:.3f} mK")
    return ok


def run(pattern, expected, label, key=str, fmt=".4f"):
    """Run a verifier: compare two captured stdout files and exit with the result.

    Reads the dysh and reference paths from ``sys.argv[1]`` and ``sys.argv[2]``, as passed by
    ``run_bench.py``.

    Parameters
    ----------
    pattern : str
        Regular expression with two groups: the key, then the numeric value.
    expected : iterable
        Keys that must be present on both sides and agree.
    label : str
        Prefix for each printed name.
    key : callable, optional
        Applied to group 1 of `pattern` to make the key. Default is `str`.
    fmt : str, optional
        Format spec for the printed values. Default is ``".4f"``.

    Raises
    ------
    SystemExit
        Always: exit code 0 if all checks pass, 1 otherwise.
    """
    dysh = extract(sys.argv[1], pattern, key)
    ref = extract(sys.argv[2], pattern, key)
    sys.exit(0 if compare(dysh, ref, expected, label, fmt=fmt) else 1)
