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
    """Return ``{key(group 1): float(group 2)}`` for every match of ``pattern`` in the file."""
    with open(path) as fh:
        return {key(m.group(1)): float(m.group(2)) for m in re.finditer(pattern, fh.read())}


def compare(dysh, ref, expected, label, fmt=".4f", atol=ATOL_K):
    """Print one PASS/FAIL line per expected key; return True if all pass.

    A key missing from either side fails, so an empty or garbled output cannot pass vacuously.
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
    """Entry point for a verifier: read the two paths from ``sys.argv`` and exit 0 on success."""
    dysh = extract(sys.argv[1], pattern, key)
    ref = extract(sys.argv[2], pattern, key)
    sys.exit(0 if compare(dysh, ref, expected, label, fmt=fmt) else 1)
