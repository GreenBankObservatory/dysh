"""
Usage: python verify.py dysh_stdout.txt reference_stdout.txt

Parses TSYS_FDNUM_N=value lines from each script's captured stdout and checks
that the mean Tsys of each of the 16 ARGUS feeds agrees within an absolute
tolerance of 1 mK.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # benchmark/workflows/
import verify_common

NFEEDS = 16

verify_common.run(r"TSYS_FDNUM_\s*(\d+)\s*[= ]+\s*([0-9.eE+\-]+)", range(NFEEDS), "TSYS_FDNUM_", key=int)
