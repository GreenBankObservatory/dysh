"""
Usage: python verify.py dysh_stdout.txt reference_stdout.txt

Parses RMS_BLUE / RMS_RED key=value lines from each script's captured stdout
and checks that values agree within an absolute tolerance of 1 mK.

1 mK is used rather than a relative tolerance because a small residual
dysh/GBTIDL difference remains in the final RMS calculation (the explicit
channel baseline windows in the benchmark already removed the larger mismatch
from frequency-region handling).
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # benchmark/workflows/
import verify_common

verify_common.run(r"RMS_(BLUE|RED)[= ]+([0-9.eE+\-]+)", ("BLUE", "RED"), "RMS_", fmt=".4e")
