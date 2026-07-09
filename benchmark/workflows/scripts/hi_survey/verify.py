"""
Usage: python verify.py dysh_stdout.txt gbtidl_stdout.txt

Parses RMS_BLUE / RMS_RED key=value lines from each script's captured stdout
and checks that values agree within rtol=0.02 (2%).
"""

import re
import sys

_PAT = re.compile(r"RMS_(BLUE|RED)[= ]+([0-9.eE+\-]+)")


def extract(path):
    vals = {}
    for m in _PAT.finditer(open(path).read()):
        vals[m.group(1)] = float(m.group(2))
    return vals


dysh = extract(sys.argv[1])
gbtidl = extract(sys.argv[2])

# The benchmark now uses explicit channel baseline windows, which removes the
# large mismatch from frequency-region handling. There is still a small
# residual dysh/GBTIDL semantics difference around the final RMS calculation,
# so keep the verify threshold at 2% until that underlying issue is fixed.
rtol = 0.02
ok = True
for key in ("BLUE", "RED"):
    d, g = dysh[key], gbtidl[key]
    rel = abs(d - g) / abs(g)
    status = "PASS" if rel < rtol else "FAIL"
    if rel >= rtol:
        ok = False
    print(f"{status} RMS_{key}: dysh={d:.4e}  gbtidl={g:.4e}  rel_diff={rel:.2%}")

sys.exit(0 if ok else 1)
