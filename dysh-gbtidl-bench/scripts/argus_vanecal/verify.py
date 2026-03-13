"""
Usage: python verify.py dysh_stdout.txt gbtidl_stdout.txt

Parses TSYS_FDNUM_N=value lines from dysh's captured stdout and
TSYS_FDNUM_N=value lines from gbtidl's captured stdout, then checks
that mean Tsys per feed agrees within rtol=0.02 (2%).
"""

import re
import sys

_PAT = re.compile(r"TSYS_FDNUM_\s*(\d+)\s*[= ]+\s*([0-9.eE+\-]+)")


def extract(path):
    vals = {}
    for m in _PAT.finditer(open(path).read()):
        vals[int(m.group(1))] = float(m.group(2))
    return vals


dysh = extract(sys.argv[1])
gbtidl = extract(sys.argv[2])

rtol = 0.02
ok = True
feeds = sorted(set(dysh) | set(gbtidl))
for fdnum in feeds:
    if fdnum not in dysh:
        print(f"FAIL FDNUM_{fdnum}: missing from dysh output")
        ok = False
        continue
    if fdnum not in gbtidl:
        print(f"FAIL FDNUM_{fdnum}: missing from gbtidl output")
        ok = False
        continue
    d, g = dysh[fdnum], gbtidl[fdnum]
    rel = abs(d - g) / abs(g)
    status = "PASS" if rel < rtol else "FAIL"
    if rel >= rtol:
        ok = False
    print(f"{status} TSYS_FDNUM_{fdnum}: dysh={d:.4f}  gbtidl={g:.4f}  rel_diff={rel:.2%}")

sys.exit(0 if ok else 1)
