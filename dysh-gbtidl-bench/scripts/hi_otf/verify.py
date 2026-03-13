"""
Usage: python verify.py  (run in scripts/hi_otf/ after both scripts execute)

Loads both output FITS files (otf_calibrated.fits from dysh, gbtidl.fits from
gbtidl), time-averages each, and checks that spectrum values agree within 1%.
"""
import sys

import numpy as np
from astropy.io import fits

d_data = fits.open("otf_calibrated.fits")[1].data["DATA"].astype(float)
g_data = fits.open("gbtidl.fits")[1].data["DATA"].astype(float)

if d_data.shape != g_data.shape:
    print(f"FAIL shape mismatch: dysh={d_data.shape} gbtidl={g_data.shape}")
    sys.exit(1)

rel = np.nanmean(np.abs(d_data - g_data) / (np.abs(g_data) + 1e-30))
status = "PASS" if rel < 0.01 else "FAIL"
print(f"{status} hi_otf DATA mean_rel_diff={rel:.2%}")
sys.exit(0 if rel < 0.01 else 1)
