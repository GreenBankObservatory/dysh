import os
import time

from dysh.fits import GBTFITSLoad
from dysh.log import logger

logger.info("START")
script_t0 = time.perf_counter()

filename = os.environ.get("DYSH_BENCH_DATA_PATH", "/home/scratch/dfrayer/DATAdemo/TGBT17A_506_11.raw.vegas")
ifnum = 0  # The 21 cm line is in the spectral window labeled 0.
fdnum = 0  # Only one feed in this data set
ref = 27  # The reference ("OFF") scan
scans = list(range(14, 27))  # The signal ("ON") scans

sdfits = GBTFITSLoad(filename)

sb0 = sdfits.getsigref(scan=scans, ref=ref, fdnum=fdnum, ifnum=ifnum, plnum=0)

out_path = os.environ.get("DYSH_BENCH_OUT_PATH", "otf_calibrated.fits")
sb0.write(out_path, overwrite=True)

print(f"DYSH_BENCH_SCRIPT_MS={(time.perf_counter() - script_t0) * 1000.0:.3f}")

logger.info("END")
