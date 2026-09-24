import os
import time

from dysh.fits import GBTFITSLoad
from dysh.log import logger

logger.info("START")
script_t0 = time.perf_counter()

filename = os.environ.get(
    "DYSH_BENCH_DATA_PATH",
    "/home/dysh/example_data/nod-KFPA/data/TGBT22A_503_02.raw.vegas",
)
sdf = GBTFITSLoad(filename)

print(f"DYSH_BENCH_SCRIPT_MS={(time.perf_counter() - script_t0) * 1000.0:.3f}")

logger.info("END")
