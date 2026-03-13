import os
import time

import numpy as np
from astropy import units as u
from astropy.time import Time

from dysh.fits import GBTFITSLoad
from dysh.log import logger
from dysh.util.weatherforecast import GBTWeatherForecast

logger.info("START")
script_t0 = time.perf_counter()

path = os.environ.get(
    "DYSH_BENCH_DATA_PATH",
    "/home/scratch/ajschmie/training/dysh/datasets/argus/TGBT22A_603_05_vanecal.raw.vegas",
)
scan1 = 10
ifnum = 0
plnum = 0
center_fdnum = 9

stage_t0 = time.perf_counter()
sdfits = GBTFITSLoad(path)
print(f"DYSH_BENCH_STAGE_MS[GBTFITSLoad]={(time.perf_counter() - stage_t0) * 1000.0:.3f}")

sky_scan = scan1 + 1
stage_t0 = time.perf_counter()
if hasattr(sdfits, "_prefetch_vanecal_metadata"):
    sdfits._prefetch_vanecal_metadata(vane_scan=scan1, sky_scan=sky_scan, ifnum=ifnum, plnum=plnum)

center = None
if hasattr(sdfits, "_get_vanecal_scan_cache"):
    scan_cache = sdfits._get_vanecal_scan_cache(
        vane_scan=scan1,
        sky_scan=sky_scan,
        ifnum=ifnum,
        plnum=plnum,
        apply_flags=True,
    )
    center = scan_cache[scan1][center_fdnum]
else:
    center = {
        "meta": sdfits.gettp(
            scan=scan1,
            ifnum=ifnum,
            plnum=plnum,
            fdnum=center_fdnum,
            calibrate=True,
            cal=False,
            apply_flags=True,
        ).timeaverage(use_wcs=False).meta
    }

center_meta = center["meta"]
twarm = center_meta["TWARM"] + 273.15
if twarm >= 370.0:
    twarm = center_meta["TAMBIENT"] + 1.5

specval = float(center_meta["OBSFREQ"]) * u.Hz
mjd = Time(center_meta["DATE-OBS"]).mjd
try:
    gbwf = GBTWeatherForecast()
    tau = gbwf.fetch(vartype="Opacity", specval=specval, mjd=mjd)[:, -1]
    tatm = gbwf.fetch(vartype="Tatm", specval=specval, mjd=mjd)[:, -1]
    airmass = 1.0 / np.sin(np.deg2rad(float(center_meta["ELEVATIO"])))
    tbg = 2.725
    tcal = (tatm - tbg) + (twarm - tatm) * np.exp(tau * airmass)
except ValueError:
    tcal = center_meta["TAMBIENT"]
print(f"DYSH_BENCH_STAGE_MS[setup_tcal]={(time.perf_counter() - stage_t0) * 1000.0:.3f}")

stage_t0 = time.perf_counter()
if hasattr(sdfits, "_get_vanecal_scan_cache"):
    sdfits._get_vanecal_scan_cache(vane_scan=scan1, sky_scan=sky_scan, ifnum=ifnum, plnum=plnum, apply_flags=True)
print(f"DYSH_BENCH_STAGE_MS[prefetch_vanecal_cache]={(time.perf_counter() - stage_t0) * 1000.0:.3f}")

feed_loop_t0 = time.perf_counter()
for i in range(16):
    stage_t0 = time.perf_counter()
    print("FDnum = ", i)
    tsys = sdfits.vanecal(scan1, fdnum=i, ifnum=ifnum, plnum=plnum, tcal=tcal)
    print(f"TSYS_FDNUM_{i}={np.mean(tsys):.4f}")
    print("   tsys = ", tsys)
    print(f"DYSH_BENCH_STAGE_MS[fdnum_{i}]={(time.perf_counter() - stage_t0) * 1000.0:.3f}")
print(f"DYSH_BENCH_STAGE_MS[feed_loop_total]={(time.perf_counter() - feed_loop_t0) * 1000.0:.3f}")

print(f"DYSH_BENCH_SCRIPT_MS={(time.perf_counter() - script_t0) * 1000.0:.3f}")

logger.info("END")
