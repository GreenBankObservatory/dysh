import sys
from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.time import Time

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))  # benchmark/
from bench_common import MarkerDTime, resolve_data
from dysh.fits import GBTFITSLoad
from dysh.log import logger
from dysh.util.weatherforecast import GBTWeatherForecast

logger.info("START")

dt = MarkerDTime(benchname="argus_vanecal")

path = resolve_data(
    "/home/scratch/ajschmie/training/dysh/datasets/argus/TGBT22A_603_05_vanecal.raw.vegas", example="otf4"
)
scan1 = 10
ifnum = 0
plnum = 0
center_fdnum = 9

sdfits = GBTFITSLoad(path)
dt.tag("GBTFITSLoad")

center_meta = (
    sdfits.gettp(
        scan=scan1,
        ifnum=ifnum,
        plnum=plnum,
        fdnum=center_fdnum,
        calibrate=True,
        cal=False,
        apply_flags=True,
    )
    .timeaverage(use_wcs=False)
    .meta
)

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
    # Weather forecasts are only reachable at GBO; fall back to ambient temperature.
    tcal = center_meta["TAMBIENT"]
dt.tag("setup_tcal")

for i in range(16):
    print("FDnum = ", i)
    tsys = sdfits.vanecal(scan1, fdnum=i, ifnum=ifnum, plnum=plnum, tcal=tcal)
    print(f"TSYS_FDNUM_{i}={np.mean(tsys):.4f}")
    print("   tsys = ", tsys)
    dt.tag(f"fdnum_{i}")
dt.span("feed_loop_total", since="setup_tcal")

dt.report()

logger.info("END")
