import sys
from pathlib import Path

import astropy.units as u

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))  # benchmark/
from bench_common import MarkerDTime, resolve_data
from dysh.fits import GBTFITSLoad
from dysh.log import logger

logger.info("START")

dt = MarkerDTime(benchname="hi_survey")
tag = dt.tag

sdfits = GBTFITSLoad(resolve_data("/home/astro-util/HIsurvey/Session02", example="survey"))
tag("GBTFITSLoad")

tp0 = sdfits.gettp(scan=299, ifnum=0, plnum=0, fdnum=0).timeaverage()
tsys0 = tp0.meta["TSYS"]
tag("gettp(299,pl=0)+timeaverage")

tp1 = sdfits.gettp(scan=299, ifnum=0, plnum=1, fdnum=0).timeaverage()
tsys1 = tp1.meta["TSYS"]
tag("gettp(299,pl=1)+timeaverage")

sigref0a = sdfits.getsigref(scan=[296], ref=295, ifnum=0, fdnum=0, plnum=0, t_sys=tsys0).timeaverage()
tag("getsigref(296/295,pl=0)+timeaverage")

sigref0b = sdfits.getsigref(scan=[298], ref=297, ifnum=0, fdnum=0, plnum=0, t_sys=tsys0).timeaverage()
tag("getsigref(298/297,pl=0)+timeaverage")

sigref0 = sigref0a.average(sigref0b)
tag("average(sigref0a,sigref0b)")

sigref1a = sdfits.getsigref(scan=[296], ref=295, ifnum=0, fdnum=0, plnum=1, t_sys=tsys1).timeaverage()
tag("getsigref(296/295,pl=1)+timeaverage")

sigref1b = sdfits.getsigref(scan=[298], ref=297, ifnum=0, fdnum=0, plnum=1, t_sys=tsys1).timeaverage()
tag("getsigref(298/297,pl=1)+timeaverage")

sigref1 = sigref1a.average(sigref1b)
tag("average(sigref1a,sigref1b)")

sigref = sigref0.average(sigref1)
tag("average(sigref0,sigref1)")

sigref_smo = sigref.smooth(kernel="gauss", width=100, decimate=0)
tag("smooth(gauss,width=100,decimate=0)")

# GBTIDL's baseline windows land on these channel pairs for this benchmark dataset.
baseline_chans = [[47, 113], [128, 136], [184, 296]]
sigref_smo.baseline(1, model="poly", include=baseline_chans, remove=True)
tag("baseline(poly,deg=1)")

stats_b = sigref_smo[2000 * u.km / u.s : 2500 * u.km / u.s].stats()
tag("stats_b")

stats_r = sigref_smo[3500 * u.km / u.s : 4000 * u.km / u.s].stats()
tag("stats_r")

print(f"RMS_BLUE={stats_b['rms']:.6e}")
print(f"RMS_RED={stats_r['rms']:.6e}")

rms = (stats_b["rms"] + stats_r["rms"]) / 2.0

cog = sigref_smo.cog(bchan=100, echan=200)
tag("cog")

dt.report()

print(cog)

logger.info("END")
