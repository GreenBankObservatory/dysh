path = GETENV('DYSH_BENCH_DATA_PATH')
if path eq '' then path = '/home/astro-util/HIsurvey/Session02'
bench_t0 = systime(/sec)
stage_t0 = bench_t0
dirin, path
print, 'GBTIDL_BENCH_STAGE_MS[GBTFITSLoad]=', (systime(/sec) - stage_t0) * 1000.0
stage_t0 = systime(/sec)

jnk='jnk'
freeze

gettp,299,plnum=0, /quiet
tsys0=!g.s.tsys
print, format='("HI_DEBUG[tp0] TSYS=", E15.6)', tsys0[0]
print, 'GBTIDL_BENCH_STAGE_MS[gettp(299,pl=0)+timeaverage]=', (systime(/sec) - stage_t0) * 1000.0
stage_t0 = systime(/sec)

gettp,299,plnum=1, /quiet
tsys1=!g.s.tsys
print, format='("HI_DEBUG[tp1] TSYS=", E15.6)', tsys1[0]
print, 'GBTIDL_BENCH_STAGE_MS[gettp(299,pl=1)+timeaverage]=', (systime(/sec) - stage_t0) * 1000.0
stage_t0 = systime(/sec)

sclear

for i=295,297,2 do begin &$
  for p=0,1 do begin &$
    if (p eq 0) then tsys=tsys0[0] else tsys=tsys1[0] &$
    getsigref,i+1,i,plnum=p,tsys=tsys, /quiet &$
    if (i eq 295 and p eq 0) then print, 'GBTIDL_BENCH_STAGE_MS[getsigref(296/295,pl=0)+timeaverage]=', (systime(/sec) - stage_t0) * 1000.0 &$
    if (i eq 295 and p eq 0) then stage_t0 = systime(/sec) &$
    if (i eq 297 and p eq 0) then print, 'GBTIDL_BENCH_STAGE_MS[getsigref(298/297,pl=0)+timeaverage]=', (systime(/sec) - stage_t0) * 1000.0 &$
    if (i eq 297 and p eq 0) then stage_t0 = systime(/sec) &$
    if (i eq 295 and p eq 1) then print, 'GBTIDL_BENCH_STAGE_MS[getsigref(296/295,pl=1)+timeaverage]=', (systime(/sec) - stage_t0) * 1000.0 &$
    if (i eq 295 and p eq 1) then stage_t0 = systime(/sec) &$
    if (i eq 297 and p eq 1) then print, 'GBTIDL_BENCH_STAGE_MS[getsigref(298/297,pl=1)+timeaverage]=', (systime(/sec) - stage_t0) * 1000.0 &$
    if (i eq 297 and p eq 1) then stage_t0 = systime(/sec) &$
    accum &$
  endfor &$
endfor

ave

setxunit,'GHz'
setx,1.401,1.412

velo
show
stats,2000,2500, ret=dbg_blue
stats,3500,4000, ret=dbg_red
stats,2815,3142, ret=dbg_line
print, format='("HI_DEBUG[after_ave] BLUE_MEAN=", E15.6, " BLUE_RMS=", E15.6, " RED_MEAN=", E15.6, " RED_RMS=", E15.6, " LINE_MEAN=", E15.6, " LINE_RMS=", E15.6)', $
  dbg_blue.mean, dbg_blue.rms, dbg_red.mean, dbg_red.rms, dbg_line.mean, dbg_line.rms

setxunit,'GHz'
setx,1.401,1.412
gsmooth,100,/decimate
print, 'GBTIDL_BENCH_STAGE_MS[smooth(gauss,width=100,decimate=0)]=', (systime(/sec) - stage_t0) * 1000.0
stage_t0 = systime(/sec)

velo
show
stats,2000,2500, ret=dbg_blue
stats,3500,4000, ret=dbg_red
stats,2815,3142, ret=dbg_line
print, format='("HI_DEBUG[after_smooth] BLUE_MEAN=", E15.6, " BLUE_RMS=", E15.6, " RED_MEAN=", E15.6, " RED_RMS=", E15.6, " LINE_MEAN=", E15.6, " LINE_RMS=", E15.6)', $
  dbg_blue.mean, dbg_blue.rms, dbg_red.mean, dbg_red.rms, dbg_line.mean, dbg_line.rms

setxunit,'GHz'
setx,1.401,1.412
show
region=[1.402, 1.4045, 1.40506, 1.4054, 1.4072, 1.4115]
reg1=xtochan(region)
print, format='("HI_DEBUG[baseline_region_chans] ", 6(I0,1X))', reg1
Nregion,reg1
nfit,1
bshape
baseline
print, 'GBTIDL_BENCH_STAGE_MS[baseline(poly,deg=1)]=', (systime(/sec) - stage_t0) * 1000.0
stage_t0 = systime(/sec)

velo
show
stats,2000,2500, ret=dbg_blue
stats,3500,4000, ret=dbg_red
stats,2815,3142, ret=dbg_line
print, format='("HI_DEBUG[after_baseline] BLUE_MEAN=", E15.6, " BLUE_RMS=", E15.6, " RED_MEAN=", E15.6, " RED_RMS=", E15.6, " LINE_MEAN=", E15.6, " LINE_RMS=", E15.6)', $
  dbg_blue.mean, dbg_blue.rms, dbg_red.mean, dbg_red.rms, dbg_line.mean, dbg_line.rms
stats,2000,2500, ret=mystats
rms1 = mystats.rms
print, 'GBTIDL_BENCH_STAGE_MS[stats_b]=', (systime(/sec) - stage_t0) * 1000.0
stage_t0 = systime(/sec)
stats,3500,4000, ret=mystats
rms2 = mystats.rms
print, 'GBTIDL_BENCH_STAGE_MS[stats_r]=', (systime(/sec) - stage_t0) * 1000.0
stage_t0 = systime(/sec)

; Verification output
print, format='("RMS_BLUE=", E15.6)', rms1
print, format='("RMS_RED=", E15.6)', rms2

rms = (rms1 + rms2) / 2

gmeasure,1, 0.5,brange=2815,erange=3142,rms=rms
gmeasure,1, 0.2,brange=2815,erange=3142,rms=rms

print, 'GBTIDL_BENCH_SCRIPT_MS=', (systime(/sec) - bench_t0) * 1000.0
exit
