path = GETENV('DYSH_BENCH_DATA_PATH')
if path eq '' then path = '/home/scratch/ajschmie/training/dysh/datasets/argus/TGBT22A_603_05_vanecal.raw.vegas'
bench_t0 = systime(/sec)
stage_t0 = systime(/sec)
dirin, path
print, 'GBTIDL_BENCH_STAGE_MS[GBTFITSLoad]=', (systime(/sec) - stage_t0) * 1000.0
.compile scripts/argus_vanecal/getatmos.pro
.compile scripts/argus_vanecal/vanecal.pro
stage_t0 = systime(/sec)
vanecal, 10
print, 'GBTIDL_BENCH_STAGE_MS[vanecal_total]=', (systime(/sec) - stage_t0) * 1000.0
print, 'GBTIDL_BENCH_SCRIPT_MS=', (systime(/sec) - bench_t0) * 1000.0
exit
