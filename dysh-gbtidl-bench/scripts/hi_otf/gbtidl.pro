filename = GETENV('DYSH_BENCH_DATA_PATH')
if filename eq '' then filename = '/home/scratch/dfrayer/DATAdemo/TGBT17A_506_11.raw.vegas'
bench_t0 = systime(/sec)
filein, filename
ref = 27
out_path = GETENV('DYSH_BENCH_OUT_PATH')
if out_path eq '' then out_path = 'gbtidl.fits'
fileout, out_path, /new
for s=14,26,1 do begin &$
    getsigref,s,ref,/avgref,/keepints &$
    keep &$
endfor
print, 'GBTIDL_BENCH_SCRIPT_MS=', (systime(/sec) - bench_t0) * 1000.0
exit
