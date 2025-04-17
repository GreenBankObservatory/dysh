To generate the file ``
using `GBTIDL`

```
vbanks = ["A", "B", "C", "D"]
ifnums = [2, 1, 0, 3]
for i=0,3,1 do begin & print,vbanks[i] & filein,"../AGBT18B_354_03.raw.vegas/AGBT18B_354_03.raw.vegas."+vbanks[i]+".fits" & fileout,"gettp_vbank_"+vbanks[i]+".fits" & for scan=6,7,1 do begin & for plnum=0,1,1 do begin & for cal_state=0,1,1 do begin gettp,scan,ifnum=ifnums[i],plnum=plnum,intnum=0,cal_state=cal_state & keep & endfor & endfor & endfor & endfor
```
