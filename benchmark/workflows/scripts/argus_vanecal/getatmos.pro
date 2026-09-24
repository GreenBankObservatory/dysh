pro getatmos,mjd=time,freq=freq,debug=debug,OUTput

;;Procedure to return zenith opacity and effective Tatm
;;from RMaddale
;; modifications by A.Schmiedeke

if (n_elements(debug) eq 0) then debug=0

my1='/users/rmaddale/bin/getForecastValues -freqList'
my2=string(freq)
my3=' -timeList '
my4=string(time)

mystr0=my1+my2+my3+my4

mystr1=mystr0+' -type Opacity'
mystr2=mystr0+' -type Tatm'

spawn,mystr1,result1
spawn,mystr2,result2

if (debug gt 0) then begin
    print,'(zenith)',result1
    print,result2
endif


tmp=strsplit(result1,/extract)
ztau=float(tmp[2])

tmp=strsplit(result2,/extract)
tatm=float(tmp[2])

OUTput=[ztau, tatm]
return
end
