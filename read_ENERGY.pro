; ENERGY RELAXATION

fi=FINDFILE('00000energy.data')
ifi=size(fi,/N_ELEMENTS)

print,fi(ifi-1)

openr,1,fi(ifi-1)
count=0L
while not(eof(1)) do begin 
   readf,1,time
   readf,1,te
   readf,1,ti
   readf,1,Ediff
   readf,1,Etot
   readf,1,xne
   readf,1,xni
   print,time,te,ti,Ediff,Etot,xne,xni
   count=count+1
endwhile
print,'Number of list elements:',count
close,1

atime=fltarr(count)
ate=fltarr(count)
ati=fltarr(count)
aEdiff=fltarr(count)
aEtot=fltarr(count)
axne=fltarr(count)
axni=fltarr(count)
aplot=fltarr(count)

openr,1,fi(ifi-1)
for i=0L,count-1 do begin 
   readf,1,time
   readf,1,te
   readf,1,ti
   readf,1,Ediff
   readf,1,Etot
   readf,1,xne
   readf,1,xni
   atime(i)=time
   ate(i)=te
   ati(i)=ti
   aEdiff(i)=Ediff
   aEtot(i)=Etot
   axne(i)=xne
   axni(i)=xni
endfor
close,1

for i=0L,count-1 do begin 
   aplot(i)=(ate(i)-ati(i))/(ate(0)-ati(0))
endfor

ydim=10000.0
xdim=10000.0

rimx=8.0
rimy=6.0

XSIZE=xdim/1000.0+rimx
YSIZE=ydim/1000.0+rimy

x0=500.0*rimx
y0=500.0*rimy


if (xdim gt 0) then begin
if (ydim gt 0) then begin

   fo='energy.eps'

   print,fo


   SET_PLOT,'PS'
   DEVICE,/color,/encapsulated,file=fo,XSIZE=XSIZE,YSIZE=YSIZE
   erase, 255      

   plot, atime, aplot, xs=1, ys=1, charsize=1.5, psym=3, /ylog, $
   POSITION=[x0,y0,x0+xdim,y0+ydim],/DEVICE, $ 
   xtitle='!8t/!8fs',$
   xticks=2,$
   xr=[atime(0),atime(count-1)], $
   ytitle='!7e!8/!8mJ',$
   yticks=2,$
   yr=[0.001,max(aplot)], $
   color=0,background=255 

;   axis, yaxis=1, charsize=1.5, ys=1, $
;   ytitle='!8n!Dp!N/!810!U8',$
;   yr=1.0e-8*[0.9*min(numb(ml,*)),1.1*max(numb(ml,*))],$
;   yticks=2, /save
;   oplot, time(ml,*), 1.0e-8*numb(ml,*), linestyle=3

endif
endif

DEVICE, /CLOSE
set_plot,'X'
end
