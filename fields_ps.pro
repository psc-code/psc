rimx=8.0
rimy=8.0

XSIZE=xdim/1000.0+rimx
YSIZE=ydim/1000.0+rimy

x0=500.0*rimx
y0=500.0*rimy


if (xdim gt 0) then begin
if (ydim gt 0) then begin


   fs=fltarr(z2-z1+1)
   for k=0L,z2-z1 do begin
      for j=0L,y2-y1 do begin
         for i=0L,x2-x1 do begin
            fs(k)=fs(k)+serie(x1+i,y1+j,z1+k)/(x2-x1+1)/(y2-y1+1)
         endfor
      endfor
   endfor
   f=smooth(fs,2)
   f=fs

   mini=min(f)
   maxi=max(f)
   rmin=mini
   rmax=maxi
   if (rmin eq 0.0) then begin
      rmin=-0.01
   endif
   if (rmax eq 0.0) then begin
      rmax=+0.01
   endif


   fo=strmid(nmfc,0,strlen(nmfc)-13)$
   +'_t'+strmid(nmfc,strlen(nmfc)-13,5)$
   +'.eps'

   print,fo

   SET_PLOT,'PS'
   DEVICE,/color,/encapsulated,file=fo,XSIZE=XSIZE,YSIZE=YSIZE
   erase, 255      


   plot, 1.0e6*z(z1:z2), f(0:z2-z1), xs=1, ys=1, charsize=1.5, psym=0, $
   POSITION=[x0,y0,x0+xdim,y0+ydim],/DEVICE, $ 
   xtitle='!8z/!7l!8m',$
   xticks=2,$
   xr=1.0e6*[z(z1),z(z2)], $
   ytitle='!8f/f!D0',$
   yticks=2,$
   yr=[rmin,rmax], $
   color=0,background=255 

   ttt='!8t='+strtrim(string(t*1.0e15),2)+'fs'
   posx='!8x='+strtrim(string(x(x1)*1.0e6),2)+'-'+strtrim(string(x(x2)*1.0e6),2)+'!7l!8m'
   posy='!8y='+strtrim(string(y(y1)*1.0e6),2)+'-'+strtrim(string(y(y2)*1.0e6),2)+'!7l!8m'
   xyouts,300,y0+ydim+2700.0,ttt,align=0.0,charsize=1.0,color=0,/DEVICE
   xyouts,300,y0+ydim+2400.0,posx,align=0.0,charsize=1.0,color=0,/DEVICE
   xyouts,300,y0+ydim+2100.0,posy,align=0.0,charsize=1.0,color=0,/DEVICE

   DEVICE, /CLOSE
   set_plot,'X'

endif
endif
