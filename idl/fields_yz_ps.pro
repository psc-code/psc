rimx=8.0
rimy=8.0

XSIZE=xdim/1000.0+rimx
YSIZE=ydim/1000.0+rimy

x0=500.0*rimx
y0=500.0*rimy


if (y2-y1 gt 1 and xdim gt 0) then begin
if (z2-z1 gt 1 and ydim gt 0) then begin

fyz=fltarr(y2-y1+1,z2-z1+1)

for xl=x1+0L,x2,stepx do begin

   fo=strmid(nmfc,0,strlen(nmfc)-13)$
   +'_t'+strmid(nmfc,strlen(nmfc)-13,5)$
   +'_x'+strtrim(string(xl),2)$
   +'.eps'

   print,fo

   SET_PLOT,'PS'
   DEVICE,/color,/encapsulated,file=fo,XSIZE=XSIZE,YSIZE=YSIZE
   erase, 255      

   for k=0L,z2-z1 do begin
      for j=0L,y2-y1 do begin
          fyz(j,k)=serie(xl,y1+j,z1+k)
      endfor
   endfor

   mini=min(fyz)
   maxi=max(fyz)
;   sfyz=smooth(fyz,5)
   sfyz=fyz

   CONTOUR,sfyz(0:y2-y1,0:z2-z1),y(y1:y2),z(z1:z2), $
   xs=1,ys=1,charsize=1.6, $
   POSITION=[x0,y0,x0+xdim,y0+ydim], $
   xtitle='!8y/!7l!8m', $
   xticks=2,$
;   xtickn=[strtrim(string(1.e6*y(y1)),2), $
;           strtrim(string(0.5e6*(y(y2)+y(y1))),2), $
;           strtrim(string(1.e6*y(y2)),2)], $
   xr=1.0e6*[y(y1),y(y2)], $
   ytitle='!8z/!7l!8m', $
   yticks=2,$
;   ytickn=[strtrim(string(1.e6*z(z1)),2),
;           strtrim(string(0.5e6*(z(z2)+z(z1))),2),
;           strtrim(string(1.e6*z(z2)),2)], $
   yr=1.0e6*[z(z1),z(z2)], $
   background=255,color=0, $ 
   /DEVICE,/NODATA

   iyz = BYTSCL(sfyz(0:y2-y1,0:z2-z1), MIN=rmin, MAX=rmax)
   TV, iyz, !X.WINDOW(0)+0.005*(!X.WINDOW(1)-!X.WINDOW(0)),!Y.WINDOW(0)+0.005*(!Y.WINDOW(1)-!Y.WINDOW(0)),$
       XSIZE=0.99*(!X.WINDOW(1) - !X.WINDOW(0)),$
       YSIZE=0.99*(!Y.WINDOW(1) - !Y.WINDOW(0)), /NORM

   if (rmax-rmin gt 0.0) then begin
      x11=(x0+xdim+500.0)/(xdim+2.0*x0)
      y11=y0/(ydim+2.0*y0)
      x22=(x0+xdim+800.0)/(xdim+2.0*x0)
      y22=(y0+ydim)/(ydim+2.0*y0)
      colorbar, division=2, NCOLORS=!D.Table_size, /right, /vertical, $
      format='(F9.4)',color=0, charsize=1.5, min=rmin, max=rmax, $
      /pscolor, POSITION=[x11,y11,x22,y22]
   endif

   ttt='!8t='+strtrim(string(t*1.0e15),2)+'fs' 
   pos='!8x='+strtrim(string(x(xl)*1.0e6),2)+'!7l!8m' 
   rmini='!8min='+strtrim(string(mini),2)
   rmaxi='!8max='+strtrim(string(maxi),2)
   xyouts,300,y0+ydim+2500.0,ttt,align=0.0,charsize=1.3,color=0,/DEVICE
   xyouts,300,y0+ydim+2000.0,pos,align=0.0,charsize=1.3,color=0,/DEVICE
   xyouts,300,y0+ydim+1500.0,rmaxi,align=0.0,charsize=1.3,color=0,/DEVICE
   xyouts,300,y0+ydim+1000.0,rmini,align=0.0,charsize=1.3,color=0,/DEVICE

   DEVICE, /CLOSE
   SET_PLOT,'X'

endfor

endif
endif
