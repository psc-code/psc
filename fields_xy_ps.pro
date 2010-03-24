rimx=8.0
rimy=8.0

XSIZE=xdim/1000.0+rimx
YSIZE=ydim/1000.0+rimy

x0=500.0*rimx
y0=500.0*rimy


if (x2-x1 gt 1 and xdim gt 0) then begin
if (y2-y1 gt 1 and ydim gt 0) then begin

fxy=fltarr(x2-x1+1,y2-y1+1)

for zl=z1+0L,z2,stepz do begin

   fo=strmid(nmfc,0,strlen(nmfc)-13)$
   +'_t'+strmid(nmfc,strlen(nmfc)-13,5)$
   +'_z'+strtrim(string(zl),2)$
   +'.eps'

   print,fo

   SET_PLOT,'PS'
   DEVICE,/color,/encapsulated,file=fo,XSIZE=XSIZE,YSIZE=YSIZE
   erase, 255      

   for k=0L,y2-y1 do begin
      for j=0L,x2-x1 do begin
          fxy(j,k)=serie(x1+j,y1+k,zl)
      endfor
   endfor

   mini=min(fxy)
   maxi=max(fxy)
;   sfxy=smooth(fxy,5)
   sfxy=fxy


   CONTOUR,sfxy(0:x2-x1,0:y2-y1),x(x1:x2),y(y1:y2), $
   xs=1,ys=1,charsize=1.6, $
   POSITION=[x0,y0,x0+xdim,y0+ydim], $
   xtitle='!8x/!7l!8m', $
   xticks=2,$
;   xtickn=[strtrim(string(1.e6*x(x1)),2), $
;           strtrim(string(0.5e6*(x(x2)+x(x1))),2), $
;           strtrim(string(1.e6*x(x2)),2)], $
   xr=1.0e6*[x(x1),x(x2)], $
   ytitle='!8y/!7l!8m', $
   yticks=2,$
;   ytickn=[strtrim(string(1.e6*y(y1)),2), $
;          strtrim(string(0.5e6*(y(y2)+y(y1))),2), $
;          strtrim(string(1.e6*y(y2)),2)], $
   yr=1.0e6*[y(y1),y(y2)], $
   background=255,color=0, $ 
   /DEVICE,/NODATA

   ixy = BYTSCL(sfxy(0:x2-x1,0:y2-y1), MIN=rmin, MAX=rmax)
   TV, ixy, !X.WINDOW(0)+0.005*(!X.WINDOW(1)-!X.WINDOW(0)),!Y.WINDOW(0)+0.005*(!Y.WINDOW(1)-!Y.WINDOW(0)),$
       XSIZE=0.99*(!X.WINDOW(1) - !X.WINDOW(0)),$
       YSIZE=0.99*(!Y.WINDOW(1) - !Y.WINDOW(0)), /NORM

   if (rmax-rmin gt 0.0) then begin
      x11=(x0+xdim+500.0)/(xdim+2.0*x0)
      y11=y0/(ydim+2.0*y0)
      x22=(x0+xdim+800.0)/(xdim+2.0*x0)
      y22=(y0+ydim)/(ydim+2.0*y0)
      colorbar, division=2, NCOLORS=!D.Table_size, /right, /vertical, $
      format='(F9.3)',color=0, charsize=1.5, min=rmin, max=rmax, $
      /pscolor, POSITION=[x11,y11,x22,y22]
   endif

   ttt='!8t='+strtrim(string(t*1.0e15),2)+'fs' 
   pos='!8z='+strtrim(string(z(zl)*1.0e6),2)+'!7l!8m' 
   rmini='!8min='+strtrim(string(mini),2)
   rmaxi='!8max='+strtrim(string(maxi),2)
   xyouts,300,y0+ydim+2500.0,ttt,align=0.0,charsize=1.3,color=0,/DEVICE
   xyouts,300,y0+ydim+2000.0,pos,align=0.0,charsize=1.3,color=0,/DEVICE
   xyouts,300,y0+ydim+1500.0,rmaxi,align=0.0,charsize=1.3,color=0,/DEVICE
   xyouts,300,y0+ydim+1000.0,rmini,align=0.0,charsize=1.3,color=0,/DEVICE


;   vv=25
;   a1r=fltarr(vv)
;   a2r=fltarr(vv)
;   f1r=fltarr(vv,vv)
;   f2r=fltarr(vv,vv)

;   f1=fltarr(x2-x1+1,y2-y1+1)
;   f2=fltarr(x2-x1+1,y2-y1+1)
;   for k=0L,y2-y1 do begin
;      for j=0L,x2-x1 do begin
;         f1(j,k)=bxt(x1+j,y1+k,zl)
;         f2(j,k)=byt(x1+j,y1+k,zl)
;      endfor
;   endfor
;   sf1=smooth(f1,10)
;   sf2=smooth(f2,10)


;   for j=0L,vv-1 do begin
;      jj=round((y2-y1+1)*j/vv+0.5*(y2-y1+1)/vv)
;      a1r(j)=y(y1+jj)
;      for i=0L,vv-1 do begin
;         ii=round((x2-x1+1)*i/vv+0.5*(x2-x1+1)/vv)
;         a2r(i)=x(x1+ii)
;         f1r(i,j)=sf1(ii,jj)
;         f2r(i,j)=sf2(ii,jj)
;      endfor
;   endfor


;   f1r(0,*)=0.0
;   f2r(0,*)=0.0
;   f1r(1,*)=0.0
;   f2r(1,*)=0.0
;   f1r(vv-1,*)=0.0
;   f2r(vv-1,*)=0.0
;   f1r(vv-2,*)=0.0
;   f2r(vv-2,*)=0.0
;   f1r(*,0)=0.0
;   f2r(*,0)=0.0
;   f1r(*,1)=0.0
;   f2r(*,1)=0.0
;   f1r(*,vv-1)=0.0
;   f2r(*,vv-1)=0.0
;   f1r(*,vv-2)=0.0
;   f2r(*,vv-2)=0.0


;   velovect,f1r,f2r,a1r,a2r,xs=5,ys=5,len=2.0,$
;   charsize=2.2,color=300,/noerase,thick=3.0,$   
;   POSITION=[x0,y0,x0+xdim,y0+ydim],/DEVICE


   DEVICE, /CLOSE
   SET_PLOT,'X'

endfor

endif
endif

