rimx=8.0
rimy=8.0

XSIZE=xdim/1000.0+rimx
YSIZE=ydim/1000.0+rimy

x0=500.0*rimx
y0=500.0*rimy

ddxp=stepx*dx*ld
dyp=stepy*dy*ld
dzp=stepz*dz*ld

xlen=size(x,/DIMENSIONS)
ylen=size(y,/DIMENSIONS)
zlen=size(z,/DIMENSIONS)


charge=fltarr(10)
charge(0)=qnp(0)
qnum=0L
for i=0L,nqcount-1 do begin
   exist=0.0
   for j=0L,qnum do begin
      if (charge(j) eq qnp(i)) then begin
         exist=1.0
      endif
   endfor
   if (exist eq 0.0) then begin
      qnum=qnum+1
      charge(qnum)=qnp(i)
   endif
endfor


if (x2 gt 1) then begin
if (z2 gt 1) then begin

for yl=y1+0L,y2,stepy do begin
for ql=0L,qnum do begin

   for i=0L,nqcount-1 do begin
      if (z(z1)-dz*ld le zp(i)) and (zp(i) le z(z1)+dzp) then begin
      if (y(yl)-dy*ld le yp(i)) and (yp(i) le y(yl)+dyp) then begin
      if (x(x1)-dx*ld le xp(i)) and (xp(i) le x(x1)+dxp) then begin
      if (pxmin le pxp(i)) and (pxp(i) le pxmax) then begin
      if (pymin le pyp(i)) and (pyp(i) le pymax) then begin
      if (pzmin le pzp(i)) and (pzp(i) le pzmax) then begin
         if (qnp(i) eq charge(ql)) then begin
            pnum=pnum+1
         endif
      endif
      endif
      endif
      endif
      endif
      endif
   endfor

   if (pnum gt 0) then begin

      xpl=fltarr(pnum)
      ypl=fltarr(pnum)
      zpl=fltarr(pnum)
      pxpl=fltarr(pnum)
      pypl=fltarr(pnum)
      pzpl=fltarr(pnum)
      qnpl=fltarr(pnum)
      mnpl=fltarr(pnum)
      lnpl=fltarr(pnum)
      wnpl=fltarr(pnum)

      pnum=0L
      for i=0L,nqcount-1 do begin
         if (z(z1)-dz*ld le zp(i)) and (zp(i) le z(z1)+dzp) then begin
         if (y(yl)-dy*ld le yp(i)) and (yp(i) le y(yl)+dyp) then begin
         if (x(x1)-dx*ld le xp(i)) and (xp(i) le x(x1)+dxp) then begin
         if (pxmin le pxp(i)) and (pxp(i) le pxmax) then begin
         if (pymin le pyp(i)) and (pyp(i) le pymax) then begin
         if (pzmin le pzp(i)) and (pzp(i) le pzmax) then begin
            if (qnp(i) eq charge(ql)) then begin
               xpl(pnum)=1.0e6*xp(i)
               ypl(pnum)=1.0e6*yp(i)
               zpl(pnum)=1.0e6*zp(i)
               pxpl(pnum)=pxp(i)
               pypl(pnum)=pyp(i)
               pzpl(pnum)=pzp(i)
               qnpl(pnum)=qnp(i)
               mnpl(pnum)=mnp(i)
               lnpl(pnum)=lnp(i)
               wnpl(pnum)=wnp(i)
               pnum=pnum+1
            endif
         endif
         endif
         endif
         endif
         endif
         endif
      endfor


      fo=strmid(nmfc,0,strlen(nmfc)-13)+'_xz'$
      +'_q'+strtrim(string(charge(ql)),2)$
      +'_t'+strmid(nmfc,strlen(nmfc)-13,5)$
      +'_y'+strtrim(string(yl),2)$
      +'.eps'

      print,fo

      SET_PLOT,'PS'
      DEVICE,/color,/encapsulated,file=fo,XSIZE=XSIZE,YSIZE=YSIZE
      erase, 255      


      plot, xpl, zpl, xs=1, ys=1, charsize=1.6, psym=3, $
      POSITION=[x0,y0,x0+xdim,y0+ydim],/DEVICE, $ 
      xtitle='!8x/!7l!8m',$
      xticks=2,$
;      xtickn=[strtrim(string(1.e6*x(x1)),2), $
;              strtrim(string(0.5e6*(x(x2)+x(x1))),2), $
;              strtrim(string(1.e6*x(x2)),2)], $
      xr=1.0e6*[x(x1),x(x2)], $
      ytitle='!8z/!7l!8m',$
      yticks=2,$
;      ytickn=[strtrim(string(1.e6*z(z1)),2), $
;              strtrim(string(0.5e6*(z(z2)+z(z1))),2), $
;              strtrim(string(1.e6*z(z2)),2)], $
      yr=1.0e6*[z(z1),z(z2)], $
      color=0,background=255 

      number=0.0
      for i=0L,pnum-1 do begin      
         number=number+wnpl(i)
      endfor
      number=n0*cori*nistep*step*dx*dy*dz*ld*ld*ld*number

      ttt='!8t='+strtrim(string(t*1.0e15),2)+'fs'
      posy='!8y='+strtrim(string((y(yl)-dy*ld)*1.0e6),2)+'-'+strtrim(string((y(yl)+dyp)*1.0e6),2)+'!7l!8m'
      momx='!8px='+strtrim(string(pxmin),2)+'-'+strtrim(string(pxmax),2)+'!8mc'
      momy='!8py='+strtrim(string(pymin),2)+'-'+strtrim(string(pymax),2)+'!8mc'
      momz='!8pz='+strtrim(string(pzmin),2)+'-'+strtrim(string(pzmax),2)+'!8mc'
      cha='!8q='+strtrim(string(charge(ql)),2)+'!8e'
      number='!8N='+strtrim(string(number),2)
      xyouts,300,y0+ydim+2700.0,ttt,align=0.0,charsize=1.0,color=0,/DEVICE
      xyouts,300,y0+ydim+2400.0,posy,align=0.0,charsize=1.0,color=0,/DEVICE
      xyouts,300,y0+ydim+2100.0,momx,align=0.0,charsize=1.0,color=0,/DEVICE
      xyouts,300,y0+ydim+1800.0,momy,align=0.0,charsize=1.0,color=0,/DEVICE
      xyouts,300,y0+ydim+1500.0,momz,align=0.0,charsize=1.0,color=0,/DEVICE
      xyouts,300,y0+ydim+1200.0,cha,align=0.0,charsize=1.0,color=0,/DEVICE
      xyouts,300,y0+ydim+900.0,number,align=0.0,charsize=1.0,color=0,/DEVICE

      DEVICE, /CLOSE
      SET_PLOT,'X'

   endif
endfor
endfor

endif
endif
