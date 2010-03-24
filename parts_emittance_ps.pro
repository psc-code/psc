dxp=stepx*dx*ld
dyp=stepy*dy*ld
dzp=stepz*dz*ld
xlen=size(x,/DIMENSIONS)
ylen=size(y,/DIMENSIONS)
zlen=size(z,/DIMENSIONS)


mass=fltarr(10)
mass(0)=mnp(0)
mnum=0L
for i=0L,nqcount-1 do begin
   exist=0.0
   for j=0L,mnum do begin
      if (mass(j) eq mnp(i)) then begin
         exist=1.0
      endif
   endfor
   if (exist eq 0.0) then begin
      mnum=mnum+1
      mass(mnum)=mnp(i)
   endif
endfor

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


if (xlen(0)-1 gt 1) then begin
if (ylen(0)-1 gt 1) then begin


for ml=0L,mnum do begin
for ql=0L,qnum do begin
for pzl=1,pznum do begin
;for pzl=1,1 do begin


   pzmin=(1.0*pzl)/mult
   pzmax=(1.0*pzl+1)/mult


ydim=8000.00
xdim=ydim*(i1x-i1n+1)/(i3x-i3n+1)

rimx=6.0
rimy=6.0

XSIZE=xdim/1000.0+rimx
YSIZE=ydim/1000.0+rimy

x0=500.0*rimx
y0=500.0*rimy


for zl=z1+0L,z2,stepz do begin
for yl=y1+0L,y2,stepy do begin


   pnum=0L
   for i=0L,nqcount-1 do begin
      if (z(zl) le zp(i)) and (zp(i) le z(zl)+dzp) then begin
      if (y(yl) le yp(i)) and (yp(i) le y(yl)+dyp) then begin
      if (mnp(i) eq mass(ml)) then begin
      if (qnp(i) eq charge(ql)) then begin
      if (pzmin le pzp(i)) and (pzp(i) le pzmax) then begin
         pnum=pnum+1
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

      pnum=0L
      for i=0L,nqcount-1 do begin
         if (z(zl) le zp(i)) and (zp(i) le z(zl)+dzp) then begin
         if (y(yl) le yp(i)) and (yp(i) le y(yl)+dyp) then begin
         if (mnp(i) eq mass(ml)) then begin
         if (qnp(i) eq charge(ql)) then begin
         if (pzmin le pzp(i)) and (pzp(i) le pzmax) then begin
            xpl(pnum)=1.0e6*xp(i)
            ypl(pnum)=1.0e6*yp(i)
            zpl(pnum)=1.0e6*zp(i)
;            pxpl(pnum)=pxp(i)/pzp(i)
;            pypl(pnum)=pyp(i)/pzp(i)
            pxpl(pnum)=pxp(i)
            pypl(pnum)=pyp(i)
            pzpl(pnum)=pzp(i)
            qnpl(pnum)=qnp(i)
            mnpl(pnum)=mnp(i)
            lnpl(pnum)=lnp(i)
            pnum=pnum+1
         endif
         endif
         endif
         endif
         endif
      endfor


      nx=120
      ny=120
      xs=1.0e6*x1*dx*ld
      xu=1.0e6*x2*dx*ld
      ys=pxmin
      yu=pxmax
      dpx=(xu-xs)/nx
      dpy=(yu-ys)/ny
      nf=fltarr(nx+2,ny+2)

      p=0L
      for l=0L,pnum-1 do begin
 
         xq=xpl(l)
         yq=pxpl(l)
 
         if (xs lt xq) and (xq lt xu) then begin
         if (ys lt yq) and (yq lt yu) then begin
 
         k1=floor((xq-xs)/dpx)
         k2=floor((yq-ys)/dpy)
 
         xq=(xq-xs)/dpx-k1
         yq=(yq-ys)/dpy-k2
 
         h1=(1.0-xq)*(1.0-yq)
         h2=xq*(1.0-yq)
         h3=(1.0-xq)*yq
         h4=xq*yq
 
         nf(k1,k2)=nf(k1,k2)+h1
         nf(k1+1,k2)=nf(k1+1,k2)+h2
         nf(k1,k2+1)=nf(k1,k2+1)+h3
         nf(k1+1,k2+1)=nf(k1+1,k2+1)+h4

         endif
         endif

      endfor


      fo=strmid(nmfc,0,strlen(nmfc)-13)+'_xxp-pz-scale'$
      +'_pz'+strtrim(string(pzmin),2)+'-'+strtrim(string(pzmax),2)$
      +'_m'+strtrim(string(mass(ml)),2)$
      +'_q'+strtrim(string(charge(ql)),2)$
      +'_t'+strmid(nmfc,strlen(nmfc)-13,5)$
      +'_y'+strtrim(string(yl),2)$
      +'_z'+strtrim(string(zl),2)$
      +'.eps'


      print,fo


      SET_PLOT,'PS'
      DEVICE,/color,/encapsulated,file=fo,XSIZE=XSIZE,YSIZE=YSIZE
      erase, 255      


      sfxy=smooth(nf,5)
      rmin=min(sfxy)
      rmax=max(sfxy)
 
      CONTOUR,sfxy(0:nx,0:ny), $
      xs=1,ys=1,charsize=1.6, $
      POSITION=[x0,y0,x0+xdim,y0+ydim], $
      xtitle='!8x/!7l!8m', $
      xticks=2,$
      xr=[xs,xu], $
      ytitle='!8p!Dx!N/!8p!Dz',$
      yticks=2,$
      yr=[ys,yu], $
      background=255,color=0, $
      /DEVICE,/NODATA

      ixy = BYTSCL(sfxy(0:nx,0:ny), MIN=rmin, MAX=rmax)
      TV, ixy, !X.WINDOW(0),!Y.WINDOW(0),$
      XSIZE=!X.WINDOW(1) - !X.WINDOW(0),$
      YSIZE=!Y.WINDOW(1) - !Y.WINDOW(0), /NORM
 
      if (rmax-rmin gt 0.0) then begin
         x11=(x0+xdim+500.0)/(xdim+2.0*x0)
         y11=y0/(ydim+2.0*y0)
         x22=(x0+xdim+800.0)/(xdim+2.0*x0)
         y22=(y0+ydim)/(ydim+2.0*y0)
         colorbar, division=2, NCOLORS=!D.Table_size, /right, /vertical, $
         format='(F8.2)',color=0, charsize=1.5, min=rmin, max=rmax, $
         /pscolor, POSITION=[x11,y11,x22,y22]
      endif

            zmin=min(zpl)
            zmax=max(zpl)
            number=n0*cori*dx*dy*dz*ld*ld*ld*pnum
            ttt='!8t='+strtrim(string(t*1.0e15),2)+'fs'
            posy='!8y='+strtrim(string(y(yl)*1.0e6),2)+'-'+strtrim(string((y(yl)+dyp)*1.0e6),2)+'!7l!8m' 
            posz='!8z='+strtrim(string(z(zl)*1.0e6),2)+'-'+strtrim(string((z(zl)+dzp)*1.0e6),2)+'!7l!8m' 
            posz='!8z='+strtrim(string(zmin),2)+'-'+strtrim(string(zmax),2)+'!7l!8m' 
            number='!8N='+strtrim(string(number),2)
            xyouts,300,y0+ydim+2500.0,ttt,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+2200.0,posy,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+1900.0,posz,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+1600.0,number,align=0.0,charsize=1.0,color=0,/DEVICE

      DEVICE, /CLOSE
      SET_PLOT,'X'


      fo=strmid(nmfc,0,strlen(nmfc)-13)+'_xxp-pz'$
      +'_pz'+strtrim(string(pzmin),2)+'-'+strtrim(string(pzmax),2)$
      +'_m'+strtrim(string(mass(ml)),2)$
      +'_q'+strtrim(string(charge(ql)),2)$
      +'_t'+strmid(nmfc,strlen(nmfc)-13,5)$
      +'_y'+strtrim(string(yl),2)$
      +'_z'+strtrim(string(zl),2)$
      +'.eps'

      print,fo


      SET_PLOT,'PS'
      DEVICE,/color,/encapsulated,file=fo,XSIZE=XSIZE,YSIZE=YSIZE
      erase, 255      


      plot, xpl, pxpl, xs=1, ys=1, charsize=1.6, psym=3, $
      POSITION=[x0,y0,x0+xdim,y0+ydim],/DEVICE, $ 
      xtitle='!8x/!7l!8m',$
      xticks=2,$
;      xtickn=[strtrim(string(1.e6*x(x1)),2), $
;              strtrim(string(0.5e6*(x(x2)+x(x1))),2), $
;              strtrim(string(1.e6*x(x2)),2)], $
      xr=1.0e6*[x(x1),x(x2)], $
      ytitle='!8p!Dx!N/!8p!Dz',$
      yticks=2,$
;      ytickn=[strtrim(string(1.e6*y(y1)),2), $
;              strtrim(string(0.5e6*(y(y2)+y(y1))),2), $
;              strtrim(string(1.e6*y(y2)),2)], $
      yr=[pxmin,pxmax], $
      color=0,background=255 

            zmin=min(zpl)
            zmax=max(zpl)
            number=n0*cori*dx*dy*dz*ld*ld*ld*pnum
            ttt='!8t='+strtrim(string(t*1.0e15),2)+'fs'
            posy='!8y='+strtrim(string(y(yl)*1.0e6),2)+'-'+strtrim(string((y(yl)+dyp)*1.0e6),2)+'!7l!8m' 
            posz='!8z='+strtrim(string(z(zl)*1.0e6),2)+'-'+strtrim(string((z(zl)+dzp)*1.0e6),2)+'!7l!8m' 
            posz='!8z='+strtrim(string(zmin),2)+'-'+strtrim(string(zmax),2)+'!7l!8m' 
            number='!8N='+strtrim(string(number),2)
            xyouts,300,y0+ydim+2500.0,ttt,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+2200.0,posy,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+1900.0,posz,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+1600.0,number,align=0.0,charsize=1.0,color=0,/DEVICE

      DEVICE, /CLOSE
      SET_PLOT,'X'

   endif

endfor
endfor


ydim=8000.00
xdim=ydim*(i2x-i2n+1)/(i3x-i3n+1)

rimx=6.0
rimy=6.0

XSIZE=xdim/1000.0+rimx
YSIZE=ydim/1000.0+rimy

x0=500.0*rimx
y0=500.0*rimy


for zl=z1+0L,z2,stepz do begin
for xl=x1+0L,x2,stepx do begin


   pnum=0L
   for i=0L,nqcount-1 do begin
      if (z(zl) le zp(i)) and (zp(i) le z(zl)+dzp) then begin
      if (x(xl) le xp(i)) and (xp(i) le x(xl)+dxp) then begin
      if (mnp(i) eq mass(ml)) then begin
      if (qnp(i) eq charge(ql)) then begin
      if (pzmin le pzp(i)) and (pzp(i) le pzmax) then begin
         pnum=pnum+1
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

      pnum=0L
      for i=0L,nqcount-1 do begin
         if (z(zl) le zp(i)) and (zp(i) le z(zl)+dzp) then begin
         if (x(xl) le xp(i)) and (xp(i) le x(xl)+dxp) then begin
         if (mnp(i) eq mass(ml)) then begin
         if (qnp(i) eq charge(ql)) then begin
         if (pzmin le pzp(i)) and (pzp(i) le pzmax) then begin
            xpl(pnum)=1.0e6*xp(i)
            ypl(pnum)=1.0e6*yp(i)
            zpl(pnum)=1.0e6*zp(i)
            pxpl(pnum)=pxp(i)/pzp(i)
            pypl(pnum)=pyp(i)/pzp(i)
            pxpl(pnum)=pxp(i)
            pypl(pnum)=pyp(i)
            pzpl(pnum)=pzp(i)
            qnpl(pnum)=qnp(i)
            mnpl(pnum)=mnp(i)
            lnpl(pnum)=lnp(i)
            pnum=pnum+1
         endif
         endif
         endif
         endif
         endif
      endfor


      nx=120
      ny=120
      xs=1.0e6*y1*dy*ld
      xu=1.0e6*y2*dy*ld
      ys=pymin
      yu=pymax
      dpx=(xu-xs)/nx
      dpy=(yu-ys)/ny
      nf=fltarr(nx+2,ny+2)

      p=0L
      for l=0L,pnum-1 do begin
 
         xq=ypl(l)
         yq=pypl(l)
 
         if (xs lt xq) and (xq lt xu) then begin
         if (ys lt yq) and (yq lt yu) then begin
 
         k1=floor((xq-xs)/dpx)
         k2=floor((yq-ys)/dpy)
 
         xq=(xq-xs)/dpx-k1
         yq=(yq-ys)/dpy-k2
 
         h1=(1.0-xq)*(1.0-yq)
         h2=xq*(1.0-yq)
         h3=(1.0-xq)*yq
         h4=xq*yq
 
         nf(k1,k2)=nf(k1,k2)+h1
         nf(k1+1,k2)=nf(k1+1,k2)+h2
         nf(k1,k2+1)=nf(k1,k2+1)+h3
         nf(k1+1,k2+1)=nf(k1+1,k2+1)+h4

         endif
         endif

      endfor


      fo=strmid(nmfc,0,strlen(nmfc)-13)+'_yyp-pz-scale'$
      +'_pz'+strtrim(string(pzmin),2)+'-'+strtrim(string(pzmax),2)$
      +'_m'+strtrim(string(mass(ml)),2)$
      +'_q'+strtrim(string(charge(ql)),2)$
      +'_t'+strmid(nmfc,strlen(nmfc)-13,5)$
      +'_x'+strtrim(string(xl),2)$
      +'_z'+strtrim(string(zl),2)$
      +'.eps'


      print,fo


      SET_PLOT,'PS'
      DEVICE,/color,/encapsulated,file=fo,XSIZE=XSIZE,YSIZE=YSIZE
      erase, 255      


      sfxy=smooth(nf,6)
      rmin=min(sfxy)
      rmax=max(sfxy)

      CONTOUR,sfxy(0:nx,0:ny), $
      xs=1,ys=1,charsize=1.6, $
      POSITION=[x0,y0,x0+xdim,y0+ydim], $
      xtitle='!8y/!7l!8m', $
      xticks=2,$
      xr=[xs,xu], $
      ytitle='!8p!Dy!N/!8p!Dz',$
      yticks=2,$
      yr=[ys,yu], $
      background=255,color=0, $
      /DEVICE,/NODATA

      ixy = BYTSCL(sfxy(0:nx,0:ny), MIN=rmin, MAX=rmax)
      TV, ixy, !X.WINDOW(0),!Y.WINDOW(0),$
      XSIZE=!X.WINDOW(1) - !X.WINDOW(0),$
      YSIZE=!Y.WINDOW(1) - !Y.WINDOW(0), /NORM

      line=fltarr(nx+2)
      fline=fltarr(nx+2)
      for i=0L,nx+1 do begin
         line(i)=1.0e6*i*dy*ld
;         fline(i)=0.07*sin(8.0*6.2831853072*line(i)/20.0)+0.025*(line(i)-20.0)
         fline(i)=0.0125*(line(i)-20.0)
      endfor
;      oplot, line[402:749], fline[402:749], color=300.0, thick=2.0

      if (rmax-rmin gt 0.0) then begin
         x11=(x0+xdim+500.0)/(xdim+2.0*x0)
         y11=y0/(ydim+2.0*y0)
         x22=(x0+xdim+800.0)/(xdim+2.0*x0)
         y22=(y0+ydim)/(ydim+2.0*y0)
         colorbar, division=2, NCOLORS=!D.Table_size, /right, /vertical, $
         format='(F8.2)',color=0, charsize=1.5, min=rmin, max=rmax, $
         /pscolor, POSITION=[x11,y11,x22,y22]
      endif

            zmin=min(zpl)
            zmax=max(zpl)
            number=n0*cori*dx*dy*dz*ld*ld*ld*pnum
            ttt='!8t='+strtrim(string(t*1.0e15),2)+'fs'
            posx='!8x='+strtrim(string(x(xl)*1.0e6),2)+'-'+strtrim(string((x(xl)+dxp)*1.0e6),2)+'!7l!8m' 
            posz='!8z='+strtrim(string(z(zl)*1.0e6),2)+'-'+strtrim(string((z(zl)+dzp)*1.0e6),2)+'!7l!8m' 
            posz='!8z='+strtrim(string(zmin),2)+'-'+strtrim(string(zmax),2)+'!7l!8m' 
            number='!8N='+strtrim(string(number),2)
            xyouts,300,y0+ydim+2500.0,ttt,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+2200.0,posx,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+1900.0,posz,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+1600.0,number,align=0.0,charsize=1.0,color=0,/DEVICE

      DEVICE, /CLOSE
      SET_PLOT,'X'


;      pnum=2000
;      ypl=fltarr(pnum)
;      pypl=fltarr(pnum)
;      ymin=15.0
;      ymax=25.0
;      dyy=(ymax-ymin)/pnum
;      for i=0L,pnum-1 do begin
;         ypl(i)=i*dyy+ymin
;         pypl(i)=0.07*sin(8.0*6.2831853072*ypl(i)/20.0)+0.03*(ypl(i)-20.0)
;      endfor
      

      fo=strmid(nmfc,0,strlen(nmfc)-13)+'_yyp-pz'$
      +'_pz'+strtrim(string(pzmin),2)+'-'+strtrim(string(pzmax),2)$
      +'_m'+strtrim(string(mass(ml)),2)$
      +'_q'+strtrim(string(charge(ql)),2)$
      +'_t'+strmid(nmfc,strlen(nmfc)-13,5)$
      +'_x'+strtrim(string(xl),2)$
      +'_z'+strtrim(string(zl),2)$
      +'.eps'

      print,fo

      SET_PLOT,'PS'
      DEVICE,/color,/encapsulated,file=fo,XSIZE=XSIZE,YSIZE=YSIZE
      erase, 255      


      plot, ypl, pypl, xs=1, ys=1, charsize=1.6, psym=3, $
      POSITION=[x0,y0,x0+xdim,y0+ydim],/DEVICE, $ 
      xtitle='!8y/!7l!8m',$
      xticks=2,$
;      xtickn=[strtrim(string(1.e6*x(x1)),2), $
;              strtrim(string(0.5e6*(x(x2)+x(x1))),2), $
;              strtrim(string(1.e6*x(x2)),2)], $
      xr=1.0e6*[y(y1),y(y2)], $
      ytitle='!8p!Dy!N/!8p!Dz',$
      yticks=2,$
;      ytickn=[strtrim(string(1.e6*y(y1)),2), $
;              strtrim(string(0.5e6*(y(y2)+y(y1))),2), $
;              strtrim(string(1.e6*y(y2)),2)], $
      yr=[pymin,pymax], $
      color=0,background=255 

            zmin=min(zpl)
            zmax=max(zpl)
            number=n0*cori*dx*dy*dz*ld*ld*ld*pnum
            ttt='!8t='+strtrim(string(t*1.0e15),2)+'fs'
            posx='!8x='+strtrim(string(x(xl)*1.0e6),2)+'-'+strtrim(string((x(xl)+dxp)*1.0e6),2)+'!7l!8m' 
            posz='!8z='+strtrim(string(z(zl)*1.0e6),2)+'-'+strtrim(string((z(zl)+dzp)*1.0e6),2)+'!7l!8m' 
            posz='!8z='+strtrim(string(zmin),2)+'-'+strtrim(string(zmax),2)+'!7l!8m' 
            number='!8N='+strtrim(string(number),2)
            xyouts,300,y0+ydim+2500.0,ttt,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+2200.0,posx,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+1900.0,posz,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+1600.0,number,align=0.0,charsize=1.0,color=0,/DEVICE

      DEVICE, /CLOSE
      SET_PLOT,'X'

   endif

endfor
endfor


ydim=8000.00
xdim=ydim*(i3x-i3n+1)/(i3x-i3n+1)

rimx=6.0
rimy=6.0

XSIZE=xdim/1000.0+rimx
YSIZE=ydim/1000.0+rimy

x0=500.0*rimx
y0=500.0*rimy


for zl=z1+0L,z2,stepz do begin
for yl=y1+0L,y2,stepy do begin
for xl=x1+0L,x2,stepx do begin


   pnum=0L
   for i=0L,nqcount-1 do begin
      if (z(zl) le zp(i)) and (zp(i) le z(zl)+dzp) then begin
      if (y(yl) le yp(i)) and (yp(i) le y(yl)+dyp) then begin
      if (x(xl) le xp(i)) and (xp(i) le x(xl)+dxp) then begin
      if (mnp(i) eq mass(ml)) then begin
      if (qnp(i) eq charge(ql)) then begin
      if (pzmin le pzp(i)) and (pzp(i) le pzmax) then begin
         pnum=pnum+1
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

      pnum=0L
      for i=0L,nqcount-1 do begin
         if (z(zl) le zp(i)) and (zp(i) le z(zl)+dzp) then begin
         if (y(yl) le yp(i)) and (yp(i) le y(yl)+dyp) then begin
         if (x(xl) le xp(i)) and (xp(i) le x(xl)+dxp) then begin
         if (mnp(i) eq mass(ml)) then begin
         if (qnp(i) eq charge(ql)) then begin
         if (pzmin le pzp(i)) and (pzp(i) le pzmax) then begin
            xpl(pnum)=1.0e6*xp(i)
            ypl(pnum)=1.0e6*yp(i)
            zpl(pnum)=1.0e6*zp(i)
            pxpl(pnum)=pxp(i)/pzp(i)
            pypl(pnum)=pyp(i)/pzp(i)
            pxpl(pnum)=pxp(i)
            pypl(pnum)=pyp(i)
            pzpl(pnum)=pzp(i)
            qnpl(pnum)=qnp(i)
            mnpl(pnum)=mnp(i)
            lnpl(pnum)=lnp(i)
            pnum=pnum+1
         endif
         endif
         endif
         endif
         endif
         endif
      endfor


      nx=120
      ny=120
      dpx=(pxmax-pxmin)/nx
      dpy=(pymax-pymin)/ny
      xs=pxmin
      xu=pxmax
      ys=pymin
      yu=pymax
      nf=fltarr(nx+2,ny+2)

      p=0L
      for l=0L,pnum-1 do begin
 
         xq=pxpl(l)
         yq=pypl(l)
 
         if (xs lt xq) and (xq lt xu) then begin
         if (ys lt yq) and (yq lt yu) then begin
 
         k1=floor((xq-xs)/dpx)
         k2=floor((yq-ys)/dpy)
 
         xq=(xq-xs)/dpx-k1
         yq=(yq-ys)/dpy-k2
 
         h1=(1.0-xq)*(1.0-yq)
         h2=xq*(1.0-yq)
         h3=(1.0-xq)*yq
         h4=xq*yq
 
         nf(k1,k2)=nf(k1,k2)+h1
         nf(k1+1,k2)=nf(k1+1,k2)+h2
         nf(k1,k2+1)=nf(k1,k2+1)+h3
         nf(k1+1,k2+1)=nf(k1+1,k2+1)+h4

         endif
         endif

      endfor


      fo=strmid(nmfc,0,strlen(nmfc)-13)+'_pxpy-pz-scale'$
      +'_pz'+strtrim(string(pzmin),2)+'-'+strtrim(string(pzmax),2)$
      +'_m'+strtrim(string(mass(ml)),2)$
      +'_q'+strtrim(string(charge(ql)),2)$
      +'_t'+strmid(nmfc,strlen(nmfc)-13,5)$
      +'_x'+strtrim(string(xl),2)$
      +'_y'+strtrim(string(yl),2)$
      +'_z'+strtrim(string(zl),2)$
      +'.eps'

      print,fo


      SET_PLOT,'PS'
      DEVICE,/color,/encapsulated,file=fo,XSIZE=XSIZE,YSIZE=YSIZE
      erase, 255      


      sfxy=smooth(nf,4)
      rmin=min(sfxy)
      rmax=max(sfxy)
 
      CONTOUR,sfxy(0:nx,0:ny), $
      xs=1,ys=1,charsize=1.6, $
      POSITION=[x0,y0,x0+xdim,y0+ydim], $
      xtitle='!8p!Dx!N/!8p!Dz',$
      xticks=2,$
      xr=[xs,xu], $
      ytitle='!8p!Dy!N/!8p!Dz',$
      yticks=2,$
      yr=[ys,yu], $
      background=255,color=0, $
      /DEVICE,/NODATA

      ixy = BYTSCL(sfxy(0:nx,0:ny), MIN=rmin, MAX=rmax)
      TV, ixy, !X.WINDOW(0),!Y.WINDOW(0),$
      XSIZE=!X.WINDOW(1) - !X.WINDOW(0),$
      YSIZE=!Y.WINDOW(1) - !Y.WINDOW(0), /NORM
 
      if (rmax-rmin gt 0.0) then begin
         x11=(x0+xdim+500.0)/(xdim+2.0*x0)
         y11=y0/(ydim+2.0*y0)
         x22=(x0+xdim+800.0)/(xdim+2.0*x0)
         y22=(y0+ydim)/(ydim+2.0*y0)
         colorbar, division=2, NCOLORS=!D.Table_size, /right, /vertical, $
         format='(F8.2)',color=0, charsize=1.5, min=rmin, max=rmax, $
         /pscolor, POSITION=[x11,y11,x22,y22]
      endif

            zmin=min(zpl)
            zmax=max(zpl)
            number=n0*cori*dx*dy*dz*ld*ld*ld*pnum
            ttt='!8t='+strtrim(string(t*1.0e15),2)+'fs'
            posx='!8x='+strtrim(string(x(xl)*1.0e6),2)+'-'+strtrim(string((x(xl)+dxp)*1.0e6),2)+'!7l!8m' 
            posy='!8y='+strtrim(string(y(yl)*1.0e6),2)+'-'+strtrim(string((y(yl)+dyp)*1.0e6),2)+'!7l!8m' 
            posz='!8z='+strtrim(string(z(zl)*1.0e6),2)+'-'+strtrim(string((z(zl)+dzp)*1.0e6),2)+'!7l!8m' 
            posz='!8z='+strtrim(string(zmin),2)+'-'+strtrim(string(zmax),2)+'!7l!8m' 
            number='!8N='+strtrim(string(number),2)
            xyouts,300,y0+ydim+2500.0,ttt,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+2200.0,posx,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+1900.0,posy,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+1600.0,posz,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+1300.0,number,align=0.0,charsize=1.0,color=0,/DEVICE

      DEVICE, /CLOSE
      SET_PLOT,'X'


      fo=strmid(nmfc,0,strlen(nmfc)-13)+'_pxpy-pz'$
      +'_pz'+strtrim(string(pzmin),2)+'-'+strtrim(string(pzmax),2)$
      +'_m'+strtrim(string(mass(ml)),2)$
      +'_q'+strtrim(string(charge(ql)),2)$
      +'_t'+strmid(nmfc,strlen(nmfc)-13,5)$
      +'_x'+strtrim(string(xl),2)$
      +'_y'+strtrim(string(yl),2)$
      +'_z'+strtrim(string(zl),2)$
      +'.eps'

      print,fo


      SET_PLOT,'PS'
      DEVICE,/color,/encapsulated,file=fo,XSIZE=XSIZE,YSIZE=YSIZE
      erase, 255      


      plot, pxpl, pypl, xs=1, ys=1, charsize=1.6, psym=3, $
      POSITION=[x0,y0,x0+xdim,y0+ydim],/DEVICE, $
      xtitle='!8p!Dx!N/!8p!Dz',$
      xticks=2,$
;      xtickn=[strtrim(string(pxmin),2), $
;              strtrim(string(pxmid),2), $
;              strtrim(string(pxmax),2)], $
      xr=[pxmin,pxmax], $
      ytitle='!8p!Dy!N/!8p!Dz',$
      yticks=2,$
;      ytickn=[strtrim(string(pymin),2), $
;              strtrim(string(pymid),2), $
;              strtrim(string(pymax),2)], $
      yr=[pymin,pymax], $
      color=0,background=255
 
            zmin=min(zpl)
            zmax=max(zpl)
            number=n0*cori*dx*dy*dz*ld*ld*ld*pnum
            ttt='!8t='+strtrim(string(t*1.0e15),2)+'fs'
            posx='!8x='+strtrim(string(x(xl)*1.0e6),2)+'-'+strtrim(string((x(xl)+dxp)*1.0e6),2)+'!7l!8m' 
            posy='!8y='+strtrim(string(y(yl)*1.0e6),2)+'-'+strtrim(string((y(yl)+dyp)*1.0e6),2)+'!7l!8m' 
            posz='!8z='+strtrim(string(z(zl)*1.0e6),2)+'-'+strtrim(string((z(zl)+dzp)*1.0e6),2)+'!7l!8m' 
            posz='!8z='+strtrim(string(zmin),2)+'-'+strtrim(string(zmax),2)+'!7l!8m' 
            number='!8N='+strtrim(string(number),2)
            xyouts,300,y0+ydim+2500.0,ttt,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+2200.0,posx,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+1900.0,posy,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+1600.0,posz,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+1300.0,number,align=0.0,charsize=1.0,color=0,/DEVICE

      DEVICE, /CLOSE
      SET_PLOT,'X'

   endif

endfor
endfor
endfor


endfor
endfor
endfor

endif
endif
end
