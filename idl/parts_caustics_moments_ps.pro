rimx=6.0
rimy=6.0

XSIZE=xdim/1000.0+rimx
YSIZE=ydim/1000.0+rimy

x0=500.0*rimx
y0=500.0*rimy

dxp=0.5*stepx*dx*ld
dyp=0.5*stepy*dy*ld
dzp=0.5*stepz*dz*ld
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


pxmin=-0.5
pxmax=0.5
pymin=-0.5
pymax=0.5


for zl=z1+0L,z2,stepz do begin
   for yl=y1+0L,y2,stepy do begin
      for xl=x1+0L,x2,stepx do begin
      for ml=0,mnum do begin

         pnum=0L
         for i=0L,nqcount-1 do begin
            if (z(zl)-dzp le zp(i)) and (zp(i) le z(zl)+dzp) then begin
            if (y(yl)-dyp le yp(i)) and (yp(i) le y(yl)+dyp) then begin
            if (x(xl)-dxp le xp(i)) and (xp(i) le x(xl)+dxp) then begin
            if (mnp(i) eq mass(ml)) then begin
               pe=0.511*mnp(i)*(sqrt(1.0+pxp(i)*pxp(i) $
                                        +pyp(i)*pyp(i) $
                                        +pzp(i)*pzp(i))-1.0)
               if (pzmin le pzp(i)) and (pzp(i) lt pzmax) then begin
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
               if (z(zl)-dzp le zp(i)) and (zp(i) le z(zl)+dzp) then begin
               if (y(yl)-dyp le yp(i)) and (yp(i) le y(yl)+dyp) then begin
               if (x(xl)-dxp le xp(i)) and (xp(i) le x(xl)+dxp) then begin
               if (mnp(i) eq mass(ml)) then begin
                  pe=0.511*mnp(i)*(sqrt(1.0+pxp(i)*pxp(i) $
                                           +pyp(i)*pyp(i) $
                                           +pzp(i)*pzp(i))-1.0)
                  if (pzmin le pzp(i)) and (pzp(i) lt pzmax) then begin
                     xpl(pnum)=1.0e6*xp(i)
                     ypl(pnum)=1.0e6*yp(i)
                     zpl(pnum)=1.0e6*zp(i)
                     pxpl(pnum)=pxp(i)/pzp(i)
                     pypl(pnum)=pyp(i)/pzp(i)
                     pzpl(pnum)=pzp(i)/pzp(i)
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


            fo=strmid(nmfc,0,strlen(nmfc)-13)+'_pxpy:pz'$
            +'_m'+strtrim(string(mass(ml)),2)$
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
;            xtickn=[strtrim(string(pxmin),2), $
;                    strtrim(string(pxmid),2), $
;                    strtrim(string(pxmax),2)], $
            xr=[pxmin,pxmax], $  
            ytitle='!8p!Dy!N/!8p!Dz',$
            yticks=2,$
;            ytickn=[strtrim(string(pymin),2), $
;                    strtrim(string(pymid),2), $
;                    strtrim(string(pymax),2)], $
            yr=[pymin,pymax], $
            color=0,background=255 

            ttt='!8t='+strtrim(string(t*1.0e15),2)+'fs'
            pos1='!8x='+strtrim(string(x(xl)*1.0e6),2)+'!7l!8m' 
            pos2='!8y='+strtrim(string(y(yl)*1.0e6),2)+'!7l!8m' 
            pos3='!8z='+strtrim(string(z(zl)*1.0e6),2)+'!7l!8m' 
            xyouts,300,y0+ydim+2500.0,ttt,align=0.0,charsize=1.3,color=0,/DEVICE
;            xyouts,300,y0+ydim+2000.0,pos1,align=0.0,charsize=1.3,color=0,/DEVICE
;            xyouts,300,y0+ydim+1500.0,pos2,align=0.0,charsize=1.3,color=0,/DEVICE
;            xyouts,300,y0+ydim+1000.0,pos3,align=0.0,charsize=1.3,color=0,/DEVICE

            DEVICE, /CLOSE


            fo=strmid(nmfc,0,strlen(nmfc)-13)+'_pxpz:pz'$
            +'_m'+strtrim(string(mass(ml)),2)$
            +'_t'+strmid(nmfc,strlen(nmfc)-13,5)$
            +'_x'+strtrim(string(xl),2)$
            +'_y'+strtrim(string(yl),2)$
            +'_z'+strtrim(string(zl),2)$
            +'.eps'

            print,fo

            SET_PLOT,'PS'
            DEVICE,/color,/encapsulated,file=fo,XSIZE=XSIZE,YSIZE=YSIZE
            erase, 255      


            plot, pxpl, pzpl, xs=1, ys=1, charsize=1.6, psym=3, $
            POSITION=[x0,y0,x0+xdim,y0+ydim],/DEVICE, $ 
            xtitle='!8p!Dx!N/!8p!Dz',$
            xticks=2,$
;            xtickn=[strtrim(string(pxmin),2), $
;                    strtrim(string(pxmid),2), $
;                    strtrim(string(pxmax),2)], $
            xr=[pxmin,pxmax], $
            ytitle='!8p!Dz!N/!8p!Dz',$
            yticks=2,$
;            ytickn=[strtrim(string(pzmin),2), $
;                    strtrim(string(pzmid),2), $
;                    strtrim(string(pzmax),2)], $
            yr=[pzmin,pzmax], $
            color=0,background=255 

            ttt='!8t='+strtrim(string(t*1.0e15),2)+'fs'
            pos1='!8x='+strtrim(string(x(xl)*1.0e6),2)+'!7l!8m' 
            pos2='!8y='+strtrim(string(y(yl)*1.0e6),2)+'!7l!8m' 
            pos3='!8z='+strtrim(string(z(zl)*1.0e6),2)+'!7l!8m' 
;            xyouts,300,y0+ydim+2500.0,ttt,align=0.0,charsize=1.3,color=0,/DEVICE
;            xyouts,300,y0+ydim+2000.0,pos1,align=0.0,charsize=1.3,color=0,/DEVICE
;            xyouts,300,y0+ydim+1500.0,pos2,align=0.0,charsize=1.3,color=0,/DEVICE
;            xyouts,300,y0+ydim+1000.0,pos3,align=0.0,charsize=1.3,color=0,/DEVICE

            DEVICE, /CLOSE


            fo=strmid(nmfc,0,strlen(nmfc)-13)+'_pypz:pz'$
            +'_m'+strtrim(string(mass(ml)),2)$
            +'_t'+strmid(nmfc,strlen(nmfc)-13,5)$
            +'_x'+strtrim(string(xl),2)$
            +'_y'+strtrim(string(yl),2)$
            +'_z'+strtrim(string(zl),2)$
            +'.eps'

            print,fo

            SET_PLOT,'PS'
            DEVICE,/color,/encapsulated,file=fo,XSIZE=XSIZE,YSIZE=YSIZE
            erase, 255      


            plot, pypl, pzpl, xs=1, ys=1, charsize=1.6, psym=3, $
            POSITION=[x0,y0,x0+xdim,y0+ydim],/DEVICE, $ 
            xtitle='!8p!Dy!N/!8p!Dz',$
            xticks=2,$
;            xtickn=[strtrim(string(pymin),2), $
;                    strtrim(string(pymid),2), $
;                    strtrim(string(pymax),2)], $
            xr=[pymin,pymax], $
            ytitle='!8p!Dz!N/!8p!Dz',$
            yticks=2,$
;            ytickn=[strtrim(string(pzmin),2), $
;                    strtrim(string(pzmid),2), $
;                    strtrim(string(pzmax),2)], $
            yr=[pzmin,pzmax], $
            color=0,background=255 

            ttt='!8t='+strtrim(string(t*1.0e15),2)+'fs'
            pos1='!8x='+strtrim(string(x(xl)*1.0e6),2)+'!7l!8m' 
            pos2='!8y='+strtrim(string(y(yl)*1.0e6),2)+'!7l!8m' 
            pos3='!8z='+strtrim(string(z(zl)*1.0e6),2)+'!7l!8m' 
;            xyouts,300,y0+ydim+2500.0,ttt,align=0.0,charsize=1.3,color=0,/DEVICE
;            xyouts,300,y0+ydim+2000.0,pos1,align=0.0,charsize=1.3,color=0,/DEVICE
;            xyouts,300,y0+ydim+1500.0,pos2,align=0.0,charsize=1.3,color=0,/DEVICE
;            xyouts,300,y0+ydim+1000.0,pos3,align=0.0,charsize=1.3,color=0,/DEVICE

            DEVICE, /CLOSE
            SET_PLOT,'X'

         endif

      endfor
      endfor
   endfor
endfor

