rimx=6.0
rimy=6.0

XSIZE=xdim/1000.0+rimx
YSIZE=ydim/1000.0+rimy

x0=500.0*rimx
y0=500.0*rimy

dxp=stepx*dx*ld
dyp=stepy*dy*ld
dzp=stepz*dz*ld


for zl=z1+0L,z2,stepz do begin
   for yl=y1+0L,y2,stepy do begin
      for xl=x1+0L,x2,stepx do begin

         pnum=0L
         for i=0L,nqcount-1 do begin
            if (z(zl)-dz*ld le zp(i)) and (zp(i) le z(zl)+dzp) then begin
            if (y(yl)-dy*ld le yp(i)) and (yp(i) le y(yl)+dyp) then begin
            if (x(xl)-dx*ld le xp(i)) and (xp(i) le x(xl)+dxp) then begin
            if (pxmin le pxp(i)) and (pxp(i) le pxmax) then begin
            if (pymin le pyp(i)) and (pyp(i) le pymax) then begin
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
            wnpl=fltarr(pnum)

            pnum=0L
            for i=0L,nqcount-1 do begin
               if (z(zl)-dz*ld le zp(i)) and (zp(i) le z(zl)+dzp) then begin
               if (y(yl)-dy*ld le yp(i)) and (yp(i) le y(yl)+dyp) then begin
               if (x(xl)-dx*ld le xp(i)) and (xp(i) le x(xl)+dxp) then begin
               if (pxmin le pxp(i)) and (pxp(i) le pxmax) then begin
               if (pymin le pyp(i)) and (pyp(i) le pymax) then begin
               if (pzmin le pzp(i)) and (pzp(i) le pzmax) then begin
                  xpl(pnum)=xp(i)
                  ypl(pnum)=yp(i)
                  zpl(pnum)=zp(i)
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
            endfor


            fo=strmid(nmfc,0,strlen(nmfc)-13)+'_pxpy'$
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
            xtitle='!8p!Dx!N/!8mc',$
            xticks=2,$
;            xtickn=[strtrim(string(pxmin),2), $
;                    strtrim(string(pxmid),2), $
;                    strtrim(string(pxmax),2)], $
            xr=[pxmin,pxmax], $  
            ytitle='!8p!Dy!N/!8mc',$
            yticks=2,$
;            ytickn=[strtrim(string(pymin),2), $
;                    strtrim(string(pymid),2), $
;                    strtrim(string(pymax),2)], $
            yr=[pymin,pymax], $
            color=0,background=255 

            number=0.0
            for i=0L,pnum-1 do begin      
               number=number+wnpl(i)
            endfor
            number=n0*cori*nistep*step*dx*dy*dz*ld*ld*ld*number

            ttt='!8t='+strtrim(string(t*1.0e15),2)+'fs'
            posx='!8x='+strtrim(string((x(xl)-dx*ld)*1.0e6),2)+'-'+strtrim(string((x(xl)+dxp)*1.0e6),2)+'!7l!8m'
            posy='!8y='+strtrim(string((y(yl)-dy*ld)*1.0e6),2)+'-'+strtrim(string((y(yl)+dyp)*1.0e6),2)+'!7l!8m'
            posz='!8z='+strtrim(string((z(zl)-dz*ld)*1.0e6),2)+'-'+strtrim(string((z(zl)+dzp)*1.0e6),2)+'!7l!8m'
            momx='!8px='+strtrim(string(pxmin),2)+'-'+strtrim(string(pxmax),2)+'!8mc'
            momy='!8py='+strtrim(string(pymin),2)+'-'+strtrim(string(pymax),2)+'!8mc'
            momz='!8pz='+strtrim(string(pzmin),2)+'-'+strtrim(string(pzmax),2)+'!8mc'
            number='!8N='+strtrim(string(number),2)
            xyouts,300,y0+ydim+2500.0,ttt,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+2200.0,posx,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+1900.0,posy,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+1600.0,posz,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+1300.0,momz,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+1000.0,number,align=0.0,charsize=1.0,color=0,/DEVICE

            DEVICE, /CLOSE


            fo=strmid(nmfc,0,strlen(nmfc)-13)+'_pxpz'$
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
            xtitle='!8p!Dx!N/!8mc',$
            xticks=2,$
;            xtickn=[strtrim(string(pxmin),2), $
;                    strtrim(string(pxmid),2), $
;                    strtrim(string(pxmax),2)], $
            xr=[pxmin,pxmax], $
            ytitle='!8p!Dz!N/!8mc',$
            yticks=2,$
;            ytickn=[strtrim(string(pzmin),2), $
;                    strtrim(string(pzmid),2), $
;                    strtrim(string(pzmax),2)], $
            yr=[pzmin,pzmax], $
            color=0,background=255 

            number=0.0
            for i=0L,pnum-1 do begin      
               number=number+wnpl(i)
            endfor
            number=n0*cori*nistep*step*dx*dy*dz*ld*ld*ld*number

            ttt='!8t='+strtrim(string(t*1.0e15),2)+'fs'
            posx='!8x='+strtrim(string((x(xl)-dx*ld)*1.0e6),2)+'-'+strtrim(string((x(xl)+dxp)*1.0e6),2)+'!7l!8m'
            posy='!8y='+strtrim(string((y(yl)-dy*ld)*1.0e6),2)+'-'+strtrim(string((y(yl)+dyp)*1.0e6),2)+'!7l!8m'
            posz='!8z='+strtrim(string((z(zl)-dz*ld)*1.0e6),2)+'-'+strtrim(string((z(zl)+dzp)*1.0e6),2)+'!7l!8m'
            momx='!8px='+strtrim(string(pxmin),2)+'-'+strtrim(string(pxmax),2)+'!8mc'
            momy='!8py='+strtrim(string(pymin),2)+'-'+strtrim(string(pymax),2)+'!8mc'
            momz='!8pz='+strtrim(string(pzmin),2)+'-'+strtrim(string(pzmax),2)+'!8mc'
            number='!8N='+strtrim(string(number),2)
            xyouts,300,y0+ydim+2500.0,ttt,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+2200.0,posx,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+1900.0,posy,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+1600.0,posz,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+1300.0,momy,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+1000.0,number,align=0.0,charsize=1.0,color=0,/DEVICE

            DEVICE, /CLOSE


            fo=strmid(nmfc,0,strlen(nmfc)-13)+'_pypz'$
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
            xtitle='!8p!Dy!N/!8mc',$
            xticks=2,$
;            xtickn=[strtrim(string(pymin),2), $
;                    strtrim(string(pymid),2), $
;                    strtrim(string(pymax),2)], $
            xr=[pymin,pymax], $
            ytitle='!8p!Dz!N/!8mc',$
            yticks=2,$
;            ytickn=[strtrim(string(pzmin),2), $
;                    strtrim(string(pzmid),2), $
;                    strtrim(string(pzmax),2)], $
            yr=[pzmin,pzmax], $
            color=0,background=255 

            number=0.0
            for i=0L,pnum-1 do begin      
               number=number+wnpl(i)
            endfor
            number=n0*cori*nistep*step*dx*dy*dz*ld*ld*ld*number

            ttt='!8t='+strtrim(string(t*1.0e15),2)+'fs'
            posx='!8x='+strtrim(string((x(xl)-dx*ld)*1.0e6),2)+'-'+strtrim(string((x(xl)+dxp)*1.0e6),2)+'!7l!8m'
            posy='!8y='+strtrim(string((y(yl)-dy*ld)*1.0e6),2)+'-'+strtrim(string((y(yl)+dyp)*1.0e6),2)+'!7l!8m'
            posz='!8z='+strtrim(string((z(zl)-dz*ld)*1.0e6),2)+'-'+strtrim(string((z(zl)+dzp)*1.0e6),2)+'!7l!8m'
            momx='!8px='+strtrim(string(pxmin),2)+'-'+strtrim(string(pxmax),2)+'!8mc'
            momy='!8py='+strtrim(string(pymin),2)+'-'+strtrim(string(pymax),2)+'!8mc'
            momz='!8pz='+strtrim(string(pzmin),2)+'-'+strtrim(string(pzmax),2)+'!8mc'
            number='!8N='+strtrim(string(number),2)
            xyouts,300,y0+ydim+2500.0,ttt,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+2200.0,posx,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+1900.0,posy,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+1600.0,posz,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+1300.0,momx,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+1000.0,number,align=0.0,charsize=1.0,color=0,/DEVICE

            DEVICE, /CLOSE
            SET_PLOT,'X'

         endif
      endfor
   endfor
endfor
