rimx=6.0
rimy=6.0

XSIZE=xdim/1000.0+rimx
YSIZE=ydim/1000.0+rimy

x0=500.0*rimx
y0=500.0*rimy

dxp=0.5*stepx*dx*ld
dyp=0.5*stepy*dy*ld
dzp=0.5*stepz*dz*ld


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


for zl=z1+0L,z2,stepz do begin
   for yl=y1+0L,y2,stepy do begin
      for xl=x1+0L,x2,stepx do begin
      for ml=0L,mnum do begin

         pnum=0L
         for i=0L,nqcount-1 do begin
            if (z(zl)-dzp le zp(i)) and (zp(i) le z(zl)+dzp) then begin
            if (y(yl)-dyp le yp(i)) and (yp(i) le y(yl)+dyp) then begin
            if (x(xl)-dxp le xp(i)) and (xp(i) le x(xl)+dxp) then begin
            if (mnp(i) eq mass(ml)) then begin
               pe=0.511*mnp(i)*(sqrt(1.0+pxp(i)*pxp(i) $
                                        +pyp(i)*pyp(i) $
                                        +pzp(i)*pzp(i))-1.0)
               if (enmin le pe) and (pe lt enmax) then begin
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
                  if (enmin le pe) and (pe lt enmax) then begin
                     xpl(pnum)=1.0e6*xp(i)
                     ypl(pnum)=1.0e6*yp(i)
                     zpl(pnum)=1.0e6*zp(i)
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

            xmin=min(xpl)
            xmax=max(xpl)
            ymin=min(ypl)
            ymax=max(ypl)
            zmin=min(zpl)
            zmax=max(zpl)

            xbin=x2-x1+1
            ybin=y2-y1+1
            zbin=z2-z1+1
            axpl=fltarr(xbin)
            aypl=fltarr(ybin)
            azpl=fltarr(zbin)
            nxpl=fltarr(xbin)
            nypl=fltarr(ybin)
            nzpl=fltarr(zbin)

            for j=0L,xbin-1 do begin
               axpl(j)=(xmax-xmin)*(j+0.5)/xbin+xmin
            endfor
            for j=0L,ybin-1 do begin
               aypl(j)=(ymax-ymin)*(j+0.5)/ybin+ymin
            endfor
            for j=0L,zbin-1 do begin
               azpl(j)=(zmax-zmin)*(j+0.5)/zbin+zmin
            endfor

            for j=0L,pnum-1 do begin
               xx=xpl(j)
               ll=floor(xbin*(xx-xmin)/(xmax-xmin))
               if (0.0 lt ll) and (ll lt xbin-1) then begin
                  nxpl(ll)=nxpl(ll)+1.0
               endif
               yy=ypl(j)
               ll=floor(ybin*(yy-ymin)/(ymax-ymin))
               if (0.0 lt ll) and (ll lt ybin-1) then begin
                  nypl(ll)=nypl(ll)+1.0
               endif
               zz=zpl(j)
               ll=floor(zbin*(zz-zmin)/(zmax-zmin))
               if (0.0 lt ll) and (ll lt zbin-1) then begin
                  nzpl(ll)=nzpl(ll)+1.0
               endif
            endfor

            nxdim=size(nxpl,/DIMENSIONS) 
            if (nxdim(0) ge 2) then begin
               snxpl=smooth(nxpl,20)
            endif else begin
               snxpl=nxpl
            endelse
            nydim=size(nypl,/DIMENSIONS) 
            if (nydim(0) ge 2) then begin
               snypl=smooth(nypl,20)
            endif else begin
               snypl=nypl
            endelse
            nzdim=size(nzpl,/DIMENSIONS) 
            if (nzdim(0) ge 2) then begin
               snzpl=smooth(nzpl,20)
            endif else begin
               snzpl=nzpl
            endelse


            nxmin=0.001
            nxmax=1.05*max(snxpl)
            nymin=0.001
            nymax=1.05*max(snypl)
            nzmin=0.001
            nzmax=1.05*max(snzpl)

;            fo=strmid(nmfc,0,strlen(nmfc)-13)+'_denx'$
;            +'_m'+strtrim(string(mass(ml)),2)$
;            +'_t'+strmid(nmfc,strlen(nmfc)-13,5)$
;            +'_x'+strtrim(string(xl),2)$
;            +'_y'+strtrim(string(yl),2)$
;            +'_z'+strtrim(string(zl),2)$
;            +'.eps'

;            print,fo
  
;            SET_PLOT,'PS'
;            DEVICE,/color,/encapsulated,file=fo,XSIZE=XSIZE,YSIZE=YSIZE
;            erase, 255      


;            plot, axpl(20:xbin-20), snxpl(20:xbin-20), $
;            xs=1, ys=1, charsize=1.6, psym=3, /ylog, $
;            POSITION=[x0,y0,x0+xdim,y0+ydim],/DEVICE, $ 
;            xtitle='!8x/!7l!8m',$
;            xticks=2,$
;            xtickn=[strtrim(string(enmin),2), $
;                    strtrim(string(enmid),2), $
;                    strtrim(string(enmax),2)], $
;            xr=[xmin,xmax], $  
;            ytitle='!8a.u', $
;            yticks=2,$
;            yr=[nxmin,nxmax], $
;            color=0,background=255 

;            ttt='!8t='+strtrim(string(t*1.0e15),2)+'fs'
;            pos1='!8x='+strtrim(string(x(xl)*1.0e6),2)+'!7l!8m' 
;            pos2='!8y='+strtrim(string(y(yl)*1.0e6),2)+'!7l!8m' 
;            pos3='!8z='+strtrim(string(z(zl)*1.0e6),2)+'!7l!8m' 
;            xyouts,300,y0+ydim+2500.0,ttt,align=0.0,charsize=1.3,color=0,/DEVICE
;            xyouts,300,y0+ydim+2000.0,pos1,align=0.0,charsize=1.3,color=0,/DEVICE
;            xyouts,300,y0+ydim+1500.0,pos2,align=0.0,charsize=1.3,color=0,/DEVICE
;            xyouts,300,y0+ydim+1000.0,pos3,align=0.0,charsize=1.3,color=0,/DEVICE

;            DEVICE, /CLOSE

            fo=strmid(nmfc,0,strlen(nmfc)-13)+'_deny'$
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


            plot, aypl(20:ybin-20), snypl(20:ybin-20), $
            xs=1, ys=1, charsize=1.6, psym=3,$
            POSITION=[x0,y0,x0+xdim,y0+ydim],/DEVICE, $ 
            xtitle='!8y/!7l!8m',$
            xticks=2,$
            xtickn=['0.0','7.5','15.0'],$
;            xtickn=[strtrim(string(enmin),2), $
;                    strtrim(string(enmid),2), $
;                    strtrim(string(enmax),2)], $
;            xr=[ymin,ymax], $  
            ytitle='!8a.u', $
            yticks=2,$
;            ytickn=['0.0','6.0','12.0'],$
;            yr=[nymin,nymax], $
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

            fo=strmid(nmfc,0,strlen(nmfc)-13)+'_denz'$
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


            plot, azpl(20:zbin-20), snzpl(20:zbin-20), $
            xs=1, ys=1, charsize=1.6, psym=3, $
            POSITION=[x0,y0,x0+xdim,y0+ydim],/DEVICE, $ 
            xtitle='!8z/!7l!8m',$
            xticks=2,$
;            xtickn=[strtrim(string(enmin),2), $
;                    strtrim(string(enmid),2), $
;                    strtrim(string(enmax),2)], $
            xr=[zmin,zmax], $  
            ytitle='!8a.u', $
            yticks=2,$
            yr=[nzmin,nzmax], $
            color=0,background=255 

            ttt='!8t='+strtrim(string(t*1.0e15),2)+'fs'
            pos1='!8x='+strtrim(string(x(xl)*1.0e6),2)+'!7l!8m' 
            pos2='!8y='+strtrim(string(y(yl)*1.0e6),2)+'!7l!8m' 
            pos3='!8z='+strtrim(string(z(zl)*1.0e6),2)+'!7l!8m' 
            xyouts,300,y0+ydim+2500.0,ttt,align=0.0,charsize=1.3,color=0,/DEVICE
            xyouts,300,y0+ydim+2000.0,pos1,align=0.0,charsize=1.3,color=0,/DEVICE
            xyouts,300,y0+ydim+1500.0,pos2,align=0.0,charsize=1.3,color=0,/DEVICE
            xyouts,300,y0+ydim+1000.0,pos3,align=0.0,charsize=1.3,color=0,/DEVICE

            DEVICE, /CLOSE
            SET_PLOT,'X'

         endif
      endfor
      endfor
   endfor
endfor
