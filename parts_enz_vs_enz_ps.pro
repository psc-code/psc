rimx=8.0
rimy=8.0

XSIZE=xdim/1000.0+rimx
YSIZE=ydim/1000.0+rimy

x0=700.0*rimx
y0=400.0*rimy

dxp=stepx*dx*ld
dyp=stepy*dy*ld
dzp=stepz*dz*ld


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
            if (z(zl)-dz*ld le zp(i)) and (zp(i) le z(zl)+dzp) then begin
            if (y(yl)-dy*ld le yp(i)) and (yp(i) le y(yl)+dyp) then begin
            if (x(xl)-dx*ld le xp(i)) and (xp(i) le x(xl)+dxp) then begin
            if (pxmin le pxp(i)) and (pxp(i) le pxmax) then begin
            if (pymin le pyp(i)) and (pyp(i) le pymax) then begin
            if (pzmin le pzp(i)) and (pzp(i) le pzmax) then begin
            if (mnp(i) eq mass(ml)) then begin
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
               if (z(zl)-dz*ld le zp(i)) and (zp(i) le z(zl)+dzp) then begin
               if (y(yl)-dy*ld le yp(i)) and (yp(i) le y(yl)+dyp) then begin
               if (x(xl)-dx*ld le xp(i)) and (xp(i) le x(xl)+dxp) then begin
               if (pxmin le pxp(i)) and (pxp(i) le pxmax) then begin
               if (pymin le pyp(i)) and (pyp(i) le pymax) then begin
               if (pzmin le pzp(i)) and (pzp(i) le pzmax) then begin
               if (mnp(i) eq mass(ml)) then begin
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
               endif
            endfor


            pbin=1000
            nnpl=fltarr(pbin)
            pl=fltarr(pbin)
            enpl=fltarr(pbin)


            if (enmin gt 0.0) then begin
               ppmin=sqrt(enmin/(0.2555*mass(ml)))
            endif
            if (enmin lt 0.0) then begin
               ppmin=-sqrt(-enmin/(0.2555*mass(ml)))
            endif
            if (enmax gt 0.0) then begin
               ppmax=sqrt(enmax/(0.2555*mass(ml)))
            endif
            if (enmax lt 0.0) then begin
               ppmax=-sqrt(-enmax/(0.2555*mass(ml)))
            endif


            for j=0L,pbin-1 do begin
               pl(j)=(ppmax-ppmin)*(j+0.5)/pbin+ppmin
            endfor
            for j=0L,pbin-1 do begin
               pe=pl(j)
               if (pe gt 0.0) then begin
                  enpl(j)=+0.2555*mass(ml)*pe*pe
               endif
               if (pe lt 0.0) then begin
                  enpl(j)=-0.2555*mass(ml)*pe*pe
               endif
            endfor


            number=n0*cori*nistep*step*dx*dy*dz*ld*ld*ld
            for j=0L,pbin-1 do begin
               nnpl(j)=0.0
            endfor
            for j=0L,pnum-1 do begin
               pe=pxpl(j)
               ll=floor(pbin*(pe-ppmin)/(ppmax-ppmin))
               if (0.0 lt ll) and (ll lt pbin-1) then begin
                  nnpl(ll)=nnpl(ll)+0.2555*mnpl(j)*wnpl(j)*pe*pe*number
               endif else begin
                  nnpl(pbin-1)=nnpl(pbin-1)+0.2555*mnpl(j)*wnpl(j)*pe*pe*number
               endelse
            endfor
            nnpl=nnpl*pbin/(enmax-enmin)
            snnpl=smooth(nnpl,10)


            nnmin=min(nnpl)
            if (nnmin eq 0.0) then begin
               nnmin=10000.0
            endif
            nnmax=max(nnpl)
            if (nnmax eq 0.0) then begin
               nnmax=20000.0
            endif


            fo=strmid(nmfc,0,strlen(nmfc)-13)+'_enx_vs_enx'$
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


            plot, enpl(20:pbin-20), snnpl(20:pbin-20), $
            xs=1, ys=1, charsize=1.6, psym=3, /ylog, $
            POSITION=[x0,y0,x0+xdim,y0+ydim],/DEVICE, $ 
            xtitle='!8E!Dx!N/MeV', $
            xticks=2,$
;            xtickn=[strtrim(string(enmin),2), $
;                    strtrim(string(enmid),2), $
;                    strtrim(string(enmax),2)], $
            xr=[enmin,enmax], $  
            ytitle='!8dE!Dx!N/dE!Dx', $
            yticks=2,$
            yr=[nnmin,nnmax], $
            color=0,background=255 

            number=0.0
            for i=0L,pnum-1 do begin      
               number=number+wnpl(i)
            endfor
            number=n0*cori*nistep*step*dx*dy*dz*ld*ld*ld*number

            energy=0.0
            for ll=0L,pbin-1 do begin
               energy=energy+1.0e6*nnpl(ll)*(enmax-enmin)/pbin
            endfor

            ttt='!8t='+strtrim(string(t*1.0e15),2)+'fs'
            posx='!8x='+strtrim(string((x(xl)-dx*ld)*1.0e6),2)+'-'+strtrim(string((x(xl)+dxp)*1.0e6),2)+'!7l!8m' 
            posy='!8y='+strtrim(string((y(yl)-dy*ld)*1.0e6),2)+'-'+strtrim(string((y(yl)+dyp)*1.0e6),2)+'!7l!8m' 
            posz='!8z='+strtrim(string((z(zl)-dz-ld)*1.0e6),2)+'-'+strtrim(string((z(zl)+dzp)*1.0e6),2)+'!7l!8m' 
            momx='!8px='+strtrim(string(pxmin),2)+'-'+strtrim(string(pxmax),2)+'!8mc'
            momy='!8py='+strtrim(string(pymin),2)+'-'+strtrim(string(pymax),2)+'!8mc'
            momz='!8pz='+strtrim(string(pzmin),2)+'-'+strtrim(string(pzmax),2)+'!8mc'
            mas='!8m='+strtrim(string(mass(ml)),2)+'!8m!De'
            number='!8N='+strtrim(string(number),2)
            energy='!8E='+strtrim(string(energy),2)+'!8eV'
            xyouts,300,y0+ydim+3300.0,ttt,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+3000.0,posx,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+2700.0,posy,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+2400.0,posz,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+2100.0,momx,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+1800.0,momy,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+1500.0,momz,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+1200.0,mas,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+900.0,number,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+600.0,energy,align=0.0,charsize=1.0,color=0,/DEVICE

            DEVICE, /CLOSE


            number=n0*cori*nistep*step*dx*dy*dz*ld*ld*ld
            for j=0L,pbin-1 do begin
               nnpl(j)=0.0
            endfor
            for j=0L,pnum-1 do begin
               pe=pypl(j)
               ll=floor(pbin*(pe-ppmin)/(ppmax-ppmin))
               if (0.0 lt ll) and (ll lt pbin-1) then begin
                  nnpl(ll)=nnpl(ll)+0.2555*mnpl(j)*wnpl(j)*pe*pe*number
               endif else begin
                  nnpl(pbin-1)=nnpl(pbin-1)+0.2555*mnpl(j)*wnpl(j)*pe*pe*number
               endelse
            endfor
            nnpl=nnpl*pbin/(enmax-enmin)
            snnpl=smooth(nnpl,10)


            nnmin=min(nnpl)
            if (nnmin eq 0.0) then begin
               nnmin=10000.0
            endif
            nnmax=max(nnpl)
            if (nnmax eq 0.0) then begin
               nnmax=20000.0
            endif


            fo=strmid(nmfc,0,strlen(nmfc)-13)+'_eny_vs_eny'$
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


            plot, enpl(20:pbin-20), snnpl(20:pbin-20), $
            xs=1, ys=1, charsize=1.6, psym=3, /ylog, $
            POSITION=[x0,y0,x0+xdim,y0+ydim],/DEVICE, $ 
            xtitle='!8E!Dy!N/MeV', $
            xticks=2,$
;            xtickn=[strtrim(string(enmin),2), $
;                    strtrim(string(enmid),2), $
;                    strtrim(string(enmax),2)], $
            xr=[enmin,enmax], $  
            ytitle='!8dE!Dy!N/dE!Dy', $
            yticks=2,$
            yr=[nnmin,nnmax], $
            color=0,background=255 

            number=0.0
            for i=0L,pnum-1 do begin      
               number=number+wnpl(i)
            endfor
            number=n0*cori*nistep*step*dx*dy*dz*ld*ld*ld*number

            energy=0.0
            for ll=0L,pbin-1 do begin
               energy=energy+1.0e6*nnpl(ll)*(enmax-enmin)/pbin
            endfor

            ttt='!8t='+strtrim(string(t*1.0e15),2)+'fs'
            posx='!8x='+strtrim(string((x(xl)-dx*ld)*1.0e6),2)+'-'+strtrim(string((x(xl)+dxp)*1.0e6),2)+'!7l!8m' 
            posy='!8y='+strtrim(string((y(yl)-dy*ld)*1.0e6),2)+'-'+strtrim(string((y(yl)+dyp)*1.0e6),2)+'!7l!8m' 
            posz='!8z='+strtrim(string((z(zl)-dz-ld)*1.0e6),2)+'-'+strtrim(string((z(zl)+dzp)*1.0e6),2)+'!7l!8m' 
            momx='!8px='+strtrim(string(pxmin),2)+'-'+strtrim(string(pxmax),2)+'!8mc'
            momy='!8py='+strtrim(string(pymin),2)+'-'+strtrim(string(pymax),2)+'!8mc'
            momz='!8pz='+strtrim(string(pzmin),2)+'-'+strtrim(string(pzmax),2)+'!8mc'
            mas='!8m='+strtrim(string(mass(ml)),2)+'!8m!De'
            number='!8N='+strtrim(string(number),2)
            energy='!8E='+strtrim(string(energy),2)+'!8eV'
            xyouts,300,y0+ydim+3300.0,ttt,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+3000.0,posx,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+2700.0,posy,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+2400.0,posz,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+2100.0,momx,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+1800.0,momy,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+1500.0,momz,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+1200.0,mas,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+900.0,number,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+600.0,energy,align=0.0,charsize=1.0,color=0,/DEVICE

            DEVICE, /CLOSE


            number=n0*cori*nistep*step*dx*dy*dz*ld*ld*ld
            for j=0L,pbin-1 do begin
               nnpl(j)=0.0
            endfor
            for j=0L,pnum-1 do begin
               pe=pzpl(j)
               ll=floor(pbin*(pe-ppmin)/(ppmax-ppmin))
               if (0.0 lt ll) and (ll lt pbin-1) then begin
                  nnpl(ll)=nnpl(ll)+0.2555*mnpl(j)*wnpl(j)*pe*pe*number
               endif else begin
                  nnpl(pbin-1)=nnpl(pbin-1)+0.2555*mnpl(j)*wnpl(j)*pe*pe*number
               endelse
            endfor
            nnpl=nnpl*pbin/(enmax-enmin)
            snnpl=smooth(nnpl,10)


            nnmin=min(nnpl)
            if (nnmin eq 0.0) then begin
               nnmin=10000.0
            endif
            nnmax=max(nnpl)
            if (nnmax eq 0.0) then begin
               nnmax=20000.0
            endif


            fo=strmid(nmfc,0,strlen(nmfc)-13)+'_enz_vs_enz'$
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


            plot, enpl(20:pbin-20), snnpl(20:pbin-20), $
            xs=1, ys=1, charsize=1.6, psym=3, /ylog, $
            POSITION=[x0,y0,x0+xdim,y0+ydim],/DEVICE, $ 
            xtitle='!8E!Dz!N/MeV', $
            xticks=2,$
;            xtickn=[strtrim(string(enmin),2), $
;                    strtrim(string(enmid),2), $
;                    strtrim(string(enmax),2)], $
            xr=[enmin,enmax], $  
            ytitle='!8dE!Dz!N/dE!Dz', $
            yticks=2,$
            yr=[nnmin,nnmax], $
            color=0,background=255 

            number=0.0
            for i=0L,pnum-1 do begin      
               number=number+wnpl(i)
            endfor
            number=n0*cori*nistep*step*dx*dy*dz*ld*ld*ld*number

            energy=0.0
            for ll=0L,pbin-1 do begin
               energy=energy+1.0e6*nnpl(ll)*(enmax-enmin)/pbin
            endfor

            ttt='!8t='+strtrim(string(t*1.0e15),2)+'fs'
            posx='!8x='+strtrim(string((x(xl)-dx*ld)*1.0e6),2)+'-'+strtrim(string((x(xl)+dxp)*1.0e6),2)+'!7l!8m' 
            posy='!8y='+strtrim(string((y(yl)-dy*ld)*1.0e6),2)+'-'+strtrim(string((y(yl)+dyp)*1.0e6),2)+'!7l!8m' 
            posz='!8z='+strtrim(string((z(zl)-dz-ld)*1.0e6),2)+'-'+strtrim(string((z(zl)+dzp)*1.0e6),2)+'!7l!8m' 
            momx='!8px='+strtrim(string(pxmin),2)+'-'+strtrim(string(pxmax),2)+'!8mc'
            momy='!8py='+strtrim(string(pymin),2)+'-'+strtrim(string(pymax),2)+'!8mc'
            momz='!8pz='+strtrim(string(pzmin),2)+'-'+strtrim(string(pzmax),2)+'!8mc'
            mas='!8m='+strtrim(string(mass(ml)),2)+'!8m!De'
            number='!8N='+strtrim(string(number),2)
            energy='!8E='+strtrim(string(energy),2)+'!8eV'
            xyouts,300,y0+ydim+3300.0,ttt,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+3000.0,posx,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+2700.0,posy,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+2400.0,posz,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+2100.0,momx,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+1800.0,momy,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+1500.0,momz,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+1200.0,mas,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+900.0,number,align=0.0,charsize=1.0,color=0,/DEVICE
            xyouts,300,y0+ydim+600.0,energy,align=0.0,charsize=1.0,color=0,/DEVICE

            DEVICE, /CLOSE
            SET_PLOT,'X'

         endif
      endfor
      endfor
   endfor
endfor
