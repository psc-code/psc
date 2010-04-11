rimx=6.0
rimy=6.0

XSIZE=xdim/1000.0+rimx
YSIZE=ydim/1000.0+rimy

x0=700.0*rimx
y0=500.0*rimy

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


for ml=0L,mnum do begin

   pnum=0L
   for i=0L,nqcount-1 do begin
      if (zmin le zp(i)) and (zp(i) lt zmax) then begin
      if (ymin le yp(i)) and (yp(i) lt ymax) then begin
      if (xmin le xp(i)) and (xp(i) lt xmax) then begin
      if (mnp(i) eq mass(ml)) then begin
         pnum=pnum+1
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
         if (zmin le zp(i)) and (zp(i) lt zmax) then begin
         if (ymin le yp(i)) and (yp(i) lt ymax) then begin
         if (xmin le xp(i)) and (xp(i) lt xmax) then begin
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
                  pnum=pnum+1
          endif
          endif
          endif
          endif
      endfor

      enbin=10000
      nnpl=fltarr(enbin)
      enpl=fltarr(enbin)

      for j=0L,enbin-1 do begin
         enpl(j)=(enmax-enmin)*(j+0.5)/enbin+enmin
      endfor

      for j=0L,pnum-1 do begin
         spe=pzpl(j)/abs(pzpl(j)+1.0e-10)
         if (spe ge 0.0) then begin
            pe=0.511*mnpl(j)*(sqrt(1.0+pxpl(j)*pxpl(j) $
                                      +pypl(j)*pypl(j) $
                                      +pzpl(j)*pzpl(j))-1.0)
         endif
         if (spe lt 0.0) then begin
            pe=-0.511*mnpl(j)*(sqrt(1.0+pxpl(j)*pxpl(j) $
                                       +pypl(j)*pypl(j) $
                                       +pzpl(j)*pzpl(j))-1.0)
         endif

         ll=floor(enbin*(pe-enmin)/(enmax-enmin))
         if (0.0 lt ll) and (ll lt enbin-1) then begin
            nnpl(ll)=nnpl(ll)+1.0*abs(pe)
;            nnpl(ll)=nnpl(ll)+1.0
         endif else begin
            nnpl(enbin-1)=nnpl(enbin-1)+1.0*abs(pe)
;            nnpl(enbin-1)=nnpl(enbin-1)+1.0
         endelse
      endfor
      snnpl=smooth(nnpl,20)


      fo=strmid(nmfc,0,strlen(nmfc)-13)+'_energy'$
      +'_m'+strtrim(string(mass(ml)),2)$
      +'_t'+strmid(nmfc,strlen(nmfc)-13,5)$
      +'.eps'

      print,fo
  
      SET_PLOT,'PS'
      DEVICE,/color,/encapsulated,file=fo,XSIZE=XSIZE,YSIZE=YSIZE
      erase, 255      


      TVLCT, [0,255,0,0], [0,0,255,0], [0,0,0,255]
      plot, enpl(20:enbin-20), snnpl(20:enbin-20), $
      xs=1, ys=1, charsize=1.6, psym=3, /ylog, $
      POSITION=[x0,y0,x0+xdim,y0+ydim],/DEVICE, $ 
      xtitle='!8E/MeV', $
      xticks=2,$
;      xtickn=[strtrim(string(enmin),2), $
;              strtrim(string(enmid),2), $
;              strtrim(string(enmax),2)], $
      xr=[enmin,enmax], $  
      ytitle='!8a.u', $
      yticks=2,$
      yr=[nnmin,nnmax], $
      color=0,background=255 

      ttt='!8t='+strtrim(string(t*1.0e15),2)+'fs'
      xyouts,300,y0+ydim+2500.0,ttt,align=0.0,charsize=1.3,color=0,/DEVICE

      DEVICE, /CLOSE
      SET_PLOT,'X'

   endif
endfor
