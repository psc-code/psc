xdim=6000
ydim=6000

rimx=6.0
rimy=6.0

XSIZE=xdim/1000.0+rimx
YSIZE=ydim/1000.0+rimy

x0=500.0*rimx
y0=500.0*rimy

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

   pnum=0L
   for i=0L,nqcount-1 do begin
      if (zmin le zp(i)) and (zp(i) lt zmax) then begin
      if (ymin le yp(i)) and (yp(i) lt ymax) then begin
      if (xmin le xp(i)) and (xp(i) lt xmax) then begin
         if (mnp(i) eq mass(ml)) then begin
         if (qnp(i) eq charge(ql)) then begin
         if (pzmin le pzp(i)) and (pzp(i) lt pzmax) then begin
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
         if (zmin le zp(i)) and (zp(i) lt zmax) then begin
         if (ymin le yp(i)) and (yp(i) lt ymax) then begin
         if (xmin le xp(i)) and (xp(i) lt xmax) then begin
            if (mnp(i) eq mass(ml)) then begin
            if (qnp(i) eq charge(ql)) then begin
            if (pzmin le pzp(i)) and (pzp(i) lt pzmax) then begin
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
         endif
      endfor

   endif


endfor
endfor

endif
endif
