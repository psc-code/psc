; 1D FIELDS

fi=FINDFILE('e*t*.data.gz')
nfi=size(fi,/N_ELEMENTS)
;loadct,5


for ifi=0,nfi-1 do begin

   print,fi(ifi)
   nmfc=fi(ifi)
   @FIELDS_EVOL


   ndim=size(serie,/DIMENSIONS) 
   mini=min(serie)
   maxi=max(serie)
   rmin=-0.5
   rmax=0.5
   if (rmin eq 0.0) then begin
      rmin=-0.01
   endif
   if (rmax eq 0.0) then begin
      rmax=+0.01
   endif


   x1=0
   x2=ndim(0)-2
   y1=0
   y2=ndim(1)-2
   z1=0
   z2=ndim(2)-2


   ydim=10000.0
   xdim=10000.0
   @fields_ps


endfor

set_plot,'X'
end
