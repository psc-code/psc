; 3D FIELDS

fi=FINDFILE('nit00000.data.gz')
nfi=size(fi,/N_ELEMENTS)
loadct,5


for ifi=0,nfi-1 do begin
   print,fi(ifi)
   nmfc=fi(ifi)
   @FIELDS_EVOL


   ndim=size(serie,/DIMENSIONS) 
   mini=min(serie)
   maxi=max(serie)
   rmin=0.00
   rmax=0.07


   x1=0
   x2=ndim(0)-1
   y1=0
   y2=ndim(1)-1
   z1=0
   z2=ndim(2)-1


   xdim=16000.0
   ydim=16000.0
   @fields_xyz_ps

endfor

set_plot,'X'
end
