; 2D FIELDS

fi=FINDFILE('net*.data.gz')

nfi=size(fi,/N_ELEMENTS)
loadct,5

	
for ifi=0,nfi-1 do begin


    print,fi(ifi)
    nmfc=fi(ifi)
    @FIELDS_EVOL

    ndim=size(serie,/DIMENSIONS)
    mini=min(serie)
    maxi=max(serie)


; data range to be plotted


   rmin=mini
   rmax=maxi
   rmin=0.5*mini
   rmax=0.5*maxi
;   rmin=-10.0
;   rmax=0.0


   x1=0
   x2=ndim(0)-2
   y1=0
   y2=ndim(1)-2
   z1=0
   z2=ndim(2)-2


   stepx=x2-x1+2
   stepx=5
   stepy=y2-y1+2
   stepy=5
   stepz=z2-z1+2
   stepz=5


   ydim=10000.0
   xdim=ydim*ndim(0)/ndim(1)
   @fields_xy_ps


   ydim=10000.0
   xdim=ydim*ndim(0)/ndim(2)
   @fields_xz_ps


   ydim=10000.0
   xdim=ydim*ndim(1)/ndim(2)
   @fields_yz_ps

endfor

set_plot,'X'
end
