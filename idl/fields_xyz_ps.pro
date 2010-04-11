rimx=10.0
rimy=10.0

XSIZE=xdim/1000.0+rimx
YSIZE=ydim/1000.0+rimy

x0=500.0*rimx
y0=500.0*rimy


fo=strmid(nmfc,0,strlen(nmfc)-13)$
+'_t'+strmid(nmfc,strlen(nmfc)-13,5)$
+'_3D.eps'

print,fo

SET_PLOT,'PS'
DEVICE,/color,/encapsulated,file=fo,XSIZE=XSIZE,YSIZE=YSIZE
erase, 255      

f=fltarr(z2-z1+1,x2-x1+1,y2-y1+1)
for k=0L,z2-z1 do begin
   for j=0L,y2-y1 do begin
      for i=0L,x2-x1 do begin
         f(k,i,j)=serie(x1+i,y1+j,z1+k)
      endfor
   endfor
endfor

nx=x2-x1+1
ny=y2-y1+1
nz=z2-z1+1

SHADE_VOLUME, f, rmax, VV, PP
CREATE_VIEW,ax=-70,az=55,xmax=nz,ymax=nx,zmax=ny,zoom=[0.6,0.6,0.6]

PLOT_3DBOX,[0,0],[0,0],[0,0],/t3d,charsize=3.0,$
xr=[0,nz-1],yr=[0,nx-1],zr=[0,ny-1],gridstyle=1,$
xticks=2,xtickname=[' ',' ',' '],xthick=1.0,$
yticks=2,ytickname=[' ',' ',' '],ythick=1.0,$
zticks=2,ztickname=[' ',' ',' '],zthick=1.0,$
background=255,color=0,/nodata

mypoly,VV,PP,ambient=0.4
TVLCT, rr,gg,bb,/GET

axis,0,0,0,xaxis=0,/noerase,/t3d,xticks=2,charsize=6.0,$ 
xrange=1.0e6*[z(z1),z(z2)],xthick=1.0,$
xtitle='!8z/!7l!8m',color=0
axis,0,0,0,yaxis=0,/noerase,/t3d,yticks=2,charsize=6.0,$ 
yrange=1.0e6*[x(x2),x(x1)],ythick=1.0,$
ytitle='!8y/!7l!8m',color=0
axis,nz-1,0,0,zaxis=0,/noerase,/t3d,zticks=2,charsize=3.5,$
zrange=1.0e6*[y(y1),y(y2)],zthick=1.0,$
ztitle='!8x/!7l!8m',color=0

ttt='!8t='+strtrim(string(t*1.0e15),2)+'fs' 
rmini='!8min='+strtrim(string(mini),2)
rmaxi='!8max='+strtrim(string(maxi),2)
xyouts,x0,y0-16000,ttt,align=1.0,charsize=5.0,color=0,/DEVICE
xyouts,x0,y0-18000,rmaxi,align=1.0,charsize=5.0,color=0,/DEVICE
xyouts,x0,y0-20000,rmini,align=1.0,charsize=5.0,color=0,/DEVICE

DEVICE, /CLOSE
SET_PLOT,'X'

;end
