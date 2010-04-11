nmfc='IMAGE.data.gz'
openr,1,nmfc,/COMPRESS


readf,1,x1,x2,y1,y2,m
readf,1,pnumx,pnumy,pnumz,pnum,nmax
readf,1,c,beta,pi
readf,1,mi,ti
readf,1,xmin,xmax,xmid
readf,1,ymin,ymax,ymid
readf,1,zmin,zmax,zmid
readf,1,pxmin,pxmax,pymin,pymax
readf,1,z0,r0
readf,1,dx,dy,dt
readf,1,er,ez
readf,1,ax,ay,az


xp=fltarr(m+1)
yp=fltarr(m+1)
zp=fltarr(m+1)
pxp=fltarr(m+1)
pyp=fltarr(m+1)
pzp=fltarr(m+1)
out_pxp=fltarr(m+1)
out_pyp=fltarr(m+1)


   datin=fltarr(8)
   for n=0L,nmax do begin

      for j=0L,m do begin
         readf,1,datin
         xp(j)=datin(0)
         yp(j)=datin(1)
         zp(j)=datin(2)
         pxp(j)=datin(3)
         pyp(j)=datin(4)
         pzp(j)=datin(5)
         out_pxp(j)=datin(6)
         out_pyp(j)=datin(7)
      endfor


; DATA OUTPUT


      xdim=20000.0
      ydim=8000.0

      rimx=6.0
      rimy=6.0

      XSIZE=xdim/1000.0+rimx
      YSIZE=ydim/1000.0+rimy

      x0=500.0*rimx
      y0=500.0*rimy


      nx=x2-x1+1
      ny=200
      xs=x1*dx
      xu=x2*dx
      ys=pxmin
      yu=pxmax
      dpx=(yu-ys)/ny
      nf=fltarr(nx+1,ny+2)


      for l=0L,pnum-1 do begin
 
         xq=xp(l)
         yq=out_pxp(l)
 
         if (xs lt xq) and (xq lt xu) then begin
         if (ys lt yq) and (yq lt yu) then begin
 
         k1=floor(xq/dx)
         k2=floor(yq/dpx+0.5*(yu-ys)/dpx)
 
         xq=xp(l)/dx-k1
         yq=(out_pxp(l)+0.5*(yu-ys))/dpx-k2
 
         h1=(1.0-xq)*(1.0-yq)
         h2=xq*(1.0-yq)
         h3=(1.0-xq)*yq
         h4=xq*yq
 
         nf(k1,k2)=nf(k1,k2)+h1
         nf(k1+1,k2)=nf(k1+1,k2)+h2
         nf(k1,k2+1)=nf(k1,k2+1)+h3
         nf(k1+1,k2+1)=nf(k1+1,k2+1)+h4

         endif
         endif

      endfor


      fo='I_xxp-pz-image'$
      +'_t'+strtrim(string(tt),2)$
      +'.eps'
      print,fo


      SET_PLOT,'PS'
      DEVICE,/color,/encapsulated,file=fo,XSIZE=XSIZE,YSIZE=YSIZE
      erase, 255      


      sfxy=nf
      rmin=min(sfxy)
      rmax=0.5*max(sfxy)
 
      CONTOUR,sfxy(0:nx,0:ny), $
      xs=1,ys=1,charsize=1.6, $
      POSITION=[x0,y0,x0+xdim,y0+ydim], $
      xtitle='!8x/!8m', $
      xticks=2,$
      xr=[xs,xu], $
      ytitle='!8p!Dx!N/!8p!Dz',$
      yticks=2,$
      yr=[ys,yu], $
      background=255,color=0, $
      /DEVICE,/NODATA

      ixy = BYTSCL(sfxy(0:nx,0:ny), MIN=rmin, MAX=rmax)
      TV, ixy, !X.WINDOW(0),!Y.WINDOW(0),$
      XSIZE=!X.WINDOW(1) - !X.WINDOW(0),$
      YSIZE=!Y.WINDOW(1) - !Y.WINDOW(0), /NORM
;      oplot, linex, out_en_pxpl, psym=3, color=300.0, thick=2.0
 
      if (rmax-rmin gt 0.0) then begin
         x11=(x0+xdim+500.0)/(xdim+2.0*x0)
         y11=y0/(ydim+2.0*y0)
         x22=(x0+xdim+800.0)/(xdim+2.0*x0)
         y22=(y0+ydim)/(ydim+2.0*y0)
         colorbar, division=2, NCOLORS=!D.Table_size, /right, /vertical, $
         color=0, charsize=1.5, min=rmin, max=rmax, $
         /pscolor, POSITION=[x11,y11,x22,y22]
      endif

      DEVICE, /CLOSE
      SET_PLOT,'X'


      xdim=20000.0
      ydim=8000.0

      rimx=6.0
      rimy=6.0

      XSIZE=xdim/1000.0+rimx
      YSIZE=ydim/1000.0+rimy

      x0=500.0*rimx
      y0=500.0*rimy


      nx=y2-y1+1
      ny=200
      xs=y1*dy
      xu=y2*dy
      ys=pymin
      yu=pymax
      dpy=(yu-ys)/ny
      nf=fltarr(nx+1,ny+2)


      for l=0L,pnum-1 do begin
 
         xq=yp(l)
         yq=out_pyp(l)
 
         if (xs lt xq) and (xq lt xu) then begin
         if (ys lt yq) and (yq lt yu) then begin
 
         k1=floor(xq/dy)
         k2=floor(yq/dpy+0.5*(yu-ys)/dpy)
 
         xq=yp(l)/dy-k1
         yq=(out_pyp(l)+0.5*(yu-ys))/dpy-k2
 
         h1=(1.0-xq)*(1.0-yq)
         h2=xq*(1.0-yq)
         h3=(1.0-xq)*yq
         h4=xq*yq
 
         nf(k1,k2)=nf(k1,k2)+h1
         nf(k1+1,k2)=nf(k1+1,k2)+h2
         nf(k1,k2+1)=nf(k1,k2+1)+h3
         nf(k1+1,k2+1)=nf(k1+1,k2+1)+h4

         endif
         endif

      endfor


      fo='I_yyp-pz-image'$
      +'_t'+strtrim(string(tt),2)$
      +'.eps'
      print,fo


      SET_PLOT,'PS'
      DEVICE,/color,/encapsulated,file=fo,XSIZE=XSIZE,YSIZE=YSIZE
      erase, 255      


      sfxy=nf
      rmin=min(sfxy)
      rmax=0.5*max(sfxy)

      CONTOUR,sfxy(0:nx,0:ny), $
      xs=1,ys=1,charsize=1.6, $
      POSITION=[x0,y0,x0+xdim,y0+ydim], $
      xtitle='!8y/!8m', $
      xticks=2,$
      xr=[xs,xu], $
      ytitle='!8p!Dy!N/!8p!Dz',$
      yticks=2,$
      yr=[ys,yu], $
      background=255,color=0, $
      /DEVICE,/NODATA

      ixy = BYTSCL(sfxy(0:nx,0:ny), MIN=rmin, MAX=rmax)
      TV, ixy, !X.WINDOW(0),!Y.WINDOW(0),$
      XSIZE=!X.WINDOW(1) - !X.WINDOW(0),$
      YSIZE=!Y.WINDOW(1) - !Y.WINDOW(0), /NORM
;      oplot, liney, out_en_pypl, psym=3, color=300.0, thick=2.0

      if (rmax-rmin gt 0.0) then begin
         x11=(x0+xdim+500.0)/(xdim+2.0*x0)
         y11=y0/(ydim+2.0*y0)
         x22=(x0+xdim+800.0)/(xdim+2.0*x0)
         y22=(y0+ydim)/(ydim+2.0*y0)
         colorbar, division=2, NCOLORS=!D.Table_size, /right, /vertical, $
         color=0, charsize=1.5, min=rmin, max=rmax, $
         /pscolor, POSITION=[x11,y11,x22,y22]
      endif

      DEVICE, /CLOSE
      SET_PLOT,'X'


      xdim=8000.0
      ydim=8000.0

      rimx=6.0
      rimy=6.0

      XSIZE=xdim/1000.0+rimx
      YSIZE=ydim/1000.0+rimy

      x0=500.0*rimx
      y0=500.0*rimy


      nx=200
      ny=200
      xs=pxmin
      xu=pxmax
      ys=pymin
      yu=pymax
      dpx=(xu-xs)/nx
      dpy=(yu-ys)/ny
      nf=fltarr(nx+2,ny+2)


      for l=0L,pnum-1 do begin
 
         xq=out_pxp(l)
         yq=out_pyp(l)
 
         if (xs lt xq) and (xq lt xu) then begin
         if (ys lt yq) and (yq lt yu) then begin
 
         k1=floor(xq/dpx+0.5*(xu-xs)/dpx)
         k2=floor(yq/dpy+0.5*(yu-ys)/dpy)
 
         xq=(out_pxp(l)+0.5*(xu-xs))/dpx-k1
         yq=(out_pyp(l)+0.5*(yu-ys))/dpy-k2
 
         h1=(1.0-xq)*(1.0-yq)
         h2=xq*(1.0-yq)
         h3=(1.0-xq)*yq
         h4=xq*yq
 
         nf(k1,k2)=nf(k1,k2)+h1
         nf(k1+1,k2)=nf(k1+1,k2)+h2
         nf(k1,k2+1)=nf(k1,k2+1)+h3
         nf(k1+1,k2+1)=nf(k1+1,k2+1)+h4

         endif
         endif

      endfor


      fo='I_pxpy-pz-image'$
      +'_t'+strtrim(string(tt),2)$
      +'.eps'
      print,fo

      SET_PLOT,'PS'
      DEVICE,/color,/encapsulated,file=fo,XSIZE=XSIZE,YSIZE=YSIZE
      erase, 255      


      sfxy=nf
      rmin=min(sfxy)
      rmax=max(sfxy)
 
      CONTOUR,sfxy(0:nx,0:ny), $
      xs=1,ys=1,charsize=1.6, $
      POSITION=[x0,y0,x0+xdim,y0+ydim], $
      xtitle='!8p!Dx!N/!8p!Dz',$
      xticks=2,$
      xr=[xs,xu], $
      ytitle='!8p!Dy!N/!8p!Dz',$
      yticks=2,$
      yr=[ys,yu], $
      background=255,color=0, $
      /DEVICE,/NODATA

      ixy = BYTSCL(sfxy(0:nx,0:ny), MIN=rmin, MAX=rmax)
      TV, ixy, !X.WINDOW(0),!Y.WINDOW(0),$
      XSIZE=!X.WINDOW(1) - !X.WINDOW(0),$
      YSIZE=!Y.WINDOW(1) - !Y.WINDOW(0), /NORM
 
      if (rmax-rmin gt 0.0) then begin
         x11=(x0+xdim+500.0)/(xdim+2.0*x0)
         y11=y0/(ydim+2.0*y0)
         x22=(x0+xdim+800.0)/(xdim+2.0*x0)
         y22=(y0+ydim)/(ydim+2.0*y0)
         colorbar, division=2, NCOLORS=!D.Table_size, /right, /vertical, $
         color=0, charsize=1.5, min=rmin, max=rmax, $
         /pscolor, POSITION=[x11,y11,x22,y22]
      endif

      DEVICE, /CLOSE
      SET_PLOT,'X'

   endfor

end
