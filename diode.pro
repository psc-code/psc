; TIME STEPS AND PARTICLE NUMBER

nmax=30                  ; number of simulated time steps
pnum=3.0e7               ; number of quasi particles in simulation

; GRID SIZE, BOX LENGTHS, RESOLUTION IN m

x1=0                     ; x-grid definition
x2=1151                  ; x-grid definition
y1=0                     ; y-grid definition
y2=1151                  ; y-grid definition

lenx=280.0e-6            ; x-length of film on rear surface of diode
leny=280.0e-6            ; y-length of film on rear surface of diode
lenz=1.00e-9             ; z-length of film on rear surface of diode
vol=lenx*leny*lenz       ; total volume of ion film at rear surface

dx=lenx/(x2-x1+1)        ; x-step increment
dy=leny/(y2-y1+1)        ; y-step increment

; ION MOMENTUM RANGE IN m_ic

pxmin=-0.2               ; px momentum range for output plots
pxmax=0.2                ; px momentum range for output plots
pymin=-0.2               ; py momentum range for output plots
pymax=0.2                ; py momentum range for output plots

; GENERAL PARAMETERS

seed=3L                  ; random number generation
pi=3.1415926534          ; number pi
cc=3.0e8                 ; speed of light

; PARAMETERS FOR LASER ACCELERATED ION SHEATH

qqs1=1.6021e-19          ; ion charge
mps1=1.6700e-27          ; ion mass
tps1=20.0*qqs1           ; ion temperature in eV

ns1_hydro=6.65e28        ; density of ion film
ns1_limit=0.6*ns1_hydro    ; density of ion film at cutoff

rrs1=50.0e-6             ; FWHM of Gaussian ion sheath
lls1=5.00e-7             ; scale length of ion sheath profile
xxs1=140.0e-6            ; x-center of ion sheath
yys1=140.0e-6            ; y-center of ion sheath
amplx=0.0000             ; x-amplitude of flow perturbation 
amply=0.0006             ; y-amplitude of flow perturbation
amplz=2.1e5              ; z-amplitude of flow perturbation

; PARAMETERS FOR THE LASER ACCELERATED ION FLOW ENVELOPE

alpha0=35.0*pi/180.0     ; RCF-target rotation angle
rr0=100.0e-6             ; FWHM of Gaussian flow envelope
xx0=140.0e-6             ; x-center of flow envelope
yy0=140.0e-6             ; y-center of flow envelope
zz0=3.00e-6              ; location of rear surface of diode
zz1=120.0e-6             ; acceleration length of ions

cr=0.11                  ; strength of radial expansion of flow envelope
cz=0.12                  ; strength of longitudinal expansion of flow envelope

tt0=0.00e-12             ; start time of simulation
tt1=2.20e-12             ; end time of simulation
dtt=(tt1-tt0)/nmax       ; time increment


; DECLARATION OF VARIABLES 


      xpl=fltarr(pnum)
      ypl=fltarr(pnum)
      zpl=fltarr(pnum)

      pxpl=fltarr(pnum)
      pypl=fltarr(pnum)
      pzpl=fltarr(pnum)

      out_pxpl=fltarr(pnum)
      out_pypl=fltarr(pnum)


; INITIAL MOMENTUM AND ENVELOPE FLOW SETUP


      mm=0L
      for ll=0L,pnum-1 do begin

; spatial distribution

         ran1=randomn(seed,/UNIFORM)
         ran2=randomn(seed,/UNIFORM)
         ran3=randomn(seed,/UNIFORM)
         ran4=randomn(seed,/UNIFORM)
         ran5=randomn(seed,/UNIFORM)
         ran6=randomn(seed,/UNIFORM)

         if (ns1_limit/ns1_hydro le ran1) then begin

; initial spatial distribution of ion sheath

            pps1=rrs1*sqrt(-alog(ran2))
            phis1=6.2831853*ran3

            xq=xxs1+pps1*cos(phis1)
            yq=yys1+pps1*sin(phis1)
            hq=lls1*exp(-(pps1/rrs1)^2)
            zq=zz0-hq*alog(ran1)

; initial ion temperature in sheath

            pp=sqrt((1.0-tps1*alog(ran4)/(mps1*cc*cc))^2-1.0)
            theta=acos(2.0*ran5-1.0)
            phi=6.2831853*ran6

            pxq=pp*sin(theta)*cos(phi)
            pyq=pp*sin(theta)*sin(phi)
            pzq=pp*cos(theta)

; impose ion flow perturbations

            pxq=pxq+amplx*sin(2.0*pi*(xq-xx0)/3.5e-6)
            pyq=pyq+amply*sin(2.0*pi*(yq-yy0)/3.5e-6)
            pzq=pzq+amplz*(zq-zz0)

; define envelope of ion flow

            rr=(xq-xx0)^2+(yq-yy0)^2
            pxq=pxq+cr*(xq-xx0)*exp(-rr/rr0^2)*(zq-zz0)/rr0^2
            pyq=pyq+cr*(yq-yy0)*exp(-rr/rr0^2)*(zq-zz0)/rr0^2
            pzq=pzq+cz*alog(zq/zz0)/alog(zz1/zz0)

; initialization of arrays

            if (0.0e-6 le xq) and (xq le lenx) then begin
            if (0.0e-6 le yq) and (yq le leny) then begin
            if (0.0e-6 le pzq) and (pzq le 1.0e2) then begin

               xpl(mm)=xq
               ypl(mm)=yq
               zpl(mm)=zq
               pxpl(mm)=pxq
               pypl(mm)=pyq
               pzpl(mm)=pzq

               mm=mm+1

            endif
            endif
            endif

         endif
      endfor


; TIME EVOLUTION OF FLOW


      q_mass=ns1_hydro*vol/mm
      for nn=0L,nmax do begin

      tt=tt0+nn*dtt
      print,tt

      print,min(pxpl(0:mm-1)),max(pxpl(0:mm-1))
      print,min(pypl(0:mm-1)),max(pypl(0:mm-1))
      print,min(pzpl(0:mm-1)),max(pzpl(0:mm-1))

      for ll=0L,mm-1 do begin

         xq=xpl(ll)
         yq=ypl(ll)
         zq=zpl(ll)
         pxq=pxpl(ll)
         pyq=pypl(ll)
         pzq=pzpl(ll)

         rr=(xq-xx0)^2+(yq-yy0)^2
         dpxq=pxq-cr*(xq-xx0)*exp(-rr/rr0^2)*(zq-zz0)/rr0^2
         dpyq=pyq-cr*(yq-yy0)*exp(-rr/rr0^2)*(zq-zz0)/rr0^2
         dpzq=pzq-cz*alog(zq/zz0)/alog(zz1/zz0)

         pp=pxq*pxq+pyq*pyq+pzq*pzq
         xq=xq+cc*dtt*pxq/sqrt(1.0+pp)
         yq=yq+cc*dtt*pyq/sqrt(1.0+pp)
         zq=zq+cc*dtt*pzq/sqrt(1.0+pp)

         rr=(xq-xx0)^2+(yq-yy0)^2
         pxq=dpxq+cr*(xq-xx0)*exp(-rr/rr0^2)*(zq-zz0)/rr0^2
         pyq=dpyq+cr*(yq-yy0)*exp(-rr/rr0^2)*(zq-zz0)/rr0^2
         pzq=dpzq+cz*alog(zq/zz0)/alog(zz1/zz0)

         xpl(ll)=xq
         ypl(ll)=yq
         zpl(ll)=zq
         pxpl(ll)=pxq
         pypl(ll)=pyq
         pzpl(ll)=pzq

         out_pxpl(ll)=(+cos(alpha0)*pxq+sin(alpha0)*pyq)/pzq
         out_pypl(ll)=(-sin(alpha0)*pxq+cos(alpha0)*pyq)/pzq

      endfor


; DATA OUTPUT

; xpx/pz projection of phase space


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
      xx=fltarr(nx+2)
      yy=fltarr(ny+2)
      nf=fltarr(nx+1,ny+2)

      for l=0L,nx+1 do begin
         xx(l)=xs+l*dx
      endfor
      for l=0L,ny+1 do begin
         yy(l)=ys+l*dpx
      endfor

      for l=0L,mm-1 do begin
 
         xq=xpl(l)
         yq=out_pxpl(l)
 
         if (xs lt xq) and (xq lt xu) then begin
         if (ys lt yq) and (yq lt yu) then begin
 
         k1=floor(xq/dx)
         k2=floor(yq/dpx+0.5*(yu-ys)/dpx)
 
         xq=xpl(l)/dx-k1
         yq=(out_pxpl(l)+0.5*(yu-ys))/dpx-k2
 
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
      rmax=0.6*max(sfxy)
 
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
 
      if (rmax-rmin gt 0.0) then begin
         x11=(x0+xdim+500.0)/(xdim+2.0*x0)
         y11=y0/(ydim+2.0*y0)
         x22=(x0+xdim+800.0)/(xdim+2.0*x0)
         y22=(y0+ydim)/(ydim+2.0*y0)
         colorbar, division=2, NCOLORS=!D.Table_size, /right, /vertical, $
         color=0, charsize=1.5, min=rmin, max=rmax, $
         /pscolor, POSITION=[x11,y11,x22,y22]
      endif

      number=q_mass*mm
      ttt='!8t='+strtrim(string(tt*1.0e15),2)+'fs'
      number='!8N='+strtrim(string(number),2)
      xyouts,300,y0+ydim+2500.0,ttt,align=0.0,charsize=1.0,color=0,/DEVICE
      xyouts,300,y0+ydim+2200.0,number,align=0.0,charsize=1.0,color=0,/DEVICE

      DEVICE, /CLOSE
      SET_PLOT,'X'


; ypy/pz projection of phase space


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
      xx=fltarr(nx+2)
      yy=fltarr(ny+2)
      nf=fltarr(nx+1,ny+2)
      ng=fltarr(nx+1,ny+2)
      nh=fltarr(nx+1,ny+2)

      for l=0L,nx+1 do begin
         xx(l)=xs+l*dy
      endfor
      for l=0L,ny+1 do begin
         yy(l)=ys+l*dpy
      endfor

      pc1=0
      pc2=0
      pc3=0
      for l=0L,mm-1 do begin
 
         yo=ypl(l)
         pyo=out_pypl(l)
         pxq=pxpl(l)
         pyq=pypl(l)
         pzq=pzpl(l)
         enq=(mps1*cc*cc/qqs1)*(sqrt(1.0+pxq^2+pyq^2+pzq^2)-1.0)
 
         if (xs lt yo) and (yo lt xu) then begin
         if (ys lt pyo) and (pyo lt yu) then begin
 
         if (8.0e6 lt enq) and (enq lt 150.0e6) then begin

         pc1=pc1+1
         k1=floor(yo/dy)
         k2=floor(pyo/dpy+0.5*(yu-ys)/dpy)
 
         xq=yo/dy-k1
         yq=(pyo+0.5*(yu-ys))/dpy-k2
 
         h1=(1.0-xq)*(1.0-yq)
         h2=xq*(1.0-yq)
         h3=(1.0-xq)*yq
         h4=xq*yq
 
         ee=3.0e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            nf(k1,k2)=nf(k1,k2)+ww*h1
            nf(k1+1,k2)=nf(k1+1,k2)+ww*h2
            nf(k1,k2+1)=nf(k1,k2+1)+ww*h3
            nf(k1+1,k2+1)=nf(k1+1,k2+1)+ww*h4
         endif
         ee=3.5e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            nf(k1,k2)=nf(k1,k2)+ww*h1
            nf(k1+1,k2)=nf(k1+1,k2)+ww*h2
            nf(k1,k2+1)=nf(k1,k2+1)+ww*h3
            nf(k1+1,k2+1)=nf(k1+1,k2+1)+ww*h4
         endif
         ee=6.0e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            nf(k1,k2)=nf(k1,k2)+ww*h1
            nf(k1+1,k2)=nf(k1+1,k2)+ww*h2
            nf(k1,k2+1)=nf(k1,k2+1)+ww*h3
            nf(k1+1,k2+1)=nf(k1+1,k2+1)+ww*h4
         endif
         ee=6.5e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            nf(k1,k2)=nf(k1,k2)+ww*h1
            nf(k1+1,k2)=nf(k1+1,k2)+ww*h2
            nf(k1,k2+1)=nf(k1,k2+1)+ww*h3
            nf(k1+1,k2+1)=nf(k1+1,k2+1)+ww*h4
         endif
         ee=8.0e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            nf(k1,k2)=nf(k1,k2)+ww*h1
            nf(k1+1,k2)=nf(k1+1,k2)+ww*h2
            nf(k1,k2+1)=nf(k1,k2+1)+ww*h3
            nf(k1+1,k2+1)=nf(k1+1,k2+1)+ww*h4
         endif
         ee=8.5e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            nf(k1,k2)=nf(k1,k2)+ww*h1
            nf(k1+1,k2)=nf(k1+1,k2)+ww*h2
            nf(k1,k2+1)=nf(k1,k2+1)+ww*h3
            nf(k1+1,k2+1)=nf(k1+1,k2+1)+ww*h4
         endif
         ee=11.0e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            nf(k1,k2)=nf(k1,k2)+ww*h1
            nf(k1+1,k2)=nf(k1+1,k2)+ww*h2
            nf(k1,k2+1)=nf(k1,k2+1)+ww*h3
            nf(k1+1,k2+1)=nf(k1+1,k2+1)+ww*h4
         endif
         ee=11.5e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            nf(k1,k2)=nf(k1,k2)+ww*h1
            nf(k1+1,k2)=nf(k1+1,k2)+ww*h2
            nf(k1,k2+1)=nf(k1,k2+1)+ww*h3
            nf(k1+1,k2+1)=nf(k1+1,k2+1)+ww*h4
         endif
         ee=12.0e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            nf(k1,k2)=nf(k1,k2)+ww*h1
            nf(k1+1,k2)=nf(k1+1,k2)+ww*h2
            nf(k1,k2+1)=nf(k1,k2+1)+ww*h3
            nf(k1+1,k2+1)=nf(k1+1,k2+1)+ww*h4
         endif
         ee=12.5e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            nf(k1,k2)=nf(k1,k2)+ww*h1
            nf(k1+1,k2)=nf(k1+1,k2)+ww*h2
            nf(k1,k2+1)=nf(k1,k2+1)+ww*h3
            nf(k1+1,k2+1)=nf(k1+1,k2+1)+ww*h4
         endif
         endif


         if (11.0e6 lt enq) and (enq lt 150.0e6) then begin

         pc2=pc2+1 
         k1=floor(yo/dy)
         k2=floor(pyo/dpy+0.5*(yu-ys)/dpy)
 
         xq=yo/dy-k1
         yq=(pyo+0.5*(yu-ys))/dpy-k2
 
         h1=(1.0-xq)*(1.0-yq)
         h2=xq*(1.0-yq)
         h3=(1.0-xq)*yq
         h4=xq*yq
 
         ee=3.0e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            ng(k1,k2)=ng(k1,k2)+ww*h1
            ng(k1+1,k2)=ng(k1+1,k2)+ww*h2
            ng(k1,k2+1)=ng(k1,k2+1)+ww*h3
            ng(k1+1,k2+1)=ng(k1+1,k2+1)+ww*h4
         endif
         ee=3.5e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            ng(k1,k2)=ng(k1,k2)+ww*h1
            ng(k1+1,k2)=ng(k1+1,k2)+ww*h2
            ng(k1,k2+1)=ng(k1,k2+1)+ww*h3
            ng(k1+1,k2+1)=ng(k1+1,k2+1)+ww*h4
         endif
         ee=6.0e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            ng(k1,k2)=ng(k1,k2)+ww*h1
            ng(k1+1,k2)=ng(k1+1,k2)+ww*h2
            ng(k1,k2+1)=ng(k1,k2+1)+ww*h3
            ng(k1+1,k2+1)=ng(k1+1,k2+1)+ww*h4
         endif
         ee=6.5e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            ng(k1,k2)=ng(k1,k2)+ww*h1
            ng(k1+1,k2)=ng(k1+1,k2)+ww*h2
            ng(k1,k2+1)=ng(k1,k2+1)+ww*h3
            ng(k1+1,k2+1)=ng(k1+1,k2+1)+ww*h4
         endif
         ee=8.0e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            ng(k1,k2)=ng(k1,k2)+ww*h1
            ng(k1+1,k2)=ng(k1+1,k2)+ww*h2
            ng(k1,k2+1)=ng(k1,k2+1)+ww*h3
            ng(k1+1,k2+1)=ng(k1+1,k2+1)+ww*h4
         endif
         ee=8.5e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            ng(k1,k2)=ng(k1,k2)+ww*h1
            ng(k1+1,k2)=ng(k1+1,k2)+ww*h2
            ng(k1,k2+1)=ng(k1,k2+1)+ww*h3
            ng(k1+1,k2+1)=ng(k1+1,k2+1)+ww*h4
         endif
         ee=11.0e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            ng(k1,k2)=ng(k1,k2)+ww*h1
            ng(k1+1,k2)=ng(k1+1,k2)+ww*h2
            ng(k1,k2+1)=ng(k1,k2+1)+ww*h3
            ng(k1+1,k2+1)=ng(k1+1,k2+1)+ww*h4
         endif
         ee=11.5e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            ng(k1,k2)=ng(k1,k2)+ww*h1
            ng(k1+1,k2)=ng(k1+1,k2)+ww*h2
            ng(k1,k2+1)=ng(k1,k2+1)+ww*h3
            ng(k1+1,k2+1)=ng(k1+1,k2+1)+ww*h4
         endif
         ee=12.0e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            ng(k1,k2)=ng(k1,k2)+ww*h1
            ng(k1+1,k2)=ng(k1+1,k2)+ww*h2
            ng(k1,k2+1)=ng(k1,k2+1)+ww*h3
            ng(k1+1,k2+1)=ng(k1+1,k2+1)+ww*h4
         endif
         ee=12.5e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            ng(k1,k2)=ng(k1,k2)+ww*h1
            ng(k1+1,k2)=ng(k1+1,k2)+ww*h2
            ng(k1,k2+1)=ng(k1,k2+1)+ww*h3
            ng(k1+1,k2+1)=ng(k1+1,k2+1)+ww*h4
         endif
         endif


         if (12.0e6 lt enq) and (enq lt 150.0e6) then begin

         pc3=pc3+1 
         k1=floor(yo/dy)
         k2=floor(pyo/dpy+0.5*(yu-ys)/dpy)
 
         xq=yo/dy-k1
         yq=(pyo+0.5*(yu-ys))/dpy-k2
 
         h1=(1.0-xq)*(1.0-yq)
         h2=xq*(1.0-yq)
         h3=(1.0-xq)*yq
         h4=xq*yq
 
         ee=3.0e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            nh(k1,k2)=nh(k1,k2)+ww*h1
            nh(k1+1,k2)=nh(k1+1,k2)+ww*h2
            nh(k1,k2+1)=nh(k1,k2+1)+ww*h3
            nh(k1+1,k2+1)=nh(k1+1,k2+1)+ww*h4
         endif
         ee=3.5e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            nh(k1,k2)=nh(k1,k2)+ww*h1
            nh(k1+1,k2)=nh(k1+1,k2)+ww*h2
            nh(k1,k2+1)=nh(k1,k2+1)+ww*h3
            nh(k1+1,k2+1)=nh(k1+1,k2+1)+ww*h4
         endif
         ee=6.0e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            nh(k1,k2)=nh(k1,k2)+ww*h1
            nh(k1+1,k2)=nh(k1+1,k2)+ww*h2
            nh(k1,k2+1)=nh(k1,k2+1)+ww*h3
            nh(k1+1,k2+1)=nh(k1+1,k2+1)+ww*h4
         endif
         ee=6.5e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            nh(k1,k2)=nh(k1,k2)+ww*h1
            nh(k1+1,k2)=nh(k1+1,k2)+ww*h2
            nh(k1,k2+1)=nh(k1,k2+1)+ww*h3
            nh(k1+1,k2+1)=nh(k1+1,k2+1)+ww*h4
         endif
         ee=8.0e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            nh(k1,k2)=nh(k1,k2)+ww*h1
            nh(k1+1,k2)=nh(k1+1,k2)+ww*h2
            nh(k1,k2+1)=nh(k1,k2+1)+ww*h3
            nh(k1+1,k2+1)=nh(k1+1,k2+1)+ww*h4
         endif
         ee=8.5e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            nh(k1,k2)=nh(k1,k2)+ww*h1
            nh(k1+1,k2)=nh(k1+1,k2)+ww*h2
            nh(k1,k2+1)=nh(k1,k2+1)+ww*h3
            nh(k1+1,k2+1)=nh(k1+1,k2+1)+ww*h4
         endif
         ee=11.0e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            nh(k1,k2)=nh(k1,k2)+ww*h1
            nh(k1+1,k2)=nh(k1+1,k2)+ww*h2
            nh(k1,k2+1)=nh(k1,k2+1)+ww*h3
            nh(k1+1,k2+1)=nh(k1+1,k2+1)+ww*h4
         endif
         ee=11.5e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            nh(k1,k2)=nh(k1,k2)+ww*h1
            nh(k1+1,k2)=nh(k1+1,k2)+ww*h2
            nh(k1,k2+1)=nh(k1,k2+1)+ww*h3
            nh(k1+1,k2+1)=nh(k1+1,k2+1)+ww*h4
         endif
         ee=12.0e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            nh(k1,k2)=nh(k1,k2)+ww*h1
            nh(k1+1,k2)=nh(k1+1,k2)+ww*h2
            nh(k1,k2+1)=nh(k1,k2+1)+ww*h3
            nh(k1+1,k2+1)=nh(k1+1,k2+1)+ww*h4
         endif
         ee=12.5e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            nh(k1,k2)=nh(k1,k2)+ww*h1
            nh(k1+1,k2)=nh(k1+1,k2)+ww*h2
            nh(k1,k2+1)=nh(k1,k2+1)+ww*h3
            nh(k1+1,k2+1)=nh(k1+1,k2+1)+ww*h4
         endif
         endif

         endif
         endif

      endfor


      fo='I_yyp-pz-image_8MeV'$
      +'_t'+strtrim(string(tt),2)$
      +'.eps'
      print,fo


      SET_PLOT,'PS'
      DEVICE,/color,/encapsulated,file=fo,XSIZE=XSIZE,YSIZE=YSIZE
      erase, 255      


      sfxy=nf
      rmin=min(sfxy)
      rmax=0.6*max(sfxy)

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

      if (rmax-rmin gt 0.0) then begin
         x11=(x0+xdim+500.0)/(xdim+2.0*x0)
         y11=y0/(ydim+2.0*y0)
         x22=(x0+xdim+800.0)/(xdim+2.0*x0)
         y22=(y0+ydim)/(ydim+2.0*y0)
         colorbar, division=2, NCOLORS=!D.Table_size, /right, /vertical, $
         color=0, charsize=1.5, min=rmin, max=rmax, $
         /pscolor, POSITION=[x11,y11,x22,y22]
      endif

      number=q_mass*pc1
      ttt='!8t='+strtrim(string(tt*1.0e15),2)+'fs'
      number='!8N='+strtrim(string(number),2)
      xyouts,300,y0+ydim+2500.0,ttt,align=0.0,charsize=1.0,color=0,/DEVICE
      xyouts,300,y0+ydim+2200.0,number,align=0.0,charsize=1.0,color=0,/DEVICE

      DEVICE, /CLOSE
      SET_PLOT,'X'


      fo='I_yyp-pz-image_11MeV'$
      +'_t'+strtrim(string(tt),2)$
      +'.eps'
      print,fo


      SET_PLOT,'PS'
      DEVICE,/color,/encapsulated,file=fo,XSIZE=XSIZE,YSIZE=YSIZE
      erase, 255      


      sfxy=ng
      rmin=min(sfxy)
      rmax=0.6*max(sfxy)

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

      if (rmax-rmin gt 0.0) then begin
         x11=(x0+xdim+500.0)/(xdim+2.0*x0)
         y11=y0/(ydim+2.0*y0)
         x22=(x0+xdim+800.0)/(xdim+2.0*x0)
         y22=(y0+ydim)/(ydim+2.0*y0)
         colorbar, division=2, NCOLORS=!D.Table_size, /right, /vertical, $
         color=0, charsize=1.5, min=rmin, max=rmax, $
         /pscolor, POSITION=[x11,y11,x22,y22]
      endif

      number=q_mass*pc2
      ttt='!8t='+strtrim(string(tt*1.0e15),2)+'fs'
      number='!8N='+strtrim(string(number),2)
      xyouts,300,y0+ydim+2500.0,ttt,align=0.0,charsize=1.0,color=0,/DEVICE
      xyouts,300,y0+ydim+2200.0,number,align=0.0,charsize=1.0,color=0,/DEVICE

      DEVICE, /CLOSE
      SET_PLOT,'X'


      fo='I_yyp-pz-image_12MeV'$
      +'_t'+strtrim(string(tt),2)$
      +'.eps'
      print,fo


      SET_PLOT,'PS'
      DEVICE,/color,/encapsulated,file=fo,XSIZE=XSIZE,YSIZE=YSIZE
      erase, 255      


      sfxy=nh
      rmin=min(sfxy)
      rmax=0.6*max(sfxy)

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

      if (rmax-rmin gt 0.0) then begin
         x11=(x0+xdim+500.0)/(xdim+2.0*x0)
         y11=y0/(ydim+2.0*y0)
         x22=(x0+xdim+800.0)/(xdim+2.0*x0)
         y22=(y0+ydim)/(ydim+2.0*y0)
         colorbar, division=2, NCOLORS=!D.Table_size, /right, /vertical, $
         color=0, charsize=1.5, min=rmin, max=rmax, $
         /pscolor, POSITION=[x11,y11,x22,y22]
      endif

      number=q_mass*pc3
      ttt='!8t='+strtrim(string(tt*1.0e15),2)+'fs'
      number='!8N='+strtrim(string(number),2)
      xyouts,300,y0+ydim+2500.0,ttt,align=0.0,charsize=1.0,color=0,/DEVICE
      xyouts,300,y0+ydim+2200.0,number,align=0.0,charsize=1.0,color=0,/DEVICE

      DEVICE, /CLOSE
      SET_PLOT,'X'


; pxpy/pz projection of phase space


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
      xx=fltarr(nx+2)
      yy=fltarr(ny+2)
      nf=fltarr(nx+2,ny+2)
      ng=fltarr(nx+2,ny+2)
      nh=fltarr(nx+2,ny+2)

      for l=0L,nx+1 do begin
         xx(l)=xs+l*dpx
      endfor
      for l=0L,ny+1 do begin
         yy(l)=ys+l*dpy
      endfor

      pc1=0
      pc2=0
      pc3=0
      for l=0L,mm-1 do begin
 
         pxo=out_pxpl(l)
         pyo=out_pypl(l)
         pxq=pxpl(l)
         pyq=pypl(l)
         pzq=pzpl(l)
         enq=(mps1*cc*cc/qqs1)*(sqrt(1.0+pxq^2+pyq^2+pzq^2)-1.0)

         if (xs lt pxo) and (pxo lt xu) then begin
         if (ys lt pyo) and (pyo lt yu) then begin

         if (8.0e6 lt enq) and (enq lt 150.0e6) then begin

         pc1=pc1+1
         k1=floor(pxo/dpx+0.5*(xu-xs)/dpx)
         k2=floor(pyo/dpy+0.5*(yu-ys)/dpy)
 
         xq=(pxo+0.5*(xu-xs))/dpx-k1
         yq=(pyo+0.5*(yu-ys))/dpy-k2
 
         h1=(1.0-xq)*(1.0-yq)
         h2=xq*(1.0-yq)
         h3=(1.0-xq)*yq
         h4=xq*yq
 
         ee=3.0e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            nf(k1,k2)=nf(k1,k2)+ww*h1
            nf(k1+1,k2)=nf(k1+1,k2)+ww*h2
            nf(k1,k2+1)=nf(k1,k2+1)+ww*h3
            nf(k1+1,k2+1)=nf(k1+1,k2+1)+ww*h4
         endif
         ee=3.5e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            nf(k1,k2)=nf(k1,k2)+ww*h1
            nf(k1+1,k2)=nf(k1+1,k2)+ww*h2
            nf(k1,k2+1)=nf(k1,k2+1)+ww*h3
            nf(k1+1,k2+1)=nf(k1+1,k2+1)+ww*h4
         endif
         ee=6.0e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            nf(k1,k2)=nf(k1,k2)+ww*h1
            nf(k1+1,k2)=nf(k1+1,k2)+ww*h2
            nf(k1,k2+1)=nf(k1,k2+1)+ww*h3
            nf(k1+1,k2+1)=nf(k1+1,k2+1)+ww*h4
         endif
         ee=6.5e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            nf(k1,k2)=nf(k1,k2)+ww*h1
            nf(k1+1,k2)=nf(k1+1,k2)+ww*h2
            nf(k1,k2+1)=nf(k1,k2+1)+ww*h3
            nf(k1+1,k2+1)=nf(k1+1,k2+1)+ww*h4
         endif
         ee=8.0e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            nf(k1,k2)=nf(k1,k2)+ww*h1
            nf(k1+1,k2)=nf(k1+1,k2)+ww*h2
            nf(k1,k2+1)=nf(k1,k2+1)+ww*h3
            nf(k1+1,k2+1)=nf(k1+1,k2+1)+ww*h4
         endif
         ee=8.5e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            nf(k1,k2)=nf(k1,k2)+ww*h1
            nf(k1+1,k2)=nf(k1+1,k2)+ww*h2
            nf(k1,k2+1)=nf(k1,k2+1)+ww*h3
            nf(k1+1,k2+1)=nf(k1+1,k2+1)+ww*h4
         endif
         ee=11.0e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            nf(k1,k2)=nf(k1,k2)+ww*h1
            nf(k1+1,k2)=nf(k1+1,k2)+ww*h2
            nf(k1,k2+1)=nf(k1,k2+1)+ww*h3
            nf(k1+1,k2+1)=nf(k1+1,k2+1)+ww*h4
         endif
         ee=11.5e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            nf(k1,k2)=nf(k1,k2)+ww*h1
            nf(k1+1,k2)=nf(k1+1,k2)+ww*h2
            nf(k1,k2+1)=nf(k1,k2+1)+ww*h3
            nf(k1+1,k2+1)=nf(k1+1,k2+1)+ww*h4
         endif
         ee=12.0e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            nf(k1,k2)=nf(k1,k2)+ww*h1
            nf(k1+1,k2)=nf(k1+1,k2)+ww*h2
            nf(k1,k2+1)=nf(k1,k2+1)+ww*h3
            nf(k1+1,k2+1)=nf(k1+1,k2+1)+ww*h4
         endif
         ee=12.5e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            nf(k1,k2)=nf(k1,k2)+ww*h1
            nf(k1+1,k2)=nf(k1+1,k2)+ww*h2
            nf(k1,k2+1)=nf(k1,k2+1)+ww*h3
            nf(k1+1,k2+1)=nf(k1+1,k2+1)+ww*h4
         endif
         endif


         if (11.0e6 lt enq) and (enq lt 150.0e6) then begin

         pc2=pc2+1 
         k1=floor(pxo/dpx+0.5*(xu-xs)/dpx)
         k2=floor(pyo/dpy+0.5*(yu-ys)/dpy)
 
         xq=(pxo+0.5*(xu-xs))/dpx-k1
         yq=(pyo+0.5*(yu-ys))/dpy-k2
 
         h1=(1.0-xq)*(1.0-yq)
         h2=xq*(1.0-yq)
         h3=(1.0-xq)*yq
         h4=xq*yq
 
         ee=3.0e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            ng(k1,k2)=ng(k1,k2)+ww*h1
            ng(k1+1,k2)=ng(k1+1,k2)+ww*h2
            ng(k1,k2+1)=ng(k1,k2+1)+ww*h3
            ng(k1+1,k2+1)=ng(k1+1,k2+1)+ww*h4
         endif
         ee=3.5e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            ng(k1,k2)=ng(k1,k2)+ww*h1
            ng(k1+1,k2)=ng(k1+1,k2)+ww*h2
            ng(k1,k2+1)=ng(k1,k2+1)+ww*h3
            ng(k1+1,k2+1)=ng(k1+1,k2+1)+ww*h4
         endif
         ee=6.0e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            ng(k1,k2)=ng(k1,k2)+ww*h1
            ng(k1+1,k2)=ng(k1+1,k2)+ww*h2
            ng(k1,k2+1)=ng(k1,k2+1)+ww*h3
            ng(k1+1,k2+1)=ng(k1+1,k2+1)+ww*h4
         endif
         ee=6.5e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            ng(k1,k2)=ng(k1,k2)+ww*h1
            ng(k1+1,k2)=ng(k1+1,k2)+ww*h2
            ng(k1,k2+1)=ng(k1,k2+1)+ww*h3
            ng(k1+1,k2+1)=ng(k1+1,k2+1)+ww*h4
         endif
         ee=8.0e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            ng(k1,k2)=ng(k1,k2)+ww*h1
            ng(k1+1,k2)=ng(k1+1,k2)+ww*h2
            ng(k1,k2+1)=ng(k1,k2+1)+ww*h3
            ng(k1+1,k2+1)=ng(k1+1,k2+1)+ww*h4
         endif
         ee=8.5e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            ng(k1,k2)=ng(k1,k2)+ww*h1
            ng(k1+1,k2)=ng(k1+1,k2)+ww*h2
            ng(k1,k2+1)=ng(k1,k2+1)+ww*h3
            ng(k1+1,k2+1)=ng(k1+1,k2+1)+ww*h4
         endif
         ee=11.0e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            ng(k1,k2)=ng(k1,k2)+ww*h1
            ng(k1+1,k2)=ng(k1+1,k2)+ww*h2
            ng(k1,k2+1)=ng(k1,k2+1)+ww*h3
            ng(k1+1,k2+1)=ng(k1+1,k2+1)+ww*h4
         endif
         ee=11.5e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            ng(k1,k2)=ng(k1,k2)+ww*h1
            ng(k1+1,k2)=ng(k1+1,k2)+ww*h2
            ng(k1,k2+1)=ng(k1,k2+1)+ww*h3
            ng(k1+1,k2+1)=ng(k1+1,k2+1)+ww*h4
         endif
         ee=12.0e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            ng(k1,k2)=ng(k1,k2)+ww*h1
            ng(k1+1,k2)=ng(k1+1,k2)+ww*h2
            ng(k1,k2+1)=ng(k1,k2+1)+ww*h3
            ng(k1+1,k2+1)=ng(k1+1,k2+1)+ww*h4
         endif
         ee=12.5e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            ng(k1,k2)=ng(k1,k2)+ww*h1
            ng(k1+1,k2)=ng(k1+1,k2)+ww*h2
            ng(k1,k2+1)=ng(k1,k2+1)+ww*h3
            ng(k1+1,k2+1)=ng(k1+1,k2+1)+ww*h4
         endif
         endif


         if (12.0e6 lt enq) and (enq lt 150.0e6) then begin

         pc3=pc3+1 
         k1=floor(pxo/dpx+0.5*(xu-xs)/dpx)
         k2=floor(pyo/dpy+0.5*(yu-ys)/dpy)
 
         xq=(pxo+0.5*(xu-xs))/dpx-k1
         yq=(pyo+0.5*(yu-ys))/dpy-k2
 
         h1=(1.0-xq)*(1.0-yq)
         h2=xq*(1.0-yq)
         h3=(1.0-xq)*yq
         h4=xq*yq
 
         ee=3.0e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            nh(k1,k2)=nh(k1,k2)+ww*h1
            nh(k1+1,k2)=nh(k1+1,k2)+ww*h2
            nh(k1,k2+1)=nh(k1,k2+1)+ww*h3
            nh(k1+1,k2+1)=nh(k1+1,k2+1)+ww*h4
         endif
         ee=3.5e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            nh(k1,k2)=nh(k1,k2)+ww*h1
            nh(k1+1,k2)=nh(k1+1,k2)+ww*h2
            nh(k1,k2+1)=nh(k1,k2+1)+ww*h3
            nh(k1+1,k2+1)=nh(k1+1,k2+1)+ww*h4
         endif
         ee=6.0e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            nh(k1,k2)=nh(k1,k2)+ww*h1
            nh(k1+1,k2)=nh(k1+1,k2)+ww*h2
            nh(k1,k2+1)=nh(k1,k2+1)+ww*h3
            nh(k1+1,k2+1)=nh(k1+1,k2+1)+ww*h4
         endif
         ee=6.5e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            nh(k1,k2)=nh(k1,k2)+ww*h1
            nh(k1+1,k2)=nh(k1+1,k2)+ww*h2
            nh(k1,k2+1)=nh(k1,k2+1)+ww*h3
            nh(k1+1,k2+1)=nh(k1+1,k2+1)+ww*h4
         endif
         ee=8.0e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            nh(k1,k2)=nh(k1,k2)+ww*h1
            nh(k1+1,k2)=nh(k1+1,k2)+ww*h2
            nh(k1,k2+1)=nh(k1,k2+1)+ww*h3
            nh(k1+1,k2+1)=nh(k1+1,k2+1)+ww*h4
         endif
         ee=8.5e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            nh(k1,k2)=nh(k1,k2)+ww*h1
            nh(k1+1,k2)=nh(k1+1,k2)+ww*h2
            nh(k1,k2+1)=nh(k1,k2+1)+ww*h3
            nh(k1+1,k2+1)=nh(k1+1,k2+1)+ww*h4
         endif
         ee=11.0e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            nh(k1,k2)=nh(k1,k2)+ww*h1
            nh(k1+1,k2)=nh(k1+1,k2)+ww*h2
            nh(k1,k2+1)=nh(k1,k2+1)+ww*h3
            nh(k1+1,k2+1)=nh(k1+1,k2+1)+ww*h4
         endif
         ee=11.5e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            nh(k1,k2)=nh(k1,k2)+ww*h1
            nh(k1+1,k2)=nh(k1+1,k2)+ww*h2
            nh(k1,k2+1)=nh(k1,k2+1)+ww*h3
            nh(k1+1,k2+1)=nh(k1+1,k2+1)+ww*h4
         endif
         ee=12.0e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            nh(k1,k2)=nh(k1,k2)+ww*h1
            nh(k1+1,k2)=nh(k1+1,k2)+ww*h2
            nh(k1,k2+1)=nh(k1,k2+1)+ww*h3
            nh(k1+1,k2+1)=nh(k1+1,k2+1)+ww*h4
         endif
         ee=12.5e6
         dee=3.0e6
         if (enq gt ee) then begin
            ww=0.1+0.9*exp(-(enq-ee)/dee)
            nh(k1,k2)=nh(k1,k2)+ww*h1
            nh(k1+1,k2)=nh(k1+1,k2)+ww*h2
            nh(k1,k2+1)=nh(k1,k2+1)+ww*h3
            nh(k1+1,k2+1)=nh(k1+1,k2+1)+ww*h4
         endif
         endif

         endif
         endif

      endfor


      fo='I_pxpy-pz-image_8MeV'$
      +'_t'+strtrim(string(tt),2)$
      +'.eps'
      print,fo

      SET_PLOT,'PS'
      DEVICE,/color,/encapsulated,file=fo,XSIZE=XSIZE,YSIZE=YSIZE
      erase, 255      


      sfxy=nf
      rmin=min(sfxy)
      rmax=0.8*max(sfxy)
 
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

      number=q_mass*pc1
      ttt='!8t='+strtrim(string(tt*1.0e15),2)+'fs'
      number='!8N='+strtrim(string(number),2)
      xyouts,300,y0+ydim+2500.0,ttt,align=0.0,charsize=1.0,color=0,/DEVICE
      xyouts,300,y0+ydim+2200.0,number,align=0.0,charsize=1.0,color=0,/DEVICE

      DEVICE, /CLOSE
      SET_PLOT,'X'


      fo='I_pxpy-pz-image_11MeV'$
      +'_t'+strtrim(string(tt),2)$
      +'.eps'
      print,fo

      SET_PLOT,'PS'
      DEVICE,/color,/encapsulated,file=fo,XSIZE=XSIZE,YSIZE=YSIZE
      erase, 255      


      sfxy=ng
      rmin=min(sfxy)
      rmax=0.8*max(sfxy)
 
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

      number=q_mass*pc2
      ttt='!8t='+strtrim(string(tt*1.0e15),2)+'fs'
      number='!8N='+strtrim(string(number),2)
      xyouts,300,y0+ydim+2500.0,ttt,align=0.0,charsize=1.0,color=0,/DEVICE
      xyouts,300,y0+ydim+2200.0,number,align=0.0,charsize=1.0,color=0,/DEVICE

      DEVICE, /CLOSE
      SET_PLOT,'X'


      fo='I_pxpy-pz-image_12MeV'$
      +'_t'+strtrim(string(tt),2)$
      +'.eps'
      print,fo

      SET_PLOT,'PS'
      DEVICE,/color,/encapsulated,file=fo,XSIZE=XSIZE,YSIZE=YSIZE
      erase, 255      


      sfxy=nh
      rmin=min(sfxy)
      rmax=0.8*max(sfxy)
 
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

      number=q_mass*pc3
      ttt='!8t='+strtrim(string(tt*1.0e15),2)+'fs'
      number='!8N='+strtrim(string(number),2)
      xyouts,300,y0+ydim+2500.0,ttt,align=0.0,charsize=1.0,color=0,/DEVICE
      xyouts,300,y0+ydim+2200.0,number,align=0.0,charsize=1.0,color=0,/DEVICE

      DEVICE, /CLOSE
      SET_PLOT,'X'


; yz projection of phase space


      xdim=20000.0
      ydim=70.0*20000.0/140.0

      rimx=6.0
      rimy=6.0

      XSIZE=xdim/1000.0+rimx
      YSIZE=ydim/1000.0+rimy

      x0=500.0*rimx
      y0=500.0*rimy


      fo='I_yz-image'$
      +'_t'+strtrim(string(tt),2)$
      +'.eps'
      print,fo


      SET_PLOT,'PS'
      DEVICE,/color,/encapsulated,file=fo,XSIZE=XSIZE,YSIZE=YSIZE
      erase, 255      

      TVLCT, [0,255,0,0], [0,0,255,0], [0,0,0,255]
      plot, 1.0e6*ypl(0:mm-1), 1.0e6*zpl(0:mm-1), xs=1, ys=1, charsize=1.6, psym=3, $
      POSITION=[x0,y0,x0+xdim,y0+ydim],/DEVICE, $ 
      xtitle='!8y/!7l!8m',$
      xticks=2,$
      xr=1.0e6*[0.0,leny], $
      ytitle='!8z/!7l!8m',$
      yticks=2,$
      yr=1.0e6*[0.0,70.0e-6], $
      color=0,background=255 
;      oplot, ypl2, zpl2, psym=3, color=1
;      oplot, ypl1, zpl1, psym=3, color=2
;      oplot, ypl0, zpl0, psym=3, color=3

      number=q_mass*mm
      ttt='!8t='+strtrim(string(tt*1.0e15),2)+'fs'
      number='!8N='+strtrim(string(number),2)
      xyouts,300,y0+ydim+2500.0,ttt,align=0.0,charsize=1.0,color=0,/DEVICE
      xyouts,300,y0+ydim+2200.0,number,align=0.0,charsize=1.0,color=0,/DEVICE

      DEVICE, /CLOSE
      SET_PLOT,'X'


; pyEz/pz projection of phase space


      xdim=8000.0
      ydim=8000.0

      rimx=6.0
      rimy=6.0

      XSIZE=xdim/1000.0+rimx
      YSIZE=ydim/1000.0+rimy

      x0=500.0*rimx
      y0=500.0*rimy


      fo='I_pyEz-pz-image'$
      +'_t'+strtrim(string(tt),2)$
      +'.eps'
      print,fo


      SET_PLOT,'PS'
      DEVICE,/color,/encapsulated,file=fo,XSIZE=XSIZE,YSIZE=YSIZE
      erase, 255      

      TVLCT, [0,255,0,0], [0,0,255,0], [0,0,0,255]
      plot, out_pypl(0:mm-1), 469.8*pzpl(0:mm-1)*pzpl(0:mm-1), $
      xs=1, ys=1, charsize=1.6, psym=3, $
      POSITION=[x0,y0,x0+xdim,y0+ydim],/DEVICE, $ 
      xtitle='!8p!Dy!N/!8p!Dz',$
      xticks=2,$
      xr=[pymin,pymax], $
      ytitle='!8E!Dz!N/!8MeV',$
      yticks=2,$
      yr=[8.00,12.0], $
      color=0,background=255 

      number=q_mass*mm
      ttt='!8t='+strtrim(string(tt*1.0e15),2)+'fs'
      number='!8N='+strtrim(string(number),2)
      xyouts,300,y0+ydim+2500.0,ttt,align=0.0,charsize=1.0,color=0,/DEVICE
      xyouts,300,y0+ydim+2200.0,number,align=0.0,charsize=1.0,color=0,/DEVICE

      DEVICE, /CLOSE
      SET_PLOT,'X'


; zpz projection of phase space


      xdim=12000.0
      ydim=8000.0

      rimx=6.0
      rimy=6.0

      XSIZE=xdim/1000.0+rimx
      YSIZE=ydim/1000.0+rimy

      x0=500.0*rimx
      y0=500.0*rimy


      fo='I_zpz-image'$
      +'_t'+strtrim(string(tt),2)$
      +'.eps'
      print,fo


      SET_PLOT,'PS'
      DEVICE,/color,/encapsulated,file=fo,XSIZE=XSIZE,YSIZE=YSIZE
      erase, 255      

      TVLCT, [0,255,0,0], [0,0,255,0], [0,0,0,255]
      plot, zpl(0:mm-1), pzpl(0:mm-1), xs=1, ys=1, charsize=1.6, psym=3, $
      POSITION=[x0,y0,x0+xdim,y0+ydim],/DEVICE, $ 
      xtitle='!8z/!7l!8m',$
      xticks=2,$
      xr=[0.0,70.0e-6], $
      ytitle='!8p!Dz!N/!8mc',$
      yticks=2,$
      yr=[0.00,0.15], $
      color=0,background=255 

      number=q_mass*mm
      ttt='!8t='+strtrim(string(tt*1.0e15),2)+'fs'
      number='!8N='+strtrim(string(number),2)
      xyouts,300,y0+ydim+2500.0,ttt,align=0.0,charsize=1.0,color=0,/DEVICE
      xyouts,300,y0+ydim+2200.0,number,align=0.0,charsize=1.0,color=0,/DEVICE

      DEVICE, /CLOSE
      SET_PLOT,'X'


; zpy projection of phase space


      xdim=12000.0
      ydim=8000.0

      rimx=6.0
      rimy=6.0

      XSIZE=xdim/1000.0+rimx
      YSIZE=ydim/1000.0+rimy

      x0=500.0*rimx
      y0=500.0*rimy


      fo='I_zpy-image'$
      +'_t'+strtrim(string(tt),2)$
      +'.eps'
      print,fo


      SET_PLOT,'PS'
      DEVICE,/color,/encapsulated,file=fo,XSIZE=XSIZE,YSIZE=YSIZE
      erase, 255      

      TVLCT, [0,255,0,0], [0,0,255,0], [0,0,0,255]
      plot, zpl(0:mm-1), pypl(0:mm-1), xs=1, ys=1, charsize=1.6, psym=3, $
      POSITION=[x0,y0,x0+xdim,y0+ydim],/DEVICE, $ 
      xtitle='!8z/!7l!8m',$
      xticks=2,$
      xr=[0.0,70.0e-6], $
      ytitle='!8p!Dy!N/!8mc',$
      yticks=2,$
      yr=[-0.05,0.05], $
      color=0,background=255 

      number=q_mass*mm
      ttt='!8t='+strtrim(string(tt*1.0e15),2)+'fs'
      number='!8N='+strtrim(string(number),2)
      xyouts,300,y0+ydim+2500.0,ttt,align=0.0,charsize=1.0,color=0,/DEVICE
      xyouts,300,y0+ydim+2200.0,number,align=0.0,charsize=1.0,color=0,/DEVICE

      DEVICE, /CLOSE
      SET_PLOT,'X'


; pxpy/py lineouts of phase space


      xdim=8000.0
      ydim=8000.0

      rimx=6.0
      rimy=6.0

      XSIZE=xdim/1000.0+rimx
      YSIZE=ydim/1000.0+rimy

      x0=500.0*rimx
      y0=500.0*rimy


      fo='I_pxpy-pz-image_lineout'$
      +'_t'+strtrim(string(tt),2)$
      +'.eps'
      print,fo


      SET_PLOT,'PS'
      DEVICE,/color,/encapsulated,file=fo,XSIZE=XSIZE,YSIZE=YSIZE
      erase, 255      

      sfxy=nf
      rmin=min(sfxy)
      rmax=max(sfxy)
      if (rmin eq 0.0) then begin
         rmin=-0.01
      endif
      if (rmax eq 0.0) then begin
         rmax=+0.01
      endif

      TVLCT, [0,255,0,0], [0,0,255,0], [0,0,0,255]
      plot, yy(0:ny), sfxy(100,0:ny), xs=1, ys=1, charsize=1.6, $
      POSITION=[x0,y0,x0+xdim,y0+ydim],/DEVICE, $ 
      xtitle='!8p!Dy!N/!8p!Dz',$
      xticks=2,$
      xr=[ys,yu], $
      ytitle='!8a.u',$
      yticks=2,$
      yr=[rmin,rmax], $
      color=0,background=255 

      number=q_mass*mm
      ttt='!8t='+strtrim(string(tt*1.0e15),2)+'fs'
      number='!8N='+strtrim(string(number),2)
      xyouts,300,y0+ydim+2500.0,ttt,align=0.0,charsize=1.0,color=0,/DEVICE
      xyouts,300,y0+ydim+2200.0,number,align=0.0,charsize=1.0,color=0,/DEVICE

      DEVICE, /CLOSE
      SET_PLOT,'X'

      endfor

end
