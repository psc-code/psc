;nmfc=' '
;print,'Data file name?'
;read,nmfc
nmfc=fi(ifi)


openr,1,nmfc,/COMPRESS
readf,1,nqcount
readf,1,i1n,i1x,i2n,i2x,i3n,i3x
readf,1,tm,np
readf,1,nmax,xnpe,ynpe
readf,1,nprparti,dnprparti,nistep
readf,1,dt,dx,dy,dz
readf,1,i0,e0,b0,j0,rho0,n0,cori,shift_z
readf,1,wl,wp,ld
readf,1,alpha,beta,eta
readf,1,vos,vt,c,qq,mm,tt,eps0


ne1=i1x-i1n+1
ne2=i2x-i2n+1
ne3=i3x-i3n+1


t=tm*dt/wl
x=fltarr(ne1+1)
y=fltarr(ne2+1)
z=fltarr(ne3+1)
for i=0L,ne1-1 do begin
    x(i)=(i1n+i)*dx*ld
endfor
for j=0L,ne2-1 do begin
    y(j)=(i2n+j)*dy*ld
endfor
for k=0L,ne3-1 do begin
    z(k)=(i3n+k)*dz*ld+shift_z*ld
endfor


datin=fltarr(11)
if (nqcount gt 0) then begin

   pnum=0L
   for i=0L,nqcount-step,step do begin
      pnum=pnum+1
   endfor

   xp=fltarr(pnum)
   yp=fltarr(pnum)
   zp=fltarr(pnum)
   pxp=fltarr(pnum)
   pyp=fltarr(pnum)
   pzp=fltarr(pnum)
   qnp=fltarr(pnum)
   mnp=fltarr(pnum)
   cnp=fltarr(pnum)
   lnp=fltarr(pnum)
   wnp=fltarr(pnum)

   pnum=0L
   for j=0L,nqcount-step,step do begin
      for i=1,step do begin
         readf,1,datin
      endfor
      xp(pnum)=datin(0)*ld
      yp(pnum)=datin(1)*ld
      zp(pnum)=datin(2)*ld+shift_z*ld
      pxp(pnum)=datin(3)
      pyp(pnum)=datin(4)
      pzp(pnum)=datin(5)
      qnp(pnum)=datin(6)
      mnp(pnum)=datin(7)
      cnp(pnum)=datin(8)
      lnp(pnum)=datin(9)
      wnp(pnum)=datin(10)
      pnum=pnum+1
   endfor
   nqcount=pnum
   print,nqcount
endif
close,1

;end
