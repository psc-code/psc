;nmfc=' '
;print,'Data file name?'
;read,nmfc
;nmfc=fi(ifi)


openr,1,nmfc,/COMPRESS
readf,1,i1n,i1x,i2n,i2x,i3n,i3x
readf,1,tm,np
readf,1,nmax,xnpe,ynpe
readf,1,nprf,dnprf
readf,1,dt,dx,dy,dz,shift_z
readf,1,i0,e0,b0,j0,rho0,n0
readf,1,wl,wp,ld
readf,1,alpha,beta,eta
readf,1,vos,vt,c,qq,mm,tt,eps0


ne1=i1x-i1n+1
ne2=i2x-i2n+1
ne3=i3x-i3n+1
serie=fltarr(ne1+1,ne2+1,ne3+1)


for k=0,ne3-1 do begin
   for j=0,ne2-1 do begin
      for i=0,ne1-1 do begin
         readf,1,we
         serie(i,j,k)=we
      endfor
   endfor
endfor


t=tm*dt/wl
x=fltarr(ne1+1)
y=fltarr(ne2+1)
z=fltarr(ne3+1)
for i=0,ne1-1 do begin
    x(i)=(i1n+i)*dx*ld
endfor
for j=0,ne2-1 do begin
    y(j)=(i2n+j)*dy*ld
endfor
for k=0,ne3-1 do begin
    z(k)=(i3n+k)*dz*ld+shift_z*ld
endfor
close,1

;end
