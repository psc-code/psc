c THIS SUBROUTINE DETERMINES A 3D SOLUTION OF MAXWELLS EQUATIONS.
c THE UPDATED E AND H-FIELDS ARE AT t=n*dt.


      subroutine PIC_msa_e(ex, ey, ez, hx, hy, hz, jxi, jyi, jzi)

      use VLA_variables
      use PSC_patch
      use OUT_poyn

      implicit none

      real(kind=8),dimension(i1mn-rd1:i1mx+rd1,
     &                       i2mn-rd2:i2mx+rd2,
     &                       i3mn-rd3:i3mx+rd3) :: ex, ey, ez
      real(kind=8),dimension(i1mn-rd1:i1mx+rd1,
     &                       i2mn-rd2:i2mx+rd2,
     &                       i3mn-rd3:i3mx+rd3) :: hx, hy, hz
      real(kind=8),dimension(i1mn-rd1:i1mx+rd1,
     &                       i2mn-rd2:i2mx+rd2,
     &                       i3mn-rd3:i3mx+rd3) :: jxi, jyi, jzi

      real(kind=8) :: cnx,cny,cnz
      real(kind=8) :: lx,ly,lz
      real(kind=8) :: jx,jy,jz

      real(kind=8) :: exh,eyh,ezh,hxh,hyh,hzh
      real(kind=8) :: exk,eyk,ezk,hxk,hyk,hzk

      integer :: i1, i2, i3

c initialization


      lx=dt/dx
      ly=dt/dy
      lz=dt/dz

      cnx=0.5*lx
      cny=0.5*ly
      cnz=0.5*lz


c energy conservation


      je=0.0
      do i3=i3mn,i3mx
         do i2=i2mn,i2mx
            do i1=i1mn,i1mx
               jx=jxi(i1,i2,i3)
               jy=jyi(i1,i2,i3)
               jz=jzi(i1,i2,i3)
               je=je+0.5*dx*dy*dz*(jx*ex(i1,i2,i3)
     &                            +jy*ey(i1,i2,i3)
     &                            +jz*ez(i1,i2,i3))
            enddo
         enddo
      enddo


c E-field propagation E^(n), B^(n), j^(n) 
c -> E^(n+0.5), B^(n), j^(n)


      ex2A=0.0
      do i3=i3mn,i3mx
         do i2=i2mn,i2mx
            do i1=i1mn,i1mx
               jx=jxi(i1,i2,i3)
               exh=ex(i1,i2,i3)
               exk=ex(i1,i2,i3)
     &             +cny*(hz(i1,i2,i3)-hz(i1,i2-1,i3))
     &             -cnz*(hy(i1,i2,i3)-hy(i1,i2,i3-1))
     &             -0.5*dt*jx
               ex2A=ex2A+dx*dy*dz*exh*exk
            enddo
         enddo
      enddo
      do i3=i3mn-1,i3mx+1
         do i2=i2mn-1,i2mx+1
            do i1=i1mn-1,i1mx+1
               jx=jxi(i1,i2,i3)
               ex(i1,i2,i3)=ex(i1,i2,i3)
     &                       +cny*(hz(i1,i2,i3)-hz(i1,i2-1,i3))
     &                       -cnz*(hy(i1,i2,i3)-hy(i1,i2,i3-1))
     &                       -0.5*dt*jx
            enddo
         enddo
      enddo

      ey2A=0.0
      do i3=i3mn,i3mx
         do i2=i2mn,i2mx
            do i1=i1mn,i1mx
               jy=jyi(i1,i2,i3)
               eyh=ey(i1,i2,i3)
               eyk=ey(i1,i2,i3)
     &             +cnz*(hx(i1,i2,i3)-hx(i1,i2,i3-1))
     &             -cnx*(hz(i1,i2,i3)-hz(i1-1,i2,i3))
     &             -0.5*dt*jy
               ey2A=ey2A+dx*dy*dz*eyh*eyk
            enddo
         enddo
      enddo
      do i3=i3mn-1,i3mx+1
         do i2=i2mn-1,i2mx+1
            do i1=i1mn-1,i1mx+1
               jy=jyi(i1,i2,i3)
               ey(i1,i2,i3)=ey(i1,i2,i3)
     &                       +cnz*(hx(i1,i2,i3)-hx(i1,i2,i3-1))
     &                       -cnx*(hz(i1,i2,i3)-hz(i1-1,i2,i3))
     &                       -0.5*dt*jy
            enddo
         enddo
      enddo

      ez2A=0.0
      do i3=i3mn,i3mx
         do i2=i2mn,i2mx
            do i1=i1mn,i1mx
               jz=jzi(i1,i2,i3)
               ezh=ez(i1,i2,i3)
               ezk=ez(i1,i2,i3)
     &             +cnx*(hy(i1,i2,i3)-hy(i1-1,i2,i3))
     &             -cny*(hx(i1,i2,i3)-hx(i1,i2-1,i3))
     &             -0.5*dt*jz
               ez2A=ez2A+dx*dy*dz*ezh*ezk
            enddo
         enddo
      enddo
      do i3=i3mn-1,i3mx+1
         do i2=i2mn-1,i2mx+1
            do i1=i1mn-1,i1mx+1
               jz=jzi(i1,i2,i3)
               ez(i1,i2,i3)=ez(i1,i2,i3)
     &                       +cnx*(hy(i1,i2,i3)-hy(i1-1,i2,i3))
     &                       -cny*(hx(i1,i2,i3)-hx(i1,i2-1,i3))
     &                       -0.5*dt*jz
            enddo
         enddo
      enddo

      end subroutine PIC_msa_e


      subroutine PIC_msa_h(ex, ey, ez, hx, hy, hz, jxi, jyi, jzi)

      use VLA_variables
      use PSC_patch
      use OUT_poyn

      implicit none

      real(kind=8),dimension(i1mn-rd1:i1mx+rd1,
     &                       i2mn-rd2:i2mx+rd2,
     &                       i3mn-rd3:i3mx+rd3) :: ex, ey, ez
      real(kind=8),dimension(i1mn-rd1:i1mx+rd1,
     &                       i2mn-rd2:i2mx+rd2,
     &                       i3mn-rd3:i3mx+rd3) :: hx, hy, hz
      real(kind=8),dimension(i1mn-rd1:i1mx+rd1,
     &                       i2mn-rd2:i2mx+rd2,
     &                       i3mn-rd3:i3mx+rd3) :: jxi, jyi, jzi

      real(kind=8) :: cnx,cny,cnz
      real(kind=8) :: lx,ly,lz
      real(kind=8) :: jx,jy,jz

      real(kind=8) :: exh,eyh,ezh,hxh,hyh,hzh
      real(kind=8) :: exk,eyk,ezk,hxk,hyk,hzk

      integer :: i1, i2, i3

c initialization


      lx=dt/dx
      ly=dt/dy
      lz=dt/dz

      cnx=0.5*lx
      cny=0.5*ly
      cnz=0.5*lz


c B-field propagation E^(n+0.5), B^(n), j^(n), m^(n+0.5)
c -> E^(n+0.5), B^(n+0.5), j^(n), m^(n+0.5)


      hx2A=0.0
      do i3=i3mn,i3mx
         do i2=i2mn,i2mx
            do i1=i1mn,i1mx
               hxh=hx(i1,i2,i3)
               hxk=hx(i1,i2,i3)
     &             -cny*(ez(i1,i2+1,i3)-ez(i1,i2,i3))
     &             +cnz*(ey(i1,i2,i3+1)-ey(i1,i2,i3))
               hx2A=hx2A+dx*dy*dz*hxh*hxk
            enddo
         enddo
      enddo
      do i3=i3mn-1,i3mx+1
         do i2=i2mn-1,i2mx+1
            do i1=i1mn-1,i1mx+1
               hx(i1,i2,i3)=hx(i1,i2,i3)
     &                      -cny*(ez(i1,i2+1,i3)-ez(i1,i2,i3))
     &                      +cnz*(ey(i1,i2,i3+1)-ey(i1,i2,i3))
            enddo
         enddo
      enddo

      hy2A=0.0
      do i3=i3mn,i3mx
         do i2=i2mn,i2mx
            do i1=i1mn,i1mx
               hyh=hy(i1,i2,i3)
               hyk=hy(i1,i2,i3)
     &             -cnz*(ex(i1,i2,i3+1)-ex(i1,i2,i3))
     &             +cnx*(ez(i1+1,i2,i3)-ez(i1,i2,i3))
               hy2A=hy2A+dx*dy*dz*hyh*hyk
            enddo
         enddo
      enddo
      do i3=i3mn-1,i3mx+1
         do i2=i2mn-1,i2mx+1
            do i1=i1mn-1,i1mx+1
               hy(i1,i2,i3)=hy(i1,i2,i3)
     &                      -cnz*(ex(i1,i2,i3+1)-ex(i1,i2,i3))
     &                      +cnx*(ez(i1+1,i2,i3)-ez(i1,i2,i3))
            enddo
         enddo
      enddo

      hz2A=0.0
      do i3=i3mn,i3mx
         do i2=i2mn,i2mx
            do i1=i1mn,i1mx
               hzh=hz(i1,i2,i3)
               hzk=hz(i1,i2,i3)
     &             -cnx*(ey(i1+1,i2,i3)-ey(i1,i2,i3))
     &             +cny*(ex(i1,i2+1,i3)-ex(i1,i2,i3))
               hz2A=hz2A+dx*dy*dz*hzh*hzk
            enddo
         enddo
      enddo
      do i3=i3mn-1,i3mx+1
         do i2=i2mn-1,i2mx+1
            do i1=i1mn-1,i1mx+1
               hz(i1,i2,i3)=hz(i1,i2,i3)
     &                      -cnx*(ey(i1+1,i2,i3)-ey(i1,i2,i3))
     &                      +cny*(ex(i1,i2+1,i3)-ex(i1,i2,i3))
            enddo
         enddo
      enddo


c energy conservation // FIXME, probably broken since ghost points haven't been set


      fluxi=0.0
      fluxo=0.0
      do i2=i2mn,i2mx
         do i1=i1mn,i1mx
            fluxi=fluxi+0.25*dx*dy
     &            *((ex(i1,i2,i3mn+1)+hy(i1,i2,i3mn))**2
     &             +(ey(i1,i2,i3mn)-hx(i1,i2,i3mn-1))**2)
            fluxo=fluxo+0.25*dx*dy
     &            *((ex(i1,i2,i3mn+1)-hy(i1,i2,i3mn))**2
     &             +(ey(i1,i2,i3mn)+hx(i1,i2,i3mn-1))**2)
         enddo
      enddo


      pox=0.0
      poy=0.0
      poz=0.0
      do i3=i3mn,i3mx
         do i2=i2mn,i2mx
            do i1=i1mn,i1mx
               pox=pox+dy*dz
     &             *(ey(i1+1,i2,i3)*hz(i1,i2,i3)
     &              -ey(i1,i2,i3)*hz(i1-1,i2,i3)
     &              -ez(i1+1,i2,i3)*hy(i1,i2,i3)
     &              +ez(i1,i2,i3)*hy(i1-1,i2,i3))
               poy=poy+dx*dz
     &             *(ez(i1,i2+1,i3)*hx(i1,i2,i3)
     &              -ez(i1,i2,i3)*hx(i1,i2-1,i3)
     &              -ex(i1,i2+1,i3)*hz(i1,i2,i3)
     &              +ex(i1,i2,i3)*hz(i1,i2-1,i3))
               poz=poz+dx*dy
     &             *(ex(i1,i2,i3+1)*hy(i1,i2,i3)
     &              -ex(i1,i2,i3)*hy(i1,i2,i3-1)
     &              -ey(i1,i2,i3+1)*hx(i1,i2,i3)
     &              +ey(i1,i2,i3)*hx(i1,i2,i3-1))
            enddo
         enddo
      enddo


      end subroutine PIC_msa_h
