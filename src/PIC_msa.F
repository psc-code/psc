c THIS SUBROUTINE DETERMINES A 3D SOLUTION OF MAXWELLS EQUATIONS.
c THE UPDATED E AND B-FIELDS ARE AT t=n*dt.


      subroutine PIC_msa

      use PIC_variables
      use VLA_variables
      use fields

      implicit none

      real(kind=8) :: cnx,cny,cnz
      real(kind=8) :: lx,ly,lz
      real(kind=8) :: jx,jy,jz

      real(kind=8) :: exh,eyh,ezh,bxh,byh,bzh
      real(kind=8) :: exk,eyk,ezk,bxk,byk,bzk

      real(kind=8) :: s_cpub
      real(kind=8) :: s_cpuh


      s_cpub=0.0d0
      s_cpuh=0.0d0
      call SERV_systime(cpug)


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
     &             +cny*(bz(i1,i2,i3)-bz(i1,i2-1,i3))
     &             -cnz*(by(i1,i2,i3)-by(i1,i2,i3-1))
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
     &                       +cny*(bz(i1,i2,i3)-bz(i1,i2-1,i3))
     &                       -cnz*(by(i1,i2,i3)-by(i1,i2,i3-1))
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
     &             +cnz*(bx(i1,i2,i3)-bx(i1,i2,i3-1))
     &             -cnx*(bz(i1,i2,i3)-bz(i1-1,i2,i3))
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
     &                       +cnz*(bx(i1,i2,i3)-bx(i1,i2,i3-1))
     &                       -cnx*(bz(i1,i2,i3)-bz(i1-1,i2,i3))
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
     &             +cnx*(by(i1,i2,i3)-by(i1-1,i2,i3))
     &             -cny*(bx(i1,i2,i3)-bx(i1,i2-1,i3))
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
     &                       +cnx*(by(i1,i2,i3)-by(i1-1,i2,i3))
     &                       -cny*(bx(i1,i2,i3)-bx(i1,i2-1,i3))
     &                       -0.5*dt*jz
            enddo
         enddo
      enddo


      call SERV_systime(cpua)
      call PIC_fex(ex)
      call PIC_fey(ex)
      call PIC_fez(ex)
      call PIC_fex(ey)
      call PIC_fey(ey)
      call PIC_fez(ey)
      call PIC_fex(ez)
      call PIC_fey(ez)
      call PIC_fez(ez)
      call SERV_systime(cpub)
      s_cpub=s_cpub+cpub-cpua


c B-field propagation E^(n+0.5), B^(n), j^(n), m^(n+0.5)
c -> E^(n+0.5), B^(n+0.5), j^(n), m^(n+0.5)


      bx2A=0.0
      do i3=i3mn,i3mx
         do i2=i2mn,i2mx
            do i1=i1mn,i1mx
               bxh=bx(i1,i2,i3)
               bxk=bx(i1,i2,i3)
     &             -cny*(ez(i1,i2+1,i3)-ez(i1,i2,i3))
     &             +cnz*(ey(i1,i2,i3+1)-ey(i1,i2,i3))
               bx2A=bx2A+dx*dy*dz*bxh*bxk
            enddo
         enddo
      enddo
      do i3=i3mn-1,i3mx+1
         do i2=i2mn-1,i2mx+1
            do i1=i1mn-1,i1mx+1
               bx(i1,i2,i3)=bx(i1,i2,i3)
     &                      -cny*(ez(i1,i2+1,i3)-ez(i1,i2,i3))
     &                      +cnz*(ey(i1,i2,i3+1)-ey(i1,i2,i3))
            enddo
         enddo
      enddo

      by2A=0.0
      do i3=i3mn,i3mx
         do i2=i2mn,i2mx
            do i1=i1mn,i1mx
               byh=by(i1,i2,i3)
               byk=by(i1,i2,i3)
     &             -cnz*(ex(i1,i2,i3+1)-ex(i1,i2,i3))
     &             +cnx*(ez(i1+1,i2,i3)-ez(i1,i2,i3))
               by2A=by2A+dx*dy*dz*byh*byk
            enddo
         enddo
      enddo
      do i3=i3mn-1,i3mx+1
         do i2=i2mn-1,i2mx+1
            do i1=i1mn-1,i1mx+1
               by(i1,i2,i3)=by(i1,i2,i3)
     &                      -cnz*(ex(i1,i2,i3+1)-ex(i1,i2,i3))
     &                      +cnx*(ez(i1+1,i2,i3)-ez(i1,i2,i3))
            enddo
         enddo
      enddo

      bz2A=0.0
      do i3=i3mn,i3mx
         do i2=i2mn,i2mx
            do i1=i1mn,i1mx
               bzh=bz(i1,i2,i3)
               bzk=bz(i1,i2,i3)
     &             -cnx*(ey(i1+1,i2,i3)-ey(i1,i2,i3))
     &             +cny*(ex(i1,i2+1,i3)-ex(i1,i2,i3))
               bz2A=bz2A+dx*dy*dz*bzh*bzk
            enddo
         enddo
      enddo
      do i3=i3mn-1,i3mx+1
         do i2=i2mn-1,i2mx+1
            do i1=i1mn-1,i1mx+1
               bz(i1,i2,i3)=bz(i1,i2,i3)
     &                      -cnx*(ey(i1+1,i2,i3)-ey(i1,i2,i3))
     &                      +cny*(ex(i1,i2+1,i3)-ex(i1,i2,i3))
            enddo
         enddo
      enddo


      call SERV_systime(cpua)
      call PIC_fex(bx)
      call PIC_fey(bx)
      call PIC_fez(bx)
      call PIC_fex(by)
      call PIC_fey(by)
      call PIC_fez(by)
      call PIC_fex(bz)
      call PIC_fey(bz)
      call PIC_fez(bz)
      call SERV_systime(cpub)
      s_cpub=s_cpub+cpub-cpua


c energy conservation


      fluxi=0.0
      fluxo=0.0
      do i2=i2mn,i2mx
         do i1=i1mn,i1mx
            fluxi=fluxi+0.25*dx*dy
     &            *((ex(i1,i2,i3mn+1)+by(i1,i2,i3mn))**2
     &             +(ey(i1,i2,i3mn)-bx(i1,i2,i3mn-1))**2)
            fluxo=fluxo+0.25*dx*dy
     &            *((ex(i1,i2,i3mn+1)-by(i1,i2,i3mn))**2
     &             +(ey(i1,i2,i3mn)+bx(i1,i2,i3mn-1))**2)
         enddo
      enddo


      pox=0.0
      poy=0.0
      poz=0.0
      do i3=i3mn,i3mx
         do i2=i2mn,i2mx
            do i1=i1mn,i1mx
               pox=pox+dy*dz
     &             *(ey(i1+1,i2,i3)*bz(i1,i2,i3)
     &              -ey(i1,i2,i3)*bz(i1-1,i2,i3)
     &              -ez(i1+1,i2,i3)*by(i1,i2,i3)
     &              +ez(i1,i2,i3)*by(i1-1,i2,i3))
               poy=poy+dx*dz
     &             *(ez(i1,i2+1,i3)*bx(i1,i2,i3)
     &              -ez(i1,i2,i3)*bx(i1,i2-1,i3)
     &              -ex(i1,i2+1,i3)*bz(i1,i2,i3)
     &              +ex(i1,i2,i3)*bz(i1,i2-1,i3))
               poz=poz+dx*dy
     &             *(ex(i1,i2,i3+1)*by(i1,i2,i3)
     &              -ex(i1,i2,i3)*by(i1,i2,i3-1)
     &              -ey(i1,i2,i3+1)*bx(i1,i2,i3)
     &              +ey(i1,i2,i3)*bx(i1,i2,i3-1))
            enddo
         enddo
      enddo


      call SERV_systime(cpuh)
      s_cpuh=s_cpuh+cpuh-cpug


      cpumess=cpumess+s_cpub
      cpucomp=cpucomp+s_cpuh-s_cpub


      end subroutine PIC_msa
