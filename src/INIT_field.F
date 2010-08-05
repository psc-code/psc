
#include "config.h"

c     ----------------------------------------------------------------------
      subroutine ALLOC_field_fortran
c     ----------------------------------------------------------------------

      use PIC_variables
      use VLA_variables

      implicit none

      call OUT_alloc_fields

      allocate(cpu(1:5))
      allocate(cpu_ary(1:5,0:npe-1))

! array allocation
! size of pml should be equal to size of em field arrays

! x

      allocate(kappax(i1mn-rd1:i1mx+rd1))
      allocate(sigmax(i1mn-rd1:i1mx+rd1))
      
      allocate(cxp(i1mn-rd1:i1mx+rd1))
      allocate(cxm(i1mn-rd1:i1mx+rd1))
     
      allocate(fbx(i1mn-rd1:i1mx+rd1))
      allocate(fcx(i1mn-rd1:i1mx+rd1))
      allocate(fdx(i1mn-rd1:i1mx+rd1))
      allocate(fex(i1mn-rd1:i1mx+rd1))

      allocate(bxp(i1mn-rd1:i1mx+rd1))
      allocate(bxm(i1mn-rd1:i1mx+rd1))
      
      allocate(gbx(i1mn-rd1:i1mx+rd1))
      allocate(gcx(i1mn-rd1:i1mx+rd1))
      allocate(gdx(i1mn-rd1:i1mx+rd1))
      allocate(gex(i1mn-rd1:i1mx+rd1))

! y
      
      allocate(kappay(i2mn-rd2:i2mx+rd2))
      allocate(sigmay(i2mn-rd2:i2mx+rd2))
      
      allocate(cyp(i2mn-rd2:i2mx+rd2))
      allocate(cym(i2mn-rd2:i2mx+rd2))
      
      allocate(fby(i2mn-rd2:i2mx+rd2))
      allocate(fcy(i2mn-rd2:i2mx+rd2))
      allocate(fdy(i2mn-rd2:i2mx+rd2))
      allocate(fey(i2mn-rd2:i2mx+rd2))

      allocate(byp(i2mn-rd2:i2mx+rd2))
      allocate(bym(i2mn-rd2:i2mx+rd2))
      
      allocate(gby(i2mn-rd2:i2mx+rd2))
      allocate(gcy(i2mn-rd2:i2mx+rd2))
      allocate(gdy(i2mn-rd2:i2mx+rd2))
      allocate(gey(i2mn-rd2:i2mx+rd2))

! z
      
      allocate(kappaz(i3mn-rd3:i3mx+rd3))
      allocate(sigmaz(i3mn-rd3:i3mx+rd3))
      
      allocate(czp(i3mn-rd3:i3mx+rd3))
      allocate(czm(i3mn-rd3:i3mx+rd3))
      
      allocate(fbz(i3mn-rd3:i3mx+rd3))
      allocate(fcz(i3mn-rd3:i3mx+rd3))
      allocate(fdz(i3mn-rd3:i3mx+rd3))
      allocate(fez(i3mn-rd3:i3mx+rd3))

      allocate(bzp(i3mn-rd3:i3mx+rd3))
      allocate(bzm(i3mn-rd3:i3mx+rd3))
      
      allocate(gbz(i3mn-rd3:i3mx+rd3))
      allocate(gcz(i3mn-rd3:i3mx+rd3))
      allocate(gdz(i3mn-rd3:i3mx+rd3))
      allocate(gez(i3mn-rd3:i3mx+rd3))

#ifdef USE_PML
      call PML_alloc_field
#endif

      end subroutine ALLOC_field_fortran

c THIS SUBROUTINE INITIALIZES THE ELECTROMAGNETIC FIELDS.


      subroutine INIT_field

      use PIC_variables
      use VLA_variables
      use fields

      implicit none
      integer :: l1,l2,l3

      real(kind=8) :: h1,h2,h3
      real(kind=8) :: p_pulse_z2

c initialization of the fields

      do i3=i3mn-rd3,i3mx+rd3
         h3=i3*dz
         do i2=i2mn-rd2,i2mx+rd2
            h2=i2*dy
            do i1=i1mn-rd1,i1mx+rd1
               h1=i1*dx
               ex(i1,i2,i3)=0.0d0
               ey(i1,i2,i3)=p_pulse_z2(h1,h2+0.5*dy,h3,-0.5*dt)
               ey(i1,i2,i3)=0.0d0
               ez(i1,i2,i3)=0.0d0
               bx(i1,i2,i3)=p_pulse_z2(h1,h2+0.5*dy,h3+0.5*dz,0.0*dt)
               bx(i1,i2,i3)=0.0d0
               by(i1,i2,i3)=0.0d0
               bz(i1,i2,i3)=0.0d0

               jxi(i1,i2,i3)=0.0d0
               jyi(i1,i2,i3)=0.0d0
               jzi(i1,i2,i3)=0.0d0

               ne(i1,i2,i3)=0.0d0
               ni(i1,i2,i3)=0.0d0
               nn(i1,i2,i3)=0.0d0
            enddo
         enddo
      enddo

      end subroutine INIT_field

C     set up various things after the main fields have been initialized

      subroutine SETUP_field

      use PIC_variables
      use VLA_variables

      implicit none
      integer :: l1,l2,l3

      real(kind=8) :: xi,yi,zi
      real(kind=8) :: x0,y0,z0,xr,yr,zr,rot

      x0=5.0*1.0e-6/ld
      y0=5.0*1.0e-6/ld
      z0=5.0*1.0e-6/ld
      rot=45.0*3.141592/180.0

c setting up charge and current densities

c FIXME
#if 0
      if (niloc.gt.0) then
         do l=1,niloc

            xi=p_niloc(11*l)
            yi=p_niloc(11*l+1)
            zi=p_niloc(11*l+2)

            xr=xi
            yr=cos(rot)*(yi-y0)-sin(rot)*(zi-z0)+y0
            zr=cos(rot)*(zi-z0)+sin(rot)*(yi-y0)+z0

            if (((r3n-1)*dz.le.zr.and.zr.le.(r3x+1)*dz).and.
     &          ((r2n-1)*dy.le.yr.and.yr.le.(r2x+1)*dy).and.
     &          ((r1n-1)*dx.le.xr.and.xr.le.(r1x+1)*dx)) then

               p_niloc(11*l+9)=0.0

            endif

         enddo
      endif
#endif



c setting up time-averaged fields

      call OUT_setup_fields

#ifdef USE_PML
      call PML_setup_field
#endif


      cpu=0.0d0
      cpu_ary=0.0d0

      end subroutine SETUP_field

#ifdef USE_PML
c     ----------------------------------------------------------------------
      subroutine PML_alloc_field
c     ----------------------------------------------------------------------

      use PIC_variables
      use VLA_variables

      implicit none

      ! additional d field  - added by ab

      allocate(dvx(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))
      allocate(dvy(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))
      allocate(dvz(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))


      ! additional h field - added by ab

      allocate(hx(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))
      allocate(hy(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))
      allocate(hz(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))

      ! electric and magnetic permittivity tensors - added by ab

      allocate(eps(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))
      allocate(mu(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))

      call OUT_alloc_fields_pml

      dvx=0.0d0
      dvy=0.0d0
      dvz=0.0d0
      hx=0.0d0
      hy=0.0d0
      hz=0.0d0
      eps = 1.0
      mu = 1.0

      end subroutine PML_alloc_field

c     ----------------------------------------------------------------------
      subroutine PML_setup_field
c     ----------------------------------------------------------------------

      use PIC_variables
      use VLA_variables
      use fields

      implicit none

      ! initialisation of additional fields - added by ab
      ! because of normalization these are equal to e and b fields

      dvx=ex
      dvy=ey
      dvz=ez
      hx=bx
      hy=by
      hz=bz
      
      ! initialisation of permittivity tensors, always set to 1.0
      ! due to scaling these are dimensionless
      
      eps = 1.0
      mu = 1.0

      call OUT_setup_fields_pml

      fluxit=0.0d0
      fluxot=0.0d0

      enEXt=0.0d0
      enEYt=0.0d0
      enEZt=0.0d0
      enBXt=0.0d0
      enBYt=0.0d0
      enBZt=0.0d0
      enHXt=0.0d0   ! ab
      enHYt=0.0d0   ! ab
      enHZt=0.0d0   ! ab
      enDXt=0.0d0   ! ab
      enDYt=0.0d0   ! ab
      enDZt=0.0d0   ! ab

      ent=0.0d0
      poxt=0.0d0
      poyt=0.0d0
      pozt=0.0d0
      jet=0.0d0

      sum_fed=0.0d0
      sum_ped=0.0d0

! deriving pml attenuation coefficients

! not moved by half step
      
      ! x

      do i1 = i1mn-rd1, i1mx+rd1
         kappax(i1) = 1.0
         sigmax(i1) = 0.0
         cxp(i1) = 2*eps0*kappax(i1)+sigmax(i1)*dt
         cxm(i1) = 2*eps0*kappax(i1)-sigmax(i1)*dt
         fbx(i1) = 2*eps0*kappax(i1)
         fcx(i1) = cxm(i1)/cxp(i1)
         fdx(i1) = 2*eps0*dt/cxp(i1)
         fex(i1) = 1.0/cxp(i1)
      end do

      ! y
      
       do i2 = i2mn-rd2, i2mx+rd2
         kappay(i2) = 1.0
         sigmay(i2) = 0.0
         cyp(i2) = 2*eps0*kappay(i2)+sigmay(i2)*dt
         cym(i2) = 2*eps0*kappay(i2)-sigmay(i2)*dt
         fby(i2) = 2*eps0*kappay(i2)
         fcy(i2) = cym(i2)/cyp(i2)
         fdy(i2) = 2*eps0*dt/cyp(i2)
         fey(i2) = 1.0/cyp(i2)
      end do

      ! z

      do i3 = i3mn-rd3, i3mx+rd3
         kappaz(i3) = 1.0
         sigmaz(i3) = 0.0
         czp(i3) = 2*eps0*kappaz(i3)+sigmaz(i3)*dt
         czm(i3) = 2*eps0*kappaz(i3)-sigmaz(i3)*dt
         fbz(i3) = 2*eps0*kappaz(i3)
         fcz(i3) = czm(i3)/czp(i3)
         fdz(i3) = 2*eps0*dt/czp(i3)
         fez(i3) = 1.0/czp(i3)
      end do


! moved by half step

      ! x

      do i1 = i1mn-rd1, i1mx+rd1
         kappax(i1) = 1.0
         sigmax(i1) = 0.0
         bxp(i1) = 2*eps0*kappax(i1)+sigmax(i1)*dt
         bxm(i1) = 2*eps0*kappax(i1)-sigmax(i1)*dt
         gbx(i1) = 2*eps0*kappax(i1)
         gcx(i1) = bxm(i1)/bxp(i1)
         gdx(i1) = 2*eps0*dt/bxp(i1)
         gex(i1) = 1.0/bxp(i1)
      end do

      ! y

      do i2 = i2mn-rd2, i2mx+rd2
         kappay(i2) = 1.0
         sigmay(i2) = 0.0
         byp(i2) = 2*eps0*kappay(i2)+sigmay(i2)*dt
         bym(i2) = 2*eps0*kappay(i2)-sigmay(i2)*dt
         gby(i2) = 2*eps0*kappay(i2)
         gcy(i2) = bym(i2)/byp(i2)
         gdy(i2) = 2*eps0*dt/byp(i2)
         gey(i2) = 1.0/byp(i2)
      end do

      ! z

      do i3 = i3mn-rd3, i3mx+rd3
         kappaz(i3) = 1.0
         sigmaz(i3) = 0.0
         bzp(i3) = 2*eps0*kappaz(i3)+sigmaz(i3)*dt
         bzm(i3) = 2*eps0*kappaz(i3)-sigmaz(i3)*dt
         gbz(i3) = 2*eps0*kappaz(i3)
         gcz(i3) = bzm(i3)/bzp(i3)
         gdz(i3) = 2*eps0*dt/bzp(i3)
         gez(i3) = 1.0/bzp(i3)
      end do

      end subroutine PML_setup_field
#endif

