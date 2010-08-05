
#include "config.h"

c     ----------------------------------------------------------------------
      subroutine ALLOC_field_fortran
c     ----------------------------------------------------------------------

      use PIC_variables
      use VLA_variables

      implicit none

      allocate(cpu(1:5))
      allocate(cpu_ary(1:5,0:npe-1))

#ifdef USE_PML
      call PML_coeff_alloc
      call PML_fields_alloc
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
               hx(i1,i2,i3)=p_pulse_z2(h1,h2+0.5*dy,h3+0.5*dz,0.0*dt)
               hx(i1,i2,i3)=0.0d0
               hy(i1,i2,i3)=0.0d0
               hz(i1,i2,i3)=0.0d0

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

#ifdef USE_PML
      call PML_fields_setup
      call PML_coeff_setup
#endif


      cpu=0.0d0
      cpu_ary=0.0d0

      end subroutine SETUP_field

