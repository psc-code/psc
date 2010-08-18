#include "config.h"

      function p_pulse_z2(x,y,z,t)

      implicit none
      real(kind=8) :: t,x,y,z
      real(kind=8) :: p_pulse_z2
      real(kind=8) :: p_pulse_z2_
      real(kind=8) :: C_p_pulse_z2

#ifdef USE_C_INTERFACE
      p_pulse_z2 = C_p_pulse_z2(x,y,z,t)
#else
      p_pulse_z2 = p_pulse_z2_(x,y,z,t)
#endif

      end function p_pulse_z2

c Laser pulse initialization (p-polarization).

      function p_pulse_z2_(x,y,z,t)

      use VLA_variables
      use PIC_variables

      implicit none
      real(kind=8) :: t,x,y,z
      real(kind=8) :: xm,ym,zm,zb
      real(kind=8) :: dxm,dym,dzm
      real(kind=8) :: xl,yl,zl
      real(kind=8) :: xr,yr,zr
      real(kind=8) :: p_pulse_z2_


c NOTE: The pulse is placed behind of the
c simulation box at a distance "zm" from the
c origin. The pulse then propagates into the 
c simulation box from the right. 


c  COORDINATE SYSTEM

c                                    ^ y       zm
c                                    |<----------------->
c                                    |
c                                    |             laser pulse
c            simulation box          |
c                                    |                | | |
c                                    |  ^   <-----    | | |
c                                    |  | ym          | | |
c                                    |  |
c          ------------------------------------------------->
c                              (i1n,i2n,i3n)=box origin    z 



c****************************************************************
c SETUP OF SHORT LASER PULSE
c****************************************************************

c dxm: width of pulse in x in m
c dym: width of pulse in y in m
c dzm: width of pulse in z in m
c xm: x-location of pulse center in m
c ym: y-location of pulse center in m
c zm: z-location of pulse center in m


      dxm=1.0*1.0e-6                       !modify for pulse setup
      dym=1.0*1.0e-6
      dzm=0.5*1.0e-6
      zb=1.0*1.0e-6
      xm=5.0*1.0e-6
      ym=5.0*1.0e-6
      zm=12.0*1.0e-6


      xm=xm/ld                              !normalization
      ym=ym/ld
      zm=zm/ld
      zb=zb/ld
      dxm=dxm/ld
      dym=dym/ld
      dzm=dzm/ld


      xl=x
      yl=y
      zl=z+t


      xr=xl-xm
      yr=yl-ym
      zr=zl-zm


      if (zl.lt.zm) then
         p_pulse_z2_=dsin(-zr)
     &           *dexp(-(xr/dxm)**2)
     &           *dexp(-(yr/dym)**2)
     &           *1.0/(1.0+dexp((abs(zr)-zb)/dzm))
      else if (zm.le.zl) then
          p_pulse_z2_=dsin(-zr)
     &           *dexp(-(xr/dxm)**2)
     &           *dexp(-(yr/dym)**2)
     &           *1.0/(1.0+dexp((abs(zr)-zb)/dzm))
      endif

      pos_z2 = thick*dz+1.5*zb+zm-lengthz/ld      ! ab

c Turn pulse off


      p_pulse_z2_=0.0d0


      end function p_pulse_z2_
