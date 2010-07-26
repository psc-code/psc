
#include "config.h"

      function p_pulse_z1(x,y,z,t)

      implicit none
      real(kind=8) :: t,x,y,z
      real(kind=8) :: p_pulse_z1
      real(kind=8) :: p_pulse_z1_
      real(kind=8) :: C_p_pulse_z1

#ifdef USE_C_INTERFACE
      p_pulse_z1 = C_p_pulse_z1(x,y,z,t)
#else
      p_pulse_z1 = p_pulse_z1_(x,y,z,t)
#endif

      end function p_pulse_z1


c Laser pulse initialization (p-polarization).

      function p_pulse_z1_(x,y,z,t)

      use VLA_variables
      use PIC_variables

      implicit none
      real(kind=8) :: t,x,y,z
      real(kind=8) :: xm,ym,zm,zb
      real(kind=8) :: dxm,dym,dzm
      real(kind=8) :: xl,yl,zl
      real(kind=8) :: xr,yr,zr
      real(kind=8) :: p_pulse_z1_


c NOTE: The pulse is placed behind of the
c simulation box at a distance "zm" from the
c origin. The pulse then propagates into the 
c simulation box from the left. 


c  COORDINATE SYSTEM

c                          zm        ^ y
c                 <----------------->|
c                                    |
c            laser pulse             |
c                                    |     simulation
c               | | |                |     box
c               | | |   ----->   ^   |
c               | | |         ym |   |
c                                |   |
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


      dxm=4.0*1.0e-6                       !modify for pulse setup
      dym=4.0*1.0e-6
      dzm=0.1*1.0e-6
      zb=300.0*1.0e-6
      xm=0.01*1.0e-6
      ym=0.01*1.0e-6
      zm=-301.0*1.0e-6


      xm=xm/ld                              !normalization
      ym=ym/ld
      zm=zm/ld
      zb=zb/ld
      dxm=dxm/ld
      dym=dym/ld
      dzm=dzm/ld


      xl=x
      yl=y
      zl=z-t


      xr=xl-xm
      yr=yl-ym
      zr=zl-zm


      if (zl.gt.zm) then
         p_pulse_z1_=dsin(zr)
     &           *dexp(-(xr/dxm)**2)
     &           *dexp(-(yr/dym)**2)
     &           *1.0/(1.0+dexp((abs(zr)-zb)/dzm))
      else if (zm.ge.zl) then
         p_pulse_z1_=dsin(zr)
     &           *dexp(-(xr/dxm)**2)
     &           *dexp(-(yr/dym)**2)
     &           *1.0/(1.0+dexp((abs(zr)-zb)/dzm))
      endif

      pos_z1 = thick*dz+1.5*zb-zm        ! ab

c Turn pulse off


      p_pulse_z1_=0.0d0


      end function p_pulse_z1_
