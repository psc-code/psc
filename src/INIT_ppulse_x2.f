c Laser pulse initialization (p-polarization).

      function p_pulse_x2(x,y,z,t)

      use VLA_variables
      use PIC_variables

      implicit none
      real(kind=8) :: t,x,y,z
      real(kind=8) :: xm,ym,zm
      real(kind=8) :: dxm,dym,dzm
      real(kind=8) :: xl,yl,zl
      real(kind=8) :: xr,yr,zr
      real(kind=8) :: p_pulse_x2


c NOTE: The pulse is placed in front of the
c simulation box at a distance "zm". The pulse 
c then propagates into the simulation box
c from the left edge. 


c  COORDINATE SYSTEM

c                                    ^ y       xm
c                                    |<----------------->
c                                    |
c                                    |             laser pulse
c            simulation box          |
c                                    |                | | |
c                                    |  ^   <-----    | | |
c                                    |  | ym          | | |
c                                    |  |
c          ------------------------------------------------->
c                              (i1n,i2n,i3n)=box origin    x 


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
      dym=5.0*1.0e-6
      dzm=5.0*1.0e-6
      xm=202.0*1.0e-6
      ym=10.0*1.0e-6
      zm=10.0*1.0e-6


      xm=xm/ld                              !normalization
      ym=ym/ld
      zm=zm/ld
      dxm=dxm/ld
      dym=dym/ld
      dzm=dzm/ld


      xl=x+t
      yl=y
      zl=z


      xr=xl-xm
      yr=yl-ym
      zr=zl-zm


      if (xl.lt.xm) then
         p_pulse_x2=dsin(-xr)
     &           *dexp(-(xr/dxm)**2)
     &           *dexp(-(yr/dym)**2)
     &           *dexp(-(zr/dzm)**2)
      else if (xm.le.xl) then
          p_pulse_x2=dsin(xr)
     &           *dexp(-(xr/dxm)**2)
     &           *dexp(-(yr/dym)**2)
     &           *dexp(-(zr/dzm)**2)
      endif

      pos_x2 = thick*dx+3.0*dxm+xm-lengthx/ld      ! ab

c Turn pulse off


      p_pulse_x2=0.0d0


      end function p_pulse_x2
