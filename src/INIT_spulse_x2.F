c Laser pulse initialization (p-polarization).

      function s_pulse_x2(x,y,z,t)

      implicit none
      real(kind=8) :: t,x,y,z
      real(kind=8) :: s_pulse_x2
      real(kind=8) :: C_s_pulse_x2

      s_pulse_x2 = C_s_pulse_x2(x,y,z,t)

      end function s_pulse_x2

      function s_pulse_x2_(x,y,z,t)

      use PSC_params

      implicit none
      real(kind=8) :: t,x,y,z
      real(kind=8) :: xm,ym,zm
      real(kind=8) :: dxm,dym,dzm
      real(kind=8) :: xl,yl,zl
      real(kind=8) :: xr,yr,zr
      real(kind=8) :: s_pulse_x2_


c NOTE: The pulse is placed in front of the
c simulation box at a distance "zm". The pulse 
c then propagates into the simulation box
c from the left edge. 


c  COORDINATE SYSTEM

c                                    ^ z       xm
c                                    |<----------------->
c                                    |
c                                    |             laser pulse
c            simulation box          |
c                                    |                | | |
c                                    |  ^   <-----    | | |
c                                    |  | zm          | | |
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


      dxm=1.5*1.0e-6                       !modify for pulse setup
      dym=1.5*1.0e-6
      dzm=1.5*1.0e-6
      xm=0.05*1.0e-6
      ym=0.05*1.0e-6
      zm=5.0*1.0e-6


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
         s_pulse_x2_=dcos(-xr)
     &           *dexp(-(xr/dxm)**2)
     &           *dexp(-(yr/dym)**2)
     &           *dexp(-(zr/dzm)**2)
      else if (xm.le.xl) then
         s_pulse_x2_=dcos(xr)
     &           *dexp(-(xr/dxm)**2)
     &           *dexp(-(yr/dym)**2)
     &           *dexp(-(zr/dzm)**2)
      endif


c Turn pulse off


      s_pulse_x2_=0.0d0


      end function s_pulse_x2_
