c Laser pulse initialization (p-polarization).

      function p_pulse_y1(x,y,z,t)

      implicit none
      real(kind=8) :: t,x,y,z
      real(kind=8) :: p_pulse_y1
      real(kind=8) :: C_p_pulse_y1

      p_pulse_y1 = C_p_pulse_y1(x,y,z,t)

      end function p_pulse_y1

c Dies ist der blaue Probe von links
c Beachte: Intensität zurück auf 0.1 setzen !!

      function p_pulse_y1_(x,y,z,t)

      use PSC_params
      use PSC_patch
      use PSC_domain
      use COEFF_pml

      implicit none
      real(kind=8) :: t,x,y,z
      real(kind=8) :: xm,ym,zm,yb
      real(kind=8) :: dxm,dym,dzm
      real(kind=8) :: xl,yl,zl
      real(kind=8) :: xr,yr,zr
      real(kind=8) :: p_pulse_y1_


c NOTE: The pulse is placed in front of the
c simulation box at a distance "zm". The pulse 
c then propagates into the simulation box
c from the left edge. 


c  COORDINATE SYSTEM

c                          ym        ^ x
c                 <----------------->|
c                                    |
c            laser pulse             |
c                                    |     simulation
c               | | |                |     box
c               | | |   ----->   ^   |
c               | | |         xm |   |
c                                |   |
c          ------------------------------------------------->
c                              (i1n,i2n,i3n)=box origin    y 



c****************************************************************
c SETUP OF SHORT LASER PULSE
c****************************************************************

c dxm: width of pulse in x in m
c dym: width of pulse in y in m
c dzm: gradient of pulse envelope in z in m
c xm: x-location of pulse center in m
c ym: y-location of pulse center in m
c zm: z-location of pulse center in m

      dxm=1.0*1.0e-6            !modify for pulse setup
      dym=0.5*1.0e-6
      dzm=1.0*1.0e-6
      yb=1.0*1.0e-6
      xm=5.0*1.0e-6
      ym=-2.0*1.0e-6
      zm=5.0*1.0e-6

      xm=xm/ld                              !normalization
      ym=ym/ld
      zm=zm/ld
      yb=yb/ld
      dxm=dxm/ld
      dym=dym/ld
      dzm=dzm/ld


      xl=x
      yl=y-t
      zl=z


      xr=xl-xm
      yr=yl-ym
      zr=zl-zm


      if (yl.gt.ym) then
         p_pulse_y1_=dsin(yr)
     &           *dexp(-(xr/dxm)**2)
     &           *1.0/(1.0+dexp((abs(yr)-yb)/dym))
     &           *dexp(-(zr/dzm)**2)
      else if (ym.ge.yl) then
         p_pulse_y1_=dsin(yr)
     &           *dexp(-(xr/dxm)**2)
     &           *1.0/(1.0+dexp((abs(yr)-yb)/dym))
     &           *dexp(-(zr/dzm)**2)
      endif

      pos_y1 = thick*dy+1.5*yb-ym      ! ab

c Turn pulse off


      p_pulse_y1_=0.0d0


      end function p_pulse_y1_
