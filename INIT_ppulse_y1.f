c Laser pulse initialization (p-polarization).

c Dies ist der blaue Probe von links
c Beachte: Intensität zurück auf 0.1 setzen !!

      function p_pulse_y1(x,y,z,t)

      use VLA_variables
      use PIC_variables

      implicit none
      real(kind=8) :: t,x,y,z
      real(kind=8) :: xm,ym,zm
      real(kind=8) :: dxm,dym,dzm
      real(kind=8) :: widthz
      real(kind=8) :: xl,yl,zl
      real(kind=8) :: xr,yr,zr
      real(kind=8) :: p_pulse_y1


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
c widthz: width of pulse in z in m
c xm: x-location of pulse center in m
c ym: y-location of pulse center in m
c zm: z-location of pulse center in m


c      dxm=10.0*1.0e-6                       !modify for pulse setup
c      dym=8.0*0.255*1.0e-6                  !10 fs Puls hier, 1 fs = 0,255 mum
c      dzm=10.0*1.0e-6                       ! vormals dxm dzm = 6.0
c      xm=10.0*1.0e-6
c      ym=0.0*1.0e-6
c      zm=10.0*1.0e-6
c      widthz=6.0*1.0e-6

      dxm=1.5*1.0e-6            !modify for pulse setup
      dym=1.5*1.0e-6
      dzm=1.5*1.0e-6
      xm=5.0*1.0e-6
      ym=-4.0*1.0e-6
      zm=5.0*1.0e-6

      xm=xm/ld                              !normalization
      ym=ym/ld
      zm=zm/ld
      dxm=dxm/ld
      dym=dym/ld
      dzm=dzm/ld
      widthz=widthz/ld


      xl=x
      yl=y-t
      zl=z


      xr=xl-xm
      yr=yl-ym
      zr=zl-zm


      if (yl.gt.ym) then
c         p_pulse_y1=0.05*dsin(2.0*yr)             ! Frequency doubled and reduced intensity
         p_pulse_y1=dsin(yr)
     &           *dexp(-(xr/dxm)**2)
     &           *dexp(-(yr/dym)**2)
     &           *dexp(-(zr/dzm)**2)
c     &           *1.0/(1.0+dexp((abs(zr)-widthz)/dzm)**2)
      else if (ym.ge.yl) then
c          p_pulse_y1=0.05*dsin(2.0*yr)
         p_pulse_y1=dsin(yr)
     &           *dexp(-(xr/dxm)**2)
     &           *dexp(-(yr/dym)**2)            ! diese Zeile auskommentieren, wenn Puls anbleiben soll
     &           *dexp(-(zr/dzm)**2)
c     &           *1.0/(1.0+dexp((abs(zr)-widthz)/dzm)**2)
      endif

!      pos_y1 = abs(ym)+3.0*dym+thick*dy
      pos_y1 = thick*dy+3.0*dym-ym      ! ab

c Turn pulse off

      p_pulse_y1=0.0d0

      end function p_pulse_y1
