
c Laser pulse initialization (p-polarization).

      function p_pulse_x1_(x,y,z,t)

      use PSC_params
      use PSC_patch
      use PSC_domain
      use COEFF_pml

      implicit none
      real(kind=8) :: t,x,y,z
      real(kind=8) :: xm,ym,zm,xb
      real(kind=8) :: dxm,dym,dzm
      real(kind=8) :: xl,yl,zl
      real(kind=8) :: xr,yr,zr
      real(kind=8) :: p_pulse_x1_


c NOTE: The pulse is placed in front of the
c simulation box at a distance "zm". The pulse 
c then propagates into the simulation box
c from the left edge. 


c  COORDINATE SYSTEM

c                          xm        ^ y
c                 <----------------->|
c                                    |
c            laser pulse             |
c                                    |     simulation
c               | | |                |     box
c               | | |   ------   ^   |
c               | | |         ym |   |
c                                |   |
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


      dxm=0.5*1.0e-6                       !modify for pulse setup
      dym=1.0*1.0e-6
      dzm=1.0*1.0e-6
      xb=1.0*1.0e-6
      xm=-2.0*1.0e-6
      ym=5.0*1.0e-6
      zm=5.0*1.0e-6


      xm=xm/ld                              !normalization
      ym=ym/ld
      zm=zm/ld
      xb=xb/ld
      dxm=dxm/ld
      dym=dym/ld
      dzm=dzm/ld


      xl=x-t
      yl=y
      zl=z


      xr=xl-xm
      yr=yl-ym
      zr=zl-zm


      if (xl.gt.xm) then
         p_pulse_x1_=dsin(xr)
     &           *1.0/(1.0+dexp((abs(xr)-xb)/dxm))
     &           *dexp(-(yr/dym)**2)
     &           *dexp(-(zr/dzm)**2)
      else if (xm.ge.xl) then
          p_pulse_x1_=dsin(xr)
     &           *1.0/(1.0+dexp((abs(xr)-xb)/dxm))
     &           *dexp(-(yr/dym)**2)
     &           *dexp(-(zr/dzm)**2)
      endif

      pos_x1 = thick*dx+1.5*xb-xm        ! ab

c Turn pulse off


      p_pulse_x1_=0.0d0


      end function p_pulse_x1_
