c Laser pulse initialization (p-polarization).

      function p_pulse_z1(x,y,z,t)

      implicit none
      real(kind=8) :: t,x,y,z
      real(kind=8) :: xm,ym,zm
      real(kind=8) :: rr,rm,drm,dtm
      real(kind=8) :: th,sph,cph
      real(kind=8) :: xl,yl,zl
      real(kind=8) :: p_pulse_z1


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

c drm: width of pulse in radial direction in m
c dtm: width of pulse in angular direction in m
c xm: x-location of geometrical pulse focus in m
c ym: y-location of geometrical pulse focus in m
c zm: z-location of geometrical pulse focus in m


      dtm=0.1*pi                           !modify for pulse setup
      drm=2.0*1.0e-6
      xm=10.0*1.0e-6
      ym=10.0*1.0e-6
      zm=27.0*1.0e-6
      rm=29.0*1.0e-6


      xm=xm/ld                              !normalization
      ym=ym/ld
      zm=zm/ld
      rm=rm/ld
      drm=drm/ld


      if (z/=zm) then
         th=atan(sqrt((x-xm)**2+(y-ym)**2)/abs(z-zm))
      else
         th=0.5*pi
      endif
      cph=(x-xm)/sqrt((x-xm)**2+(y-ym)**2)
      sph=(y-ym)/sqrt((x-xm)**2+(y-ym)**2)
        

      xl=x+t*sin(th)*cph
      yl=y+t*sin(th)*sph
      zl=z-t*cos(th)
      rr=sqrt((xl-xm)**2+(yl-ym)**2+(zl-zm)**2)


      if (rr.lt.rm) then
         p_pulse_z1=dcos(rr-rm)
     &           *dexp(-((rr-rm)/drm)**2)
     &           *dexp(-(th/dtm)**2)
      else if (rr.ge.rm) then
         p_pulse_z1=dcos(rr-rm)
     &           *dexp(-((rr-rm)/drm)**2)
     &           *dexp(-(th/dtm)**2)
      endif


c Turn pulse off


c      p_pulse_z1=0.0d0


      end function p_pulse_z1
