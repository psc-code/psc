c DENSITY PROFILE SETUP. TO GUARANTEE THE CORRECT NORMALIZATION
c IN THE CODE THE CONDIDTION 0.0 <= INIT_den <= 1.0 MUST BE
c FULLFILLED. TO LOCATE THE DENSITY PROFILE IN THE SIMULATION
c BOX A COORDINATE FRAME IS ADOPTED WHOSE ORIGIN IS THE LOWER
c LEFT CORNER OF THE SIMULATION BOX. THE DENSITY PROFILE IS
c DEFINED RELATIVE TO THIS FRAME. THE CODE USES THE OVERLAP 
c BETWEEN DENSITY PROFILE AND SIMULATION BOX.


      function INIT_den(x,y,z)

      use VLA_variables
      use PIC_variables

      implicit none 
      real(kind=8) x,y,z,xr,yr,zr,rot
      real(kind=8) x0,y0,z0,Lx,Ly,Lz
      real(kind=8) widthx,widthy,widthz
      real(kind=8) argx,argy,argz
      real(kind=8) INIT_den
      real(kind=8) xmin,xmax,ymin,ymax,zmin,zmax   ! added by ab


c x0: location of density center in x in m
c y0: location of density center in y in m
c z0: location of density center in z in m
c Lx: gradient of density profile in x in m
c Ly: gradient of density profile in y in m
c Lz: gradient of density profile in z in m
c widthx: width of transverse density profile in m
c widthy: width of transverse density profile in m
c widthz: width of longitudinal density profile in m


      x0=5.0*1.0e-6
      y0=5.0*1.0e-6
      z0=5.0*1.0e-6
      Lx=0.1*1.0e-6
      Ly=1.5*1.0e-6
      Lz=1.5*1.0e-6
      widthx=1.0*1.0e-6
      widthy=0.0*1.0e-6
      widthz=0.0*1.0e-6
      rot=0.0*3.141592/180.0


c NORMALIZATION


      x0=x0/ld
      y0=y0/ld
      z0=z0/ld
      Lx=Lx/ld
      Ly=Ly/ld
      Lz=Lz/ld
      widthx=widthx/ld
      widthy=widthy/ld
      widthz=widthz/ld


      xr=x
      yr=dcos(rot)*(y-y0)-sin(rot)*(z-z0)+y0
      zr=dcos(rot)*(z-z0)+sin(rot)*(y-y0)+z0

      argx=(abs(xr-x0)-widthx)/Lx
      argy=(abs(yr-y0)-widthy)/Ly
      argz=(abs(zr-z0)-widthz)/Lz
      if(argx>200.0) argx=200.0
      if(argy>200.0) argy=200.0
      if(argz>200.0) argz=200.0

c      INIT_den=1.0/(1.0+dexp(argx))
c     &         /(1.0+dexp(argy))/(1.0+dexp(argz))
      INIT_den=dexp(-argy*argy)*dexp(-argz*argz)

          
c Switch density off


c      INIT_den=0.0


      end function INIT_den
