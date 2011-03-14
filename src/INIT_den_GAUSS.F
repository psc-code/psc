c DENSITY PROFILE SETUP. TO GUARANTEE THE CORRECT NORMALIZATION
c IN THE CODE THE CONDIDTION 0.0 <= INIT_den <= 1.0 MUST BE
c FULLFILLED. TO LOCATE THE DENSITY PROFILE IN THE SIMULATION
c BOX A COORDINATE FRAME IS ADOPTED WHOSE ORIGIN IS THE LOWER
c LEFT CORNER OF THE SIMULATION BOX. THE DENSITY PROFILE IS
c DEFINED RELATIVE TO THIS FRAME. THE CODE USES THE OVERLAP 
c BETWEEN DENSITY PROFILE AND SIMULATION BOX.


      function INIT_den(x,y,z)

      implicit none 
      real(kind=8) x,y,z
      real(kind=8) x0,y0,z0,Lx,Ly,Lz
      real(kind=8) rrx,rry,rrz
      real(kind=8) INIT_den


c x0: location of density center in x in m
c y0: location of density center in y in m
c z0: location of density center in z in m
c Lx: gradient of density profile in x in m
c Ly: gradient of density profile in y in m
c Lz: gradient of density profile in z in m


      x0=5.0*1.0e-6
      y0=5.0*1.0e-6
      z0=25.0*1.0e-6
      Lx=2.5*1.0e-6
      Ly=2.5*1.0e-6 
      Lz=10.0*1.0e-6 


c NORMALIZATION


      x0=x0/ld
      y0=y0/ld
      z0=z0/ld
      Lx=Lx/ld
      Ly=Ly/ld
      Lz=Lz/ld


      rrx=(x-x0)*(x-x0)/(Lx*Lx)
      rry=(y-y0)*(y-y0)/(Ly*Ly)
      rrz=(z-z0)*(z-z0)/(Lz*Lz)


      INIT_den=exp(-rrx)*exp(-rry)*exp(-rrz)


      end function INIT_den
