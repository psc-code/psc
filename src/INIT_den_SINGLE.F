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
      real(kind=8) INIT_den

      INIT_den=0.0d0
      IF (int(z/dz)==int(0.5*(i3x+i3n))) then
         IF (int(y/dy)==int(0.5*(i2x+i2n))) then
            IF (int(x/dx)==int(0.5*(i1x+i1n))) then
               INIT_den=1.0d0
            endif
         endif
      endif


      end function INIT_den
