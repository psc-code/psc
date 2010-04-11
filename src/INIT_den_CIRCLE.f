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
      real(kind=8) x,y,z
      real(kind=8) x1,y1,z1
      real(kind=8) x2,y2,z2
      real(kind=8) w1,w2
      real(kind=8) L1,L2
      real(kind=8) r1,r2
      real(kind=8) phi,phi1,phi2
      real(kind=8) dphi1,dphi2
      real(kind=8) Lphi1,Lphi2
      real(kind=8) arg1,arg2,arg3,arg4
      real(kind=8) INIT_den

c LOCATION OF CIRCLES
c x1: location of density center in x in m of circle 1
c y1: location of density center in y in m of circle 1
c z1: location of density center in z in m of circle 1
c r1: radius of circle 1 in m
c w1: width of circle 1 in m
c L1: density gradient of circle 1 in m
c dphi1: sector size of circle 1 in degree
c x2: location of density center in x in m of circle 2
c y2: location of density center in y in m of circle 2
c z2: location of density center in z in m of circle 2
c w2: width of circle 2 in m
c L2: density gradient of circle 2 in m
c r2: radius of circle 2 in m
c dphi2: sector size of circle 2 in degree

      x1=15.0*1.0e-6/ld
      y1=15.0*1.0e-6/ld
      z1=5060.0*1.0e-6/ld
      r1=4000.0*1.0e-6/ld
      w1=2.5*1.0e-6/ld
      L1=0.1*1.0e-6/ld
      phi1=180.0*3.141592/180.0
      dphi1=90.0*3.141592/180.0
      Lphi1=0.1*3.141592/180.0

      x2=15.0*1.0e-6/ld
      y2=15.0*1.0e-6/ld
      z2=30.0*1.0e-6/ld
      r2=25.0*1.0e-6/ld 
      w2=2.5*1.0e-6/ld 
      L2=0.1*1.0e-6/ld 
      phi2=180.0*3.141592/180.0
      dphi2=90.0*3.141592/180.0
      Lphi2=0.1*3.141592/180.0

c DENSITY SHAPE FUNCTION

      arg1=(abs(sqrt((x-x1)**2+(y-y1)**2+(z-z1)**2)-r1)-w1)/L1
      arg2=(abs(sqrt((x-x2)**2+(y-y2)**2+(z-z2)**2)-r2)-w2)/L2

      phi=acos((z-z1)/sqrt((x-x1)**2+(y-y1)**2+(z-z1)**2))
      arg3=(abs(phi1-phi)-dphi1)/Lphi1
      phi=acos((z-z2)/sqrt((x-x2)**2+(y-y2)**2+(z-z2)**2))
      arg4=(abs(phi2-phi)-dphi2)/Lphi2

      if (arg1.gt.1.0e2) arg1=1.0e2
      if (arg2.gt.1.0e2) arg2=1.0e2
      if (arg3.gt.1.0e2) arg3=1.0e2
      if (arg4.gt.1.0e2) arg4=1.0e2

      INIT_den=1.0/(1.0+dexp(arg1))/(1.0+dexp(arg3))
     &         +1.0/(1.0+dexp(arg2))/(1.0+dexp(arg4))

      end function INIT_den
