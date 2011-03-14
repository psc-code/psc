
#include "config.h"

      module VLA_variables
      implicit none

      real(kind=8) dx2A,dy2A,dz2A,bx2A,by2A,bz2A  ! ab  
      real(kind=8) dx2B,dy2B,dz2B,bx2B,by2B,bz2B  ! ab
      real(kind=8) poyx1i,poyx1o,poyx2i,poyx2o                 ! ab
      real(kind=8) poyy1i,poyy1o,poyy2i,poyy2o                 ! ab
      real(kind=8) poyz1i,poyz1o,poyz2i,poyz2o                 ! ab
      real(kind=8) poyni, poyno                 ! ab

c nstart: time step counter for restarting the code
c shift_c: time step counter for moving the frame in the code
c shift_l: step counter for moving the frame in the code
c nnp: number of timesteps for one full cycle
c np:  number of timesteps for time-averaging


      integer nstart,shift_c,shift_l,np,nnp

      real(kind=8) fluxi,fluxo
      real(kind=8) ex2A,ey2A,ez2A,hx2A,hy2A,hz2A
      real(kind=8) ex2B,ey2B,ez2B,hx2B,hy2B,hz2B
      real(kind=8) pox,poy,poz,je
      real(kind=8) p2A,p2B

c pi: the number pi
c shift_z: density function shift in dimensionless units

      real(kind=8) pi
      real(kind=8) shift_z

      integer n

      end module VLA_variables

c     ======================================================================
      subroutine PSC_set_timestep(n_)
c     ======================================================================

      use VLA_variables, only: n
      implicit none

      integer :: n_

      n = n_

      end subroutine PSC_set_timestep

