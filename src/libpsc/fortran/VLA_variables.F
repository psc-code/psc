
      module VLA_variables
      implicit none

c pi: the number pi

      real(kind=8) pi

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

