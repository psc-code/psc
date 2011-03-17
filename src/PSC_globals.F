
c     ======================================================================
      module PSC_globals
c     ======================================================================
      implicit none

      real(kind=8) :: cori       ! inverse particles / cell
      real(kind=8) :: alpha, eta ! normalization constants
      
      end module PSC_globals

c     ======================================================================
      subroutine PSC_set_globals(cori_, alpha_, eta_)
c     ======================================================================
      use PSC_globals
      implicit none

      real(kind=8) :: cori_, alpha_, eta_

      cori = cori_
      alpha = alpha_
      eta = eta_

      end subroutine PSC_set_globals

