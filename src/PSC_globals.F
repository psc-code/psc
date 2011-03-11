
c     ======================================================================
      module PSC_globals
c     ======================================================================
      implicit none

      real(kind=8) :: cori
      
      end module PSC_globals

c     ======================================================================
      subroutine PSC_set_globals(cori_)
c     ======================================================================
      use PSC_globals
      implicit none

      real(kind=8) :: cori_

      cori = cori_

      end subroutine PSC_set_globals
