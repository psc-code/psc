
#include "config.h"

c     ======================================================================
      subroutine PIC_sort
c     ======================================================================
c     THIS SUBROUTINE SORTS THE PARTICLES TO CELLS


      use PIC_variables
      use VLA_variables

      implicit none

      real(kind=8) :: s_cpuh


      s_cpuh=0.0d0
      call SERV_systime(cpug)

      call PIC_randomize(niloc, p_niloc)
      call PIC_find_cell_indices(niloc, p_niloc)
#ifdef USE_C_INTERFACE
      call C_sort(niloc, p_niloc, 99)
#else
      call PIC_sort_1(niloc, p_niloc)
#endif

      call SERV_systime(cpuh)
      s_cpuh=s_cpuh+cpuh-cpug

      cpucomp=cpucomp+s_cpuh

      end subroutine PIC_sort
