
#include "config.h"

      subroutine ALLOC_particles(niloc_)

      use PIC_variables
      use VLA_variables

      implicit none

      integer :: niloc_
      real(kind=8) :: niloc_overfac

      niloc_overfac = 1.2

      niloc = niloc_
      nialloc=nint(niloc_overfac*niloc+10)
      allocate(p_niloc(0:11*nialloc+10))

#ifdef USE_C_INTERFACE
      call C_set_particle_info(niloc, p_niloc)
#endif

      end subroutine ALLOC_particles
 
