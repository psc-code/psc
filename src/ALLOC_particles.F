
#include "config.h"

c     ======================================================================
      subroutine ALLOC_particles(niloc_)
c     ======================================================================

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
      call C_ALLOC_particles_cb(p_niloc)
#endif

      end subroutine ALLOC_particles
 
c     ======================================================================
      subroutine REALLOC_particles(niloc_new)
c     ======================================================================
c     enlarge particle array to size niloc_new
c     reallocating if necessary

      use PIC_variables
      use VLA_variables

      implicit none

      integer :: niloc_new

      integer :: la,le
      character*(5) :: node
      real(kind=8) :: niloc_overfac

      niloc_overfac = 1.2

      if (niloc_new.gt.nialloc) then
         call SERV_labelgen(mpe,node)
         call SERV_systime(cpuc)
         open(11,file=trim(data_out)//'/'//node//'ENLARGE',
     &        access='sequential',form='unformatted')
         do l=0,11*niloc+10,100
            le=min(l+99,11*niloc+10)
            write(11) (p_niloc(la),la=l,le)
         enddo
         close(11)

         nialloc=int(niloc_overfac*niloc_new+10)
         deallocate(p_niloc)
         allocate(p_niloc(0:11*nialloc+10))

         open(11,file=trim(data_out)//'/'//node//'ENLARGE',
     &        access='sequential',form='unformatted')
         do l=0,11*niloc+10,100
            le=min(l+99,11*niloc+10)
            read(11) (p_niloc(la),la=l,le)
         enddo
         close(11)
         call SERV_systime(cpud)
         cpuinou=cpuinou+cpud-cpuc
      endif

      end subroutine REALLOC_particles
