c SYSTEM TIME IN SECONDS.

      subroutine SERV_systime(cpu)
      implicit none
      include 'mpif.h'

      real(kind=8) :: cpu

      cpu=MPI_Wtime()

      end subroutine SERV_systime
