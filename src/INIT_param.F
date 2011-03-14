
c     ----------------------------------------------------------------------
      subroutine INIT_basic
c     ----------------------------------------------------------------------

      use VLA_variables, only: pi
      use PSC_domain, only: npe, mpe

      implicit none
      include 'mpif.h'

      integer :: ierr

c PARAMETERS THAT SHOULD NOT BE CHANGED UNDER NORMAL CONDITIONS

c pi: pi

      pi=3.1415927

      call MPI_COMM_SIZE(MPI_COMM_WORLD,npe,ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,mpe,ierr)

      end subroutine INIT_basic

