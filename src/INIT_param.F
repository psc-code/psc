
#include "config.h"

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

c PARAMETER INITIALIZATION 

c     ----------------------------------------------------------------------
      subroutine INIT_param_fortran
c     ----------------------------------------------------------------------
c     setup processor map

      use VLA_variables
      use PSC_domain

      implicit none

      if (use_pml) call PML_coeff_init_param

      end subroutine INIT_param_fortran

