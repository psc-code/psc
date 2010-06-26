
#include <config.h>

c VLI.F is the init program and VLA.F the main program!

#ifndef USE_C_INTERFACE
      program VLA

      use VLA_variables
      use PIC_variables

      implicit none
      include 'mpif.h'

      call MPI_INIT(info)
      call VLA_main
      call MPI_FINALIZE(info)

      end
#endif


c The real VLA main

      subroutine VLA_main

      use VLA_variables
      use PIC_variables

      implicit none

      call INIT_basic
      call INIT_param
      call INIT_param_fortran
      include './SERV_consist.F'

      call SERV_read
c      call INIT_MCC

      call PSC_driver
      
      return
      end
