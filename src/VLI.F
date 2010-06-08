
#include "config.h"

c VLI.f is the init program and VLA.f the main program!

#ifndef USE_C_INTERFACE
      program VLI

      use VLA_variables
      use PIC_variables

      implicit none
      include 'mpif.h'

      call MPI_INIT(info)
      call VLI_main
      call MPI_FINALIZE(info)

      end
#endif


c The real VLI main

      subroutine VLI_main

      use VLA_variables
      use PIC_variables

      implicit none
      include 'mpif.h'

      call MPI_COMM_SIZE(MPI_COMM_WORLD,npe,info)
      call MPI_COMM_RANK(MPI_COMM_WORLD,mpe,info)

      call INIT_param
      call SERV_systime(cpus)
      include './SERV_consist.F'


      call INIT_idistr
      call ALLOC_field
#ifdef USE_C_INTERFACE
      call C_INIT_field
#else
      call INIT_field
#endif
      call SETUP_field
c      call INIT_MCC

      call PSC_driver

      return
      end
