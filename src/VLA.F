
#include <config.h>

c VLI.f is the init program and VLA.f the main program!

      program VLA

      use VLA_variables
      use PIC_variables

      implicit none
      include 'mpif.h'


c INITIALIZATION


      call MPI_INIT(info)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,npe,info)
      call MPI_COMM_RANK(MPI_COMM_WORLD,mpe,info)


      call INIT_param
      call SERV_systime(cpus)
      include './SERV_consist.F'


      call SERV_read
c      call INIT_MCC

      call PSC_driver

      call MPI_FINALIZE(info)
      stop

      end
