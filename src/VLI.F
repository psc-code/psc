
#include "config.h"

c VLI.F is the init program and VLA.F the main program!

#ifndef USE_C_INTERFACE
c     ----------------------------------------------------------------------
      program VLI
c     ----------------------------------------------------------------------

      use VLA_variables
      use PIC_variables

      implicit none
      include 'mpif.h'

      integer :: part_label_off,niloc_new
      integer :: rd1n,rd1x,rd2n,rd2x,rd3n,rd3x

      call MPI_INIT(info)

      call INIT_basic
      call INIT_param
      call INIT_param_fortran
      include './SERV_consist.F'

      ! I don't understand why the partition may calculate different
      ! boundary sizes (rd..) to be filled with particles, but let's
      ! keep it as-is for now.
      call INIT_partition(part_label_off,rd1n,rd1x,rd2n,rd2x,rd3n,rd3x,
     *     niloc_new)
      call ALLOC_particles(niloc_new)
      call INIT_idistr(part_label_off,rd1n,rd1x,rd2n,rd2x,rd3n,rd3x)
      call ALLOC_field
      call ALLOC_field_fortran
      call INIT_field
      call SETUP_field

c      call INIT_MCC

      call PSC_driver

      call MPI_FINALIZE(info)

      end
#endif
