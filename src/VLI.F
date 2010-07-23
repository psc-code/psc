
#include "config.h"

c VLI.F is the init program and VLA.F the main program!

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

      integer :: part_label_off,niloc_new
      integer :: rd1n,rd1x,rd2n,rd2x,rd3n,rd3x

#ifdef USE_C_INTERFACE
      call SERV_systime(cpus)

      call C_init_partition
      call C_init_vars(dt, dx, dy, dz,
     &     i1mn, i2mn, i3mn, i1mx, i2mx, i3mx, rd1, rd2, rd3,
     &     i1n, i2n, i3n, i1x, i2x, i3x,
     &     ne, ni, nn, jxi, jyi, jzi, ex, ey, ez, bx, by, bz,
     &     cori, alpha, beta, eta, wl, wp,
     &     n,
     &     99)
      call C_init_particles
#else
      call INIT_basic
      call INIT_param
      call INIT_param_fortran

      call SERV_systime(cpus)
      include './SERV_consist.F'

      ! I don't understand why the partition may calculate different
      ! boundary sizes (rd..) to be filled with particles, but let's
      ! keep it as-is for now.
      call INIT_partition(part_label_off,rd1n,rd1x,rd2n,rd2x,rd3n,rd3x,
     *     niloc_new)
      call ALLOC_particles(niloc_new)
      call INIT_idistr(part_label_off,rd1n,rd1x,rd2n,rd2x,rd3n,rd3x)
#endif

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
