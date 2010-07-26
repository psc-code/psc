
#include <config.h>

c     this contains the main driver for the code, shared between VLI and VLA

      subroutine PSC_driver

      use VLA_variables
      use PIC_variables

      implicit none
      include 'mpif.h'

      call OUT_param
      call MPI_BARRIER(MPI_COMM_WORLD,info)

      do n=nstart,nmax

         cpuinou=0.0d0
         cpumess=0.0d0
         cpucomp=0.0d0

#ifdef USE_C_INTERFACE
         call C_init_vars(dt, dx, dy, dz,
     &        i1mn, i2mn, i3mn, i1mx, i2mx, i3mx, rd1, rd2, rd3,
     &        i1n, i2n, i3n, i1x, i2x, i3x,
     &        ne, ni, nn, jxi, jyi, jzi, ex, ey, ez, bx, by, bz,
     &        cori, alpha, beta, eta, wl, wp,
     &        n,
     &        99)
#endif

c         call PIC_bpulse      ! particle injection
c         call PIC_sort        ! particle randomization
c         call PIC_ionize      ! field ionization
c         call MCC_impact      ! impact ionization
c         call PIC_bin_coll    ! binary collisions

         call OUT_field       ! field output at t=n*dt
         call OUT_part        ! particle output at t=n*dt
c         call OUT_poyc        ! energy conservation

#ifdef USE_PML
         call PIC_pml_msa     ! field propagation n*dt -> (n+0.5)*dt       !ab
         call PIC_move_part   ! particle propagation n*dt -> (n+1.0)*dt  
         call PIC_pml_msb     ! field propagation (n+0.5)*dt -> (n+1.0)*dt !ab
#else
         call PIC_msa         ! field propagation n*dt -> (n+0.5)*dt
         call PIC_move_part   ! particle propagation n*dt -> (n+1.0)*dt
         call PIC_msb         ! field propagation (n+0.5)*dt -> (n+1.0)*dt
#endif

         call MPI_BARRIER(MPI_COMM_WORLD,info)

         call SERV_systime(cpuf)
         include 'SERV_cput.F'
         include 'OUT_cput.F'


         if (cpue.ge.cpum) then
            call SERV_write(n)
            call MPI_BARRIER(MPI_COMM_WORLD,info)
 
           if (mpe.eq.0) then
               open(11,file='./okfile',
     &              access='sequential',form='formatted')
               write(11,*) 'Partial job has finished regularly!'
               close(11)
            endif
 
            call MPI_FINALIZE(info)
            stop
          endif

      enddo


c END OF TIME LOOP

      call SERV_write(nmax)
      call MPI_BARRIER(MPI_COMM_WORLD,info)

      if (mpe.eq.0) then
         open(11,file='./okfile',
     &        access='sequential',form='formatted')
         write(11,*) 'Partial job has finished regularly!'
         close(11)
         open(11,file='./endfile',
     &        access='sequential',form='formatted')
         write(11,*) 'Job has finished regularly!'
         close(11)
      endif

      return
      
      end subroutine PSC_driver
