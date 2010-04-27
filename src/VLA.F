
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
      call OUT_param
#ifdef USE_C_INTERFACE
      call C_init_vars(dt, dx, dy, dz,
     &     i1mn, i2mn, i3mn, i1mx, i2mx, i3mx, rd1, rd2, rd3,
     &     ne, ni, nn, jxi, jyi, jzi, ex, ey, ez, bx, by, bz,
     &     cori, alpha, eta, wl, wp,
     &     99)
#endif
      call MPI_BARRIER(MPI_COMM_WORLD,info)


c BEGINNING OF TIME LOOP


      do n=nstart,nmax

         cpuinou=0.0d0
         cpumess=0.0d0
         cpucomp=0.0d0

c         call PIC_bpulse      ! particle injection
         call PIC_sort        ! particle randomization       by ab
c         call PIC_ionize       ! field ionization
c         call PIC_bin_coll    ! binary collisions
c         call MCC_impact      ! impact ionization

         call OUT_field       ! field output at t=n*dt
         call OUT_part        ! particle output at t=n*dt
         call OUT_poyc        ! energy conservation

         call PIC_pml_msa         ! field propagation n*dt -> (n+0.5)*dt
         call PIC_move_part   ! particle propagation n*dt -> (n+1.0)*dt  by ab
         call PIC_pml_msb         ! field propagation (n+0.5)*dt -> (n+1.0)*dt

c         call PIC_msa         ! field propagation n*dt -> (n+0.5)*dt
c         call PIC_move_part   ! particle propagation n*dt -> (n+1.0)*dt
c         call PIC_msb         ! field propagation (n+0.5)*dt -> (n+1.0)*dt


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

      call MPI_FINALIZE(info)
      stop

      end
