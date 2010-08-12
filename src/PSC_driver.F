
#include <config.h>

c     this contains the main driver for the code, shared between VLI and VLA

      subroutine PSC_driver

      use VLA_variables
      use PIC_variables
      use FIELDS

      implicit none
      include 'mpif.h'

      call OUT_param
      call MPI_BARRIER(MPI_COMM_WORLD,info)

      do n=nstart,nmax

         cpuinou=0.0d0
         cpumess=0.0d0
         cpucomp=0.0d0

c         call PIC_bpulse      ! particle injection
c         call PIC_randomize(niloc, p_niloc)
c         call PIC_find_cell_indices(niloc, p_niloc)
c         call PIC_sort(niloc, p_niloc)
c         call PIC_ionize      ! field ionization
c         call MCC_impact      ! impact ionization
c         call PIC_bin_coll    ! binary collisions

         call OUT_field       ! field output at t=n*dt
         call OUT_part        ! particle output at t=n*dt
c         call OUT_poyc        ! energy conservation

         ! field propagation n*dt -> (n+0.5)*dt
         if (use_pml) then
            call PIC_pml_msa
         else
            call PIC_msa(ex, ey, ez, bx, by, bz, jxi, jyi, jzi)
         endif

         call PIC_move_part   ! particle propagation n*dt -> (n+1.0)*dt  

         ! field propagation (n+0.5)*dt -> (n+1.0)*dt
         if (use_pml) then
            call PIC_pml_msb
         else
            call PIC_msb(ex, ey, ez, bx, by, bz, jxi, jyi, jzi)    
         endif

         call MPI_BARRIER(MPI_COMM_WORLD,info)

         call SERV_systime(cpuf)
         include 'SERV_cput.F'
         include 'OUT_cput.F'


         if (cpue.ge.cpum) then
            call SERV_write(n, niloc, p_niloc,
     &           ne, ni, nn, jxi, jyi, jzi, ex, ey, ez, hx, hy, hz)
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

      call SERV_write(nmax, niloc, p_niloc,
     &     ne, ni, nn, jxi, jyi, jzi, ex, ey, ez, hx, hy, hz)
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
