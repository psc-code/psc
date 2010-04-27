
#include <config.h>


c THIS SUBROUTINE IS THE PARTICLE MOVER.


      subroutine PIC_move_part_yz

      use PIC_variables
      use VLA_variables

      implicit none

      real(kind=8) :: s_cpub
      real(kind=8) :: s_cpuh

      s_cpub=0.0d0
      s_cpuh=0.0d0


c      call C_dump_particles(niloc, p_niloc)
c      call C_dump_fields(ex)
c      call C_dump_fields(ey)
c      call C_dump_fields(ez)
c      call C_dump_fields(bx)
c      call C_dump_fields(by)
c      call C_dump_fields(bz)
      call SERV_systime(cpug)

#ifdef USE_C_INTERFACE
      call C_push_part_yz(cori, alpha, eta, 
     &     p2A, p2B, niloc, p_niloc, 99)
#else
      call PIC_push_part_yz(niloc, p_niloc, cori, alpha, eta, 
     &     p2A, p2B, dt, dx, dy, dz,
     &     i1mn, i2mn, i3mn, i1mx, i2mx, i3mx, rd1, rd2, rd3,
     &     ne, ni, nn, jxi, jyi, jzi, ex, ey, ez, bx, by, bz)
#endif
     
      call SERV_systime(cpua)

c      call C_dump_particles(niloc, p_niloc)
c      call C_dump_fields(ne)
c      call C_dump_fields(ni)
c      call C_dump_fields(nn)
c      call C_dump_fields(jxi)
c      call C_dump_fields(jyi)
c      call C_dump_fields(jzi)

      call PIC_fay(jxi)
      call PIC_faz(jxi)
      call PIC_fay(jyi)
      call PIC_faz(jyi)
      call PIC_fay(jzi)
      call PIC_faz(jzi)

      call PIC_fey(jxi)
      call PIC_fez(jxi)
      call PIC_fey(jyi)
      call PIC_fez(jyi)
      call PIC_fey(jzi)
      call PIC_fez(jzi)

      call PIC_fay(ne)
      call PIC_faz(ne)
      call PIC_fay(ni)
      call PIC_faz(ni)
      call PIC_fay(nn)
      call PIC_faz(nn)

      call PIC_fey(ne)
      call PIC_fez(ne)
      call PIC_fey(ni)
      call PIC_fez(ni)
      call PIC_fey(nn)
      call PIC_fez(nn)

      call PIC_pey
      call PIC_pez
      call SERV_systime(cpub)
      s_cpub=s_cpub+cpub-cpua


      call SERV_systime(cpuh)
      s_cpuh=s_cpuh+cpuh-cpug


      cpumess=cpumess+s_cpub
      cpucomp=cpucomp+s_cpuh-s_cpub


      end subroutine PIC_move_part_yz


