c THIS SUBROUTINE IS THE PARTICLE MOVER.


      subroutine PIC_move_part_z

      use PIC_variables
      use VLA_variables

      implicit none

      real(kind=8) :: s_cpub
      real(kind=8) :: s_cpuh

      s_cpub=0.0d0
      s_cpuh=0.0d0
      call SERV_systime(cpug)

      call PIC_push_part_z(niloc, p_niloc, cori, alpha, eta, 
     &     wl, wp,
     &     p2A, p2B, n, dt, dx, dy, dz,
     &     i1mn, i2mn, i3mn, i1mx, i2mx, i3mx, rd1, rd2, rd3,
     &     ne, ni, nn, jxi, jyi, jzi, ex, ey, ez, bx, by, bz)

      call SERV_systime(cpua)
      call PIC_faz(jxi)
      call PIC_faz(jyi)
      call PIC_faz(jzi)

      call PIC_fez(jxi)
      call PIC_fez(jyi)
      call PIC_fez(jzi)

      call PIC_faz(ne)
      call PIC_faz(ni)
      call PIC_faz(nn)

      call PIC_fez(ne)
      call PIC_fez(ni)
      call PIC_fez(nn)

      call PIC_pez
      call SERV_systime(cpub)
      s_cpub=s_cpub+cpub-cpua


      call SERV_systime(cpuh)
      s_cpuh=s_cpuh+cpuh-cpug


      cpumess=cpumess+s_cpub
      cpucomp=cpucomp+s_cpuh-s_cpub


      end subroutine PIC_move_part_z
