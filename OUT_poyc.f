c THIS SUBROUTINE CHECKS THE POYNTING LAW IN 3D.


      subroutine OUT_poyc

      use PIC_variables
      use VLA_variables

      implicit none

      character*5 label,node

      real(kind=8) :: s_cpud
      real(kind=8) :: s_cpuh


      s_cpud=0.0d0
      s_cpuh=0.0d0
      call SERV_systime(cpug)


c Time-averaging


      if (n.eq.tmnvp) then

         fluxit=fluxit+dt*fluxi
         fluxot=fluxot+dt*fluxo

         enEXt=enEXt+ex2B-ex2A
         enEYt=enEYt+ey2B-ey2A
         enEZt=enEZt+ez2B-ez2A
         enBXt=enBXt+bx2B-bx2A
         enBYt=enBYt+by2B-by2A
         enBZt=enBZt+bz2B-bz2A
         enHXt=enHXt+hx2B-hx2A
         enHYt=enHYt+hy2B-hy2A
         enHZt=enHZt+hz2B-hz2A

!         ent=ent+ex2B-ex2A+ey2B-ey2A+ez2B-ez2A
!     &          +bx2B-bx2A+by2B-by2A+bz2B-bz2A
         ent=ent+ex2B-ex2A+ey2B-ey2A+ez2B-ez2A
     &          +hx2B-hx2A+hy2B-hy2A+hz2B-hz2A

         poxt=poxt+dt*pox
         poyt=poyt+dt*poy
         pozt=pozt+dt*poz

         poynit = poynit+dt*poyni          ! added by ab
         poynot = poynot+dt*poyno          ! added by ab

         jet=jet+dx*dy*dz*(p2B)
c         jet=jet+dt*je

      endif


      if ((n.gt.tmnvp).and.(n.lt.tmxvp)) then
         fluxit=fluxit+dt*fluxi
         fluxot=fluxot+dt*fluxo

         enEXt=enEXt+ex2B-ex2A
         enEYt=enEYt+ey2B-ey2A
         enEZt=enEZt+ez2B-ez2A
         enBXt=enBXt+bx2B-bx2A
         enBYt=enBYt+by2B-by2A
         enBZt=enBZt+bz2B-bz2A
         enHXt=enHXt+hx2B-hx2A
         enHYt=enHYt+hy2B-hy2A
         enHZt=enHZt+hz2B-hz2A

!         ent=ent+ex2B-ex2A+ey2B-ey2A+ez2B-ez2A
!     &          +bx2B-bx2A+by2B-by2A+bz2B-bz2A
         ent=ent+ex2B-ex2A+ey2B-ey2A+ez2B-ez2A
     &          +hx2B-hx2A+hy2B-hy2A+hz2B-hz2A

         poxt=poxt+dt*pox
         poyt=poyt+dt*poy
         pozt=pozt+dt*poz

         poynit = poynit+dt*poyni          ! added by ab
         poynot = poynot+dt*poyno          ! added by ab

         jet=jet+dx*dy*dz*(p2B)
c         jet=jet+dt*je

      endif


      if (n.eq.tmxvp) then
         tmnvp=n+1
         tmxvp=n+np

         fluxit=fluxit+dt*fluxi
         fluxot=fluxot+dt*fluxo

         enEXt=enEXt+ex2B-ex2A
         enEYt=enEYt+ey2B-ey2A
         enEZt=enEZt+ez2B-ez2A
         enBXt=enBXt+bx2B-bx2A
         enBYt=enBYt+by2B-by2A
         enBZt=enBZt+bz2B-bz2A
         enHXt=enHXt+hx2B-hx2A
         enHYt=enHYt+hy2B-hy2A
         enHZt=enHZt+hz2B-hz2A

!         ent=ent+ex2B-ex2A+ey2B-ey2A+ez2B-ez2A
!     &          +bx2B-bx2A+by2B-by2A+bz2B-bz2A
         ent=ent+ex2B-ex2A+ey2B-ey2A+ez2B-ez2A
     &          +hx2B-hx2A+hy2B-hy2A+hz2B-hz2A

         poxt=poxt+dt*pox
         poyt=poyt+dt*poy
         pozt=pozt+dt*poz

         poynit = poynit+dt*poyni          ! added by ab
         poynot = poynot+dt*poyno          ! added by ab

         jet=jet+dx*dy*dz*(p2B)
c         jet=jet+dt*je

         call SERV_labelgen(mpe,node)
         call SERV_labelgen(n,label)

         call SERV_systime(cpuc)

         open(11,file=trim(data_out)//'/'//node//'poynting'//label,
     &        access='sequential',form='unformatted')

         write(11) n
         write(11) i1mn
         write(11) i1mx
         write(11) i2mn
         write(11) i2mx
         write(11) i3mn
         write(11) i3mx
         write(11) fluxit
         write(11) fluxot
         write(11) enEXt
         write(11) enEYt
         write(11) enEZt
         write(11) enBXt
         write(11) enBYt
         write(11) enBZt
         write(11) enHXt
         write(11) enHYt
         write(11) enHZt
         write(11) ent
         write(11) poxt
         write(11) poyt
         write(11) pozt
         write(11) jet
         write(11) poynit               ! ab
         write(11) poynot               ! ab

         close(11)

         call SERV_systime(cpud)
         s_cpud=s_cpud+cpud-cpuc

         fluxit=0.0
         fluxot=0.0

         enEXt=0.0
         enEYt=0.0
         enEZt=0.0
         enBXt=0.0
         enBYt=0.0
         enBZt=0.0
         enHXt=0.0
         enHYt=0.0
         enHZt=0.0

         ent=0.0
         poxt=0.0
         poyt=0.0
         pozt=0.0
         jet=0.0

         poynit=0.0          ! ab
         poynot=0.0          ! ab

      endif


      call SERV_systime(cpuh)
      s_cpuh=s_cpuh+cpuh-cpug

      cpuinou=cpuinou+s_cpud
      cpucomp=cpucomp+s_cpuh-s_cpud


      end subroutine OUT_poyc
