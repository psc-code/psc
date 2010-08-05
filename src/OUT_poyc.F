
#include "config.h"

c     ======================================================================
      module out_poyn
c     ======================================================================
      implicit none

      real(kind=8) fluxit,fluxot
      real(kind=8) ent,poxt,poyt,pozt,jet
      real(kind=8) enEXt,enEYt,enEZt
      real(kind=8) enHXt,enHYt,enHZt

#ifdef USE_PML
      real(kind=8) enHXt,enHYt,enHZt            ! ab
      real(kind=8) enDXt,enDYt,enDZt            ! ab
      real(kind=8) poynit, poynot               ! ab
#endif

      end module out_poyn

c     ----------------------------------------------------------------------
      subroutine OUT_poyn_setup
c     ----------------------------------------------------------------------
      use out_poyn

      implicit none
      
      fluxit=0.0d0
      fluxot=0.0d0

      ent=0.0d0
      poxt=0.0d0
      poyt=0.0d0
      pozt=0.0d0
      jet=0.0d0

      enEXt=0.0d0
      enEYt=0.0d0
      enEZt=0.0d0
      enHXt=0.0d0
      enHYt=0.0d0
      enHZt=0.0d0

#ifdef USE_PML
      enHXt=0.0d0   ! ab
      enHYt=0.0d0   ! ab
      enHZt=0.0d0   ! ab
      enDXt=0.0d0   ! ab
      enDYt=0.0d0   ! ab
      enDZt=0.0d0   ! ab
#endif

      end subroutine OUT_poyn_setup

c     ----------------------------------------------------------------------
      subroutine OUT_poyn_serv_write
c     ----------------------------------------------------------------------

      use PIC_variables
      use VLA_variables
      use out_poyn

      write(10) fluxit,fluxot
      write(10) ent,poxt,poyt,pozt,jet
      write(10) enEXt,enEYt,enEZt
      write(10) enHXt,enHYt,enHZt
#ifdef USE_PML
      write(10) enHXt,enHYt,enHZt
#endif

      end subroutine OUT_poyn_serv_write

c     ----------------------------------------------------------------------
      subroutine OUT_poyn_serv_read
c     ----------------------------------------------------------------------

      use PIC_variables
      use VLA_variables
      use out_poyn
      read(10) fluxit,fluxot
      read(10) ent,poxt,poyt,pozt,jet
      read(10) enEXt,enEYt,enEZt
      read(10) enHXt,enHYt,enHZt
#ifdef USE_PML      
      read(10) enHXt,enHYt,enHZt
#endif

      end subroutine OUT_poyn_serv_read

c THIS SUBROUTINE CHECKS THE POYNTING LAW IN 3D.

c     ----------------------------------------------------------------------
      subroutine OUT_poyc
c     ----------------------------------------------------------------------

      use PIC_variables
      use VLA_variables
      use out_poyn

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
         enHXt=enHXt+hx2B-hx2A
         enHYt=enHYt+hy2B-hy2A
         enHZt=enHZt+hz2B-hz2A
#ifdef USE_PML
         enHXt=enHXt+hx2B-hx2A
         enHYt=enHYt+hy2B-hy2A
         enHZt=enHZt+hz2B-hz2A
#endif

#ifdef USE_PML
         ent=ent+ex2B-ex2A+ey2B-ey2A+ez2B-ez2A
     &          +hx2B-hx2A+hy2B-hy2A+hz2B-hz2A
#else
         ent=ent+ex2B-ex2A+ey2B-ey2A+ez2B-ez2A
     &          +hx2B-hx2A+hy2B-hy2A+hz2B-hz2A
#endif

         poxt=poxt+dt*pox
         poyt=poyt+dt*poy
         pozt=pozt+dt*poz

#ifdef USE_PML
         poynit = poynit+dt*poyni          ! added hy ab
         poynot = poynot+dt*poyno          ! added hy ab
#endif

         jet=jet+dx*dy*dz*(p2B)
c         jet=jet+dt*je

      endif


      if ((n.gt.tmnvp).and.(n.lt.tmxvp)) then
         fluxit=fluxit+dt*fluxi
         fluxot=fluxot+dt*fluxo

         enEXt=enEXt+ex2B-ex2A
         enEYt=enEYt+ey2B-ey2A
         enEZt=enEZt+ez2B-ez2A
         enHXt=enHXt+hx2B-hx2A
         enHYt=enHYt+hy2B-hy2A
         enHZt=enHZt+hz2B-hz2A
#ifdef USE_PML
         enHXt=enHXt+hx2B-hx2A
         enHYt=enHYt+hy2B-hy2A
         enHZt=enHZt+hz2B-hz2A
#endif

#ifdef USE_PML
         ent=ent+ex2B-ex2A+ey2B-ey2A+ez2B-ez2A
     &          +hx2B-hx2A+hy2B-hy2A+hz2B-hz2A
#else
         ent=ent+ex2B-ex2A+ey2B-ey2A+ez2B-ez2A
     &          +hx2B-hx2A+hy2B-hy2A+hz2B-hz2A
#endif

         poxt=poxt+dt*pox
         poyt=poyt+dt*poy
         pozt=pozt+dt*poz

#ifdef USE_PML
         poynit = poynit+dt*poyni          ! added hy ab
         poynot = poynot+dt*poyno          ! added hy ab
#endif

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
         enHXt=enHXt+hx2B-hx2A
         enHYt=enHYt+hy2B-hy2A
         enHZt=enHZt+hz2B-hz2A
#ifdef USE_PML
         enHXt=enHXt+hx2B-hx2A
         enHYt=enHYt+hy2B-hy2A
         enHZt=enHZt+hz2B-hz2A
#endif

#ifdef USE_PML
         ent=ent+ex2B-ex2A+ey2B-ey2A+ez2B-ez2A
     &          +hx2B-hx2A+hy2B-hy2A+hz2B-hz2A
#else
         ent=ent+ex2B-ex2A+ey2B-ey2A+ez2B-ez2A
     &          +hx2B-hx2A+hy2B-hy2A+hz2B-hz2A
#endif

         poxt=poxt+dt*pox
         poyt=poyt+dt*poy
         pozt=pozt+dt*poz

#ifdef USE_PML
         poynit = poynit+dt*poyni          ! added hy ab
         poynot = poynot+dt*poyno          ! added hy ab
#endif

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
         write(11) enHXt
         write(11) enHYt
         write(11) enHZt
#ifdef USE_PML
         write(11) enHXt
         write(11) enHYt
         write(11) enHZt
#endif
         write(11) ent
         write(11) poxt
         write(11) poyt
         write(11) pozt
         write(11) jet
#ifdef USE_PML
         write(11) poynit               ! ab
         write(11) poynot               ! ab
#endif

         close(11)

         call SERV_systime(cpud)
         s_cpud=s_cpud+cpud-cpuc

         fluxit=0.0
         fluxot=0.0

         enEXt=0.0
         enEYt=0.0
         enEZt=0.0
         enHXt=0.0
         enHYt=0.0
         enHZt=0.0
#ifdef USE_PML
         enHXt=0.0
         enHYt=0.0
         enHZt=0.0
#endif

         ent=0.0
         poxt=0.0
         poyt=0.0
         pozt=0.0
         jet=0.0

#ifdef USE_PML
         poynit=0.0          ! ab
         poynot=0.0          ! ab
#endif

      endif


      call SERV_systime(cpuh)
      s_cpuh=s_cpuh+cpuh-cpug

      cpuinou=cpuinou+s_cpud
      cpucomp=cpucomp+s_cpuh-s_cpud


      end subroutine OUT_poyc
