
c     ======================================================================
      module OUT_poyn
c     ======================================================================
      implicit none

c tmnvp: time step counter for starting time averaging of poynting flux
c tmxvp: time step counter for ending time averaging of poynting flux

      integer :: tmnvpo,tmxvpo,tmnvp,tmxvp

      real(kind=8) fluxit,fluxot
      real(kind=8) ent,poxt,poyt,pozt,jet
      real(kind=8) enEXt,enEYt,enEZt
      real(kind=8) enHXt,enHYt,enHZt

      real(kind=8) enBXt,enBYt,enBZt            ! ab
      real(kind=8) enDXt,enDYt,enDZt            ! ab
      real(kind=8) poynit, poynot               ! ab

      real(kind=8) fluxi,fluxo
      real(kind=8) ex2A,ey2A,ez2A,hx2A,hy2A,hz2A
      real(kind=8) ex2B,ey2B,ez2B,hx2B,hy2B,hz2B
      real(kind=8) pox,poy,poz,je
      real(kind=8) p2A,p2B

      end module out_poyn

c     ----------------------------------------------------------------------
      subroutine OUT_POYN_alloc
c     ----------------------------------------------------------------------
      use OUT_params
      use OUT_poyn

      implicit none
      
c tmnvp: starting time averaging of poynting flux
c tmxvp: ending time averaging of poynting flux

      tmnvp=0*nnp+1         
      tmxvp=0*nnp+np

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

      enBXt=0.0d0   ! ab
      enBYt=0.0d0   ! ab
      enBZt=0.0d0   ! ab
      enDXt=0.0d0   ! ab
      enDYt=0.0d0   ! ab
      enDZt=0.0d0   ! ab

      end subroutine OUT_POYN_alloc

c     ----------------------------------------------------------------------
      subroutine OUT_POYN_serv_write
c     ----------------------------------------------------------------------

      use OUT_POYN

      write(10) fluxit,fluxot
      write(10) ent,poxt,poyt,pozt,jet
      write(10) enEXt,enEYt,enEZt
      write(10) enHXt,enHYt,enHZt
      write(10) enBXt,enBYt,enBZt

      end subroutine OUT_POYN_serv_write

c     ----------------------------------------------------------------------
      subroutine OUT_POYN_serv_read
c     ----------------------------------------------------------------------

      use OUT_POYN

      read(10) fluxit,fluxot
      read(10) ent,poxt,poyt,pozt,jet
      read(10) enEXt,enEYt,enEZt
      read(10) enHXt,enHYt,enHZt
      read(10) enBXt,enBYt,enBZt

      end subroutine OUT_POYN_serv_read

c THIS SUBROUTINE CHECKS THE POYNTING LAW IN 3D.

c     ----------------------------------------------------------------------
      subroutine OUT_poyc
c     ----------------------------------------------------------------------

      use VLA_variables
      use OUT_params
      use OUT_POYN
      use PSC_domain
      use PSC_patch
      use COEFF_pml

      implicit none

      character*5 label,node


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
         enBXt=enBXt+bx2B-bx2A
         enBYt=enBYt+by2B-by2A
         enBZt=enBZt+bz2B-bz2A

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
         enHXt=enHXt+hx2B-hx2A
         enHYt=enHYt+hy2B-hy2A
         enHZt=enHZt+hz2B-hz2A
         enBXt=enBXt+bx2B-bx2A
         enBYt=enBYt+by2B-by2A
         enBZt=enBZt+bz2B-bz2A

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
         enHXt=enHXt+hx2B-hx2A
         enHYt=enHYt+hy2B-hy2A
         enHZt=enHZt+hz2B-hz2A
         enBXt=enBXt+bx2B-bx2A
         enBYt=enBYt+by2B-by2A
         enBZt=enBZt+bz2B-bz2A

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
         write(11) enBXt
         write(11) enBYt
         write(11) enBZt
         write(11) ent
         write(11) poxt
         write(11) poyt
         write(11) pozt
         write(11) jet
         write(11) poynit               ! ab
         write(11) poynot               ! ab

         close(11)

         fluxit=0.0
         fluxot=0.0

         enEXt=0.0
         enEYt=0.0
         enEZt=0.0
         enHXt=0.0
         enHYt=0.0
         enHZt=0.0
         enBXt=0.0
         enBYt=0.0
         enBZt=0.0

         ent=0.0
         poxt=0.0
         poyt=0.0
         pozt=0.0
         jet=0.0

         poynit=0.0          ! ab
         poynot=0.0          ! ab

      endif


      end subroutine OUT_poyc
