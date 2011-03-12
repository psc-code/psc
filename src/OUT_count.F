c PARTICLE OUTPUT!


      subroutine OUT_count(niloc, p_niloc)

      use PIC_variables, only: nprparti, dnprparti, nistep
      use VLA_variables
      use PSC_patch

      implicit none

      integer      :: niloc
      real(kind=8) :: p_niloc(0:*)

      integer :: count_in
      integer :: count_ou

      character*5 label,node

      integer :: pli

      real(kind=4) xq,yq,zq
      real(kind=4) x0,y0,z0,xr,yr,zr,rot
      real(kind=4) pxq,pyq,pzq
      real(kind=4) qnq,mnq,cnq,lnq,wnq

      real(kind=8) :: s_cpub
      real(kind=8) :: s_cpud
      real(kind=8) :: s_cpuh


      s_cpub=0.0d0
      s_cpud=0.0d0
      s_cpuh=0.0d0
      call SERV_systime(cpug)


      x0=5.0*1.0e-6/ld
      y0=5.0*1.0e-6/ld
      z0=5.0*1.0e-6/ld
      rot=45.0*3.141592/180.0


c OUTPUT OF PARTICLES AT SELECTED TIMESTEPS


      if (n.eq.nprparti) then
         nprparti=nprparti+dnprparti

         call SERV_labelgen(mpe,node)
         call SERV_labelgen(n,label)

         call SERV_systime(cpuc)

         open(11,file=trim(data_out)
     &        //'/'//node//'count'//label,
     &        access='sequential',form='unformatted')

         count_in=0
         count_ou=0

         do l=1,niloc,nistep

            xq=p_niloc(11*l+0)
            yq=p_niloc(11*l+1)
            zq=p_niloc(11*l+2)
            qnq=p_niloc(11*l+6)
            pli=p_niloc(11*l+9)

            xr=xq
            yr=cos(rot)*(yq-y0)-sin(rot)*(zq-z0)+y0
            zr=cos(rot)*(zq-z0)+sin(rot)*(yq-y0)+z0

            if (((r3n-1)*dz.le.zr.and.zr.le.(r3x+1)*dz).and.
     &         ((r2n-1)*dy.le.yr.and.yr.le.(r2x+1)*dy).and.
     &         ((r1n-1)*dx.le.xr.and.xr.le.(r1x+1)*dx)) then

               if (pli.ne.0) then
                  count_in=count_in+wnq
                  p_niloc(11*l+9)=0.0
               endif

            else

               if (pli.eq.0) then
                  count_ou=count_ou+wnq
                  p_niloc(11*l+9)=1.0
               endif

            endif

         enddo
              
         write(11) count_in, count_ou, n*dt

         close(11)

         call SERV_systime(cpud)
         s_cpud=s_cpud+cpud-cpuc

      endif


      call SERV_systime(cpuh)
      s_cpuh=s_cpuh+cpuh-cpug

      cpumess=cpumess+s_cpub
      cpuinou=cpuinou+s_cpud
      cpucomp=cpucomp+s_cpuh-s_cpud-s_cpub


      end subroutine OUT_count





