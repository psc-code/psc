c PARTICLE OUTPUT!


      subroutine OUT_count(niloc, p_niloc)

      use OUT_params
      use OUT_particles
      use VLA_variables, only: data_out, n
      use PSC_domain, only: mpe
      use PSC_params
      use PSC_patch

      implicit none

      integer      :: niloc
      real(kind=8) :: p_niloc(0:*)

      integer :: count_in
      integer :: count_ou

      character*5 label,node

      integer :: l,pli

      real(kind=4) xq,yq,zq
      real(kind=4) x0,y0,z0,xr,yr,zr,rot
      real(kind=4) pxq,pyq,pzq
      real(kind=4) qnq,mnq,cnq,lnq,wnq


      x0=5.0*1.0e-6/ld
      y0=5.0*1.0e-6/ld
      z0=5.0*1.0e-6/ld
      rot=45.0*3.141592/180.0


c OUTPUT OF PARTICLES AT SELECTED TIMESTEPS


      if (n.eq.nprparti) then
         nprparti=nprparti+dnprparti

         call SERV_labelgen(mpe,node)
         call SERV_labelgen(n,label)

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

      endif

      end subroutine OUT_count





