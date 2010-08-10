c PARTICLE OUTPUT!


      subroutine OUT_part

      use PIC_variables
      use VLA_variables

      implicit none

      character*5 label,node

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
      rot=0.0*3.141592/180.0


c OUTPUT OF PARTICLES AT SELECTED TIMESTEPS


      if (n.eq.nprparti) then
         nprparti=nprparti+dnprparti

         call SERV_labelgen(mpe,node)
         call SERV_labelgen(n,label)

         call SERV_systime(cpuc)

         open(11,file=trim(data_out)
     &        //'/'//node//'electron'//label,
     &        access='sequential',form='unformatted')
         open(12,file=trim(data_out)
     &        //'/'//node//'ion'//label,
     &        access='sequential',form='unformatted')
         open(13,file=trim(data_out)
     &        //'/'//node//'atom'//label,
     &        access='sequential',form='unformatted')


         write(11) cori,shift_z
         write(11) i1n,i1x,i2n,i2x,i3n,i3x
         write(11) np
         write(11) nmax,xnpe,ynpe
         write(11) nprparti,dnprparti,nistep
         write(11) dt,dx,dy,dz
         write(11) i0,e0,b0,j0,rho0,n0,cori,shift_z
         write(11) wl,wp,ld
         write(11) alpha,beta,eta
         write(11) vos,vt,cc,qq,mm,tt,eps0

         write(12) cori,shift_z
         write(12) i1n,i1x,i2n,i2x,i3n,i3x
         write(12) np
         write(12) nmax,xnpe,ynpe
         write(12) nprparti,dnprparti,nistep
         write(12) dt,dx,dy,dz
         write(12) i0,e0,b0,j0,rho0,n0,cori,shift_z
         write(12) wl,wp,ld
         write(12) alpha,beta,eta
         write(12) vos,vt,cc,qq,mm,tt,eps0

         write(13) cori,shift_z
         write(13) i1n,i1x,i2n,i2x,i3n,i3x
         write(13) np
         write(13) nmax,xnpe,ynpe
         write(13) nprparti,dnprparti,nistep
         write(13) dt,dx,dy,dz
         write(13) i0,e0,b0,j0,rho0,n0,cori,shift_z
         write(13) wl,wp,ld
         write(13) alpha,beta,eta
         write(13) vos,vt,cc,qq,mm,tt,eps0


         do l=1,niloc,nistep

            xq=p_niloc(11*l+0)
            yq=p_niloc(11*l+1)
            zq=p_niloc(11*l+2)
            qnq=p_niloc(11*l+6)
            pli=p_niloc(11*l+9)

            xr=xq
            yr=cos(rot)*(yq-y0)-sin(rot)*(zq-z0)+y0
            zr=cos(rot)*(zq-z0)+sin(rot)*(yq-y0)+z0

            if (plin.le.pli.and.pli.le.plix) then

               if (((r3n-1)*dz.le.zr.and.zr.le.(r3x+1)*dz).and.
     &             ((r2n-1)*dy.le.yr.and.yr.le.(r2x+1)*dy).and.
     &             ((r1n-1)*dx.le.xr.and.xr.le.(r1x+1)*dx)) then

                  if (qnq<0.0) then
                     xq=p_niloc(11*l)
                     yq=p_niloc(11*l+1)
                     zq=p_niloc(11*l+2)
                     pxq=p_niloc(11*l+3)
                     pyq=p_niloc(11*l+4)
                     pzq=p_niloc(11*l+5)
                     qnq=p_niloc(11*l+6)
                     mnq=p_niloc(11*l+7)
                     cnq=p_niloc(11*l+8)
                     lnq=p_niloc(11*l+9)
                     wnq=p_niloc(11*l+10)
                     write(11) xq,yq,zq,pxq,pyq,pzq,qnq,mnq,cnq,lnq,wnq
                  endif

                  if (qnq>0.0) then
                     xq=p_niloc(11*l)
                     yq=p_niloc(11*l+1)
                     zq=p_niloc(11*l+2)
                     pxq=p_niloc(11*l+3)
                     pyq=p_niloc(11*l+4)
                     pzq=p_niloc(11*l+5)
                     qnq=p_niloc(11*l+6)
                     mnq=p_niloc(11*l+7)
                     cnq=p_niloc(11*l+8)
                     lnq=p_niloc(11*l+9)
                     wnq=p_niloc(11*l+10)
                     write(12) xq,yq,zq,pxq,pyq,pzq,qnq,mnq,cnq,lnq,wnq
                  endif

                  if (qnq==0.0) then
                     xq=p_niloc(11*l)
                     yq=p_niloc(11*l+1)
                     zq=p_niloc(11*l+2)
                     pxq=p_niloc(11*l+3)
                     pyq=p_niloc(11*l+4)
                     pzq=p_niloc(11*l+5)
                     qnq=p_niloc(11*l+6)
                     mnq=p_niloc(11*l+7)
                     cnq=p_niloc(11*l+8)
                     lnq=p_niloc(11*l+9)
                     wnq=p_niloc(11*l+10)
                     write(13) xq,yq,zq,pxq,pyq,pzq,qnq,mnq,cnq,lnq,wnq
                  endif

               endif

            endif
         enddo

         call SERV_systime(cpud)
         s_cpud=s_cpud+cpud-cpuc

         close(11)
         close(12)
         close(13)

      endif

      call SERV_systime(cpuh)
      s_cpuh=s_cpuh+cpuh-cpug

      cpumess=cpumess+s_cpub
      cpuinou=cpuinou+s_cpud
      cpucomp=cpucomp+s_cpuh-s_cpud-s_cpub


      end subroutine OUT_part





