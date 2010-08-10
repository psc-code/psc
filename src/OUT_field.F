
#include "config.h"

c FIELD OUTPUT!

      subroutine OUT_field

      implicit none

#ifdef USE_C_INTERFACE
      call C_out_field(99)
#else
      call OUT_field_1
#endif

      end subroutine OUT_field

c legacy output routine

      subroutine OUT_field_1


      use PIC_variables
      use VLA_variables

      implicit none
      integer :: i1a,i1e,i2a,i2e,i3a,i3e

      character*(5) node,label

      real(kind=8) :: exs,eys,ezs
      real(kind=8) :: bxs,bys,bzs
      real(kind=8) :: nes,nis,nns
      real(kind=8) :: jxis,jyis,jzis
#ifdef USE_PML
      real(kind=8) :: hxs,hys,hzs
      real(kind=8) :: ex2s,ey2s,ez2s
      real(kind=8) :: bx2s,by2s,bz2s
      real(kind=8) :: hx2s,hy2s,hz2s
      real(kind=8) :: jxexis,jyeyis,jzezis
      real(kind=8) :: poyxs,poyys,poyzs
#endif
      real(kind=8) :: x0,y0,z0
      real(kind=8) :: xr,yr,zr,rot

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


c CALCULATION OF TIME AVERAGED FIELDS


      if (n.eq.tmnvf) then

         net=0.0
         nit=0.0
         nnt=0.0
         jxit=0.0
         jyit=0.0
         jzit=0.0

         ext=0.0
         eyt=0.0
         ezt=0.0
         bxt=0.0
         byt=0.0
         bzt=0.0
#ifdef USE_PML
         hxt=0.0
         hyt=0.0
         hzt=0.0

         jxexit=0.0
         jyeyit=0.0
         jzezit=0.0
         poyxt=0.0
         poyyt=0.0
         poyzt=0.0

         ex2t=0.0
         ey2t=0.0
         ez2t=0.0
         bx2t=0.0
         by2t=0.0
         bz2t=0.0
         hx2t=0.0
         hy2t=0.0
         hz2t=0.0
#endif
         do i3=i3mn,i3mx
            do i2=i2mn,i2mx
               do i1=i1mn,i1mx

                  nes=ne(i1,i2,i3)
                  nis=ni(i1,i2,i3)
                  nns=nn(i1,i2,i3)
                  jxis=0.5*(jxi(i1,i2,i3)+jxi(i1-1,i2,i3))
                  jyis=0.5*(jyi(i1,i2,i3)+jyi(i1,i2-1,i3))
                  jzis=0.5*(jzi(i1,i2,i3)+jzi(i1,i2,i3-1))
                  exs=0.5*(ex(i1,i2,i3)+ex(i1-1,i2,i3))
                  eys=0.5*(ey(i1,i2,i3)+ey(i1,i2-1,i3))
                  ezs=0.5*(ez(i1,i2,i3)+ez(i1,i2,i3-1))
                  bxs=0.25*(bx(i1,i2,i3)+bx(i1,i2,i3-1)
     &                      +bx(i1,i2-1,i3)+bx(i1,i2-1,i3-1))
                  bys=0.25*(by(i1,i2,i3)+by(i1,i2,i3-1)
     &                      +by(i1-1,i2,i3)+by(i1-1,i2,i3-1))
                  bzs=0.25*(bz(i1,i2,i3)+bz(i1,i2-1,i3)
     &                      +bz(i1-1,i2,i3)+bz(i1-1,i2-1,i3))
#ifdef USE_PML
                  hxs=0.25*(hx(i1,i2,i3)+hx(i1,i2,i3-1)
     &                      +hx(i1,i2-1,i3)+hx(i1,i2-1,i3-1))
                  hys=0.25*(hy(i1,i2,i3)+hy(i1,i2,i3-1)
     &                      +hy(i1-1,i2,i3)+hy(i1-1,i2,i3-1))
                  hzs=0.25*(hz(i1,i2,i3)+hz(i1,i2-1,i3)
     &                      +hz(i1-1,i2,i3)+hz(i1-1,i2-1,i3))
                  jxexis=exs*jxis
                  jyeyis=eys*jyis
                  jzezis=ezs*jzis
!                  poyxs=eys*bzs-ezs*bys
!                  poyys=ezs*bxs-exs*bzs
!                  poyzs=exs*bys-eys*bxs
                  poyxs=eys*hzs-ezs*hys
                  poyys=ezs*hxs-exs*hzs
                  poyzs=exs*hys-eys*hxs
                  ex2s=exs*exs
                  ey2s=eys*eys
                  ez2s=ezs*ezs
                  bx2s=bxs*bxs
                  by2s=bys*bys
                  bz2s=bzs*bzs
                  hx2s=hxs*hxs
                  hy2s=hys*hys
                  hz2s=hzs*hzs
#endif
                  net(i1,i2,i3)=net(i1,i2,i3)+0.5*nes/np
                  nit(i1,i2,i3)=nit(i1,i2,i3)+0.5*nis/np
                  nnt(i1,i2,i3)=nnt(i1,i2,i3)+0.5*nns/np
                  jxit(i1,i2,i3)=jxit(i1,i2,i3)+0.5*jxis/np
                  jyit(i1,i2,i3)=jyit(i1,i2,i3)+0.5*jyis/np
                  jzit(i1,i2,i3)=jzit(i1,i2,i3)+0.5*jzis/np
                  ext(i1,i2,i3)=ext(i1,i2,i3)+0.5*exs/np
                  eyt(i1,i2,i3)=eyt(i1,i2,i3)+0.5*eys/np
                  ezt(i1,i2,i3)=ezt(i1,i2,i3)+0.5*ezs/np
                  bxt(i1,i2,i3)=bxt(i1,i2,i3)+0.5*bxs/np
                  byt(i1,i2,i3)=byt(i1,i2,i3)+0.5*bys/np
                  bzt(i1,i2,i3)=bzt(i1,i2,i3)+0.5*bzs/np
#ifdef USE_PML
                  hxt(i1,i2,i3)=hxt(i1,i2,i3)+0.5*hxs/np
                  hyt(i1,i2,i3)=hyt(i1,i2,i3)+0.5*hys/np
                  hzt(i1,i2,i3)=hzt(i1,i2,i3)+0.5*hzs/np
                  jxexit(i1,i2,i3)=jxexit(i1,i2,i3)+0.5*jxexis/np
                  jyeyit(i1,i2,i3)=jyeyit(i1,i2,i3)+0.5*jyeyis/np
                  jzezit(i1,i2,i3)=jzezit(i1,i2,i3)+0.5*jzezis/np
                  poyxt(i1,i2,i3)=poyxt(i1,i2,i3)+0.5*poyxs/np
                  poyyt(i1,i2,i3)=poyyt(i1,i2,i3)+0.5*poyys/np
                  poyzt(i1,i2,i3)=poyzt(i1,i2,i3)+0.5*poyzs/np
                  ex2t(i1,i2,i3)=ex2t(i1,i2,i3)+0.5*ex2s/np
                  ey2t(i1,i2,i3)=ey2t(i1,i2,i3)+0.5*ey2s/np
                  ez2t(i1,i2,i3)=ez2t(i1,i2,i3)+0.5*ez2s/np
                  bx2t(i1,i2,i3)=bx2t(i1,i2,i3)+0.5*bx2s/np
                  by2t(i1,i2,i3)=by2t(i1,i2,i3)+0.5*by2s/np
                  bz2t(i1,i2,i3)=bz2t(i1,i2,i3)+0.5*bz2s/np
                  hx2t(i1,i2,i3)=hx2t(i1,i2,i3)+0.5*hx2s/np
                  hy2t(i1,i2,i3)=hy2t(i1,i2,i3)+0.5*hy2s/np
                  hz2t(i1,i2,i3)=hz2t(i1,i2,i3)+0.5*hz2s/np
#endif
               enddo
            enddo
         enddo

      endif


      if ((tmnvf.lt.n).and.(n.lt.tmxvf)) then

         do i3=i3mn,i3mx
            do i2=i2mn,i2mx
               do i1=i1mn,i1mx

                  nes=ne(i1,i2,i3)
                  nis=ni(i1,i2,i3)
                  nns=nn(i1,i2,i3)
                  jxis=0.5*(jxi(i1,i2,i3)+jxi(i1-1,i2,i3))
                  jyis=0.5*(jyi(i1,i2,i3)+jyi(i1,i2-1,i3))
                  jzis=0.5*(jzi(i1,i2,i3)+jzi(i1,i2,i3-1))
                  exs=0.5*(ex(i1,i2,i3)+ex(i1-1,i2,i3))
                  eys=0.5*(ey(i1,i2,i3)+ey(i1,i2-1,i3))
                  ezs=0.5*(ez(i1,i2,i3)+ez(i1,i2,i3-1))
                  bxs=0.25*(bx(i1,i2,i3)+bx(i1,i2,i3-1)
     &                      +bx(i1,i2-1,i3)+bx(i1,i2-1,i3-1))
                  bys=0.25*(by(i1,i2,i3)+by(i1,i2,i3-1)
     &                      +by(i1-1,i2,i3)+by(i1-1,i2,i3-1))
                  bzs=0.25*(bz(i1,i2,i3)+bz(i1,i2-1,i3)
     &                      +bz(i1-1,i2,i3)+bz(i1-1,i2-1,i3))
#ifdef USE_PML
                  hxs=0.25*(hx(i1,i2,i3)+hx(i1,i2,i3-1)
     &                      +hx(i1,i2-1,i3)+hx(i1,i2-1,i3-1))
                  hys=0.25*(hy(i1,i2,i3)+hy(i1,i2,i3-1)
     &                      +hy(i1-1,i2,i3)+hy(i1-1,i2,i3-1))
                  hzs=0.25*(hz(i1,i2,i3)+hz(i1,i2-1,i3)
     &                      +hz(i1-1,i2,i3)+hz(i1-1,i2-1,i3))
                  jxexis=exs*jxis
                  jyeyis=eys*jyis
                  jzezis=ezs*jzis
!                  poyxs=eys*bzs-ezs*bys
!                  poyys=ezs*bxs-exs*bzs
!                  poyzs=exs*bys-eys*bxs
                  poyxs=eys*hzs-ezs*hys
                  poyys=ezs*hxs-exs*hzs
                  poyzs=exs*hys-eys*hxs
                  ex2s=exs*exs
                  ey2s=eys*eys
                  ez2s=ezs*ezs
                  bx2s=bxs*bxs
                  by2s=bys*bys
                  bz2s=bzs*bzs
                  hx2s=hxs*hxs
                  hy2s=hys*hys
                  hz2s=hzs*hzs
#endif
                  net(i1,i2,i3)=net(i1,i2,i3)+nes/np
                  nit(i1,i2,i3)=nit(i1,i2,i3)+nis/np
                  nnt(i1,i2,i3)=nnt(i1,i2,i3)+nns/np
                  jxit(i1,i2,i3)=jxit(i1,i2,i3)+jxis/np
                  jyit(i1,i2,i3)=jyit(i1,i2,i3)+jyis/np
                  jzit(i1,i2,i3)=jzit(i1,i2,i3)+jzis/np
                  ext(i1,i2,i3)=ext(i1,i2,i3)+exs/np
                  eyt(i1,i2,i3)=eyt(i1,i2,i3)+eys/np
                  ezt(i1,i2,i3)=ezt(i1,i2,i3)+ezs/np
                  bxt(i1,i2,i3)=bxt(i1,i2,i3)+bxs/np
                  byt(i1,i2,i3)=byt(i1,i2,i3)+bys/np
                  bzt(i1,i2,i3)=bzt(i1,i2,i3)+bzs/np
#ifdef USE_PML
                  hxt(i1,i2,i3)=hxt(i1,i2,i3)+hxs/np
                  hyt(i1,i2,i3)=hyt(i1,i2,i3)+hys/np
                  hzt(i1,i2,i3)=hzt(i1,i2,i3)+hzs/np
                  jxexit(i1,i2,i3)=jxexit(i1,i2,i3)+jxexis/np
                  jyeyit(i1,i2,i3)=jyeyit(i1,i2,i3)+jyeyis/np
                  jzezit(i1,i2,i3)=jzezit(i1,i2,i3)+jzezis/np
                  poyxt(i1,i2,i3)=poyxt(i1,i2,i3)+poyxs/np
                  poyyt(i1,i2,i3)=poyyt(i1,i2,i3)+poyys/np
                  poyzt(i1,i2,i3)=poyzt(i1,i2,i3)+poyzs/np
                  ex2t(i1,i2,i3)=ex2t(i1,i2,i3)+ex2s/np
                  ey2t(i1,i2,i3)=ey2t(i1,i2,i3)+ey2s/np
                  ez2t(i1,i2,i3)=ez2t(i1,i2,i3)+ez2s/np
                  bx2t(i1,i2,i3)=bx2t(i1,i2,i3)+bx2s/np
                  by2t(i1,i2,i3)=by2t(i1,i2,i3)+by2s/np
                  bz2t(i1,i2,i3)=bz2t(i1,i2,i3)+bz2s/np
                  hx2t(i1,i2,i3)=hx2t(i1,i2,i3)+hx2s/np
                  hy2t(i1,i2,i3)=hy2t(i1,i2,i3)+hy2s/np
                  hz2t(i1,i2,i3)=hz2t(i1,i2,i3)+hz2s/np
#endif
               enddo
            enddo
         enddo

      endif


      if (n.eq.tmxvf) then
         tmnvf=n+1
         tmxvf=n+np

         do i3=i3mn,i3mx
            do i2=i2mn,i2mx
               do i1=i1mn,i1mx

                  nes=ne(i1,i2,i3)
                  nis=ni(i1,i2,i3)
                  nns=nn(i1,i2,i3)
                  jxis=0.5*(jxi(i1,i2,i3)+jxi(i1-1,i2,i3))
                  jyis=0.5*(jyi(i1,i2,i3)+jyi(i1,i2-1,i3))
                  jzis=0.5*(jzi(i1,i2,i3)+jzi(i1,i2,i3-1))
                  exs=0.5*(ex(i1,i2,i3)+ex(i1-1,i2,i3))
                  eys=0.5*(ey(i1,i2,i3)+ey(i1,i2-1,i3))
                  ezs=0.5*(ez(i1,i2,i3)+ez(i1,i2,i3-1))
                  bxs=0.25*(bx(i1,i2,i3)+bx(i1,i2,i3-1)
     &                      +bx(i1,i2-1,i3)+bx(i1,i2-1,i3-1))
                  bys=0.25*(by(i1,i2,i3)+by(i1,i2,i3-1)
     &                      +by(i1-1,i2,i3)+by(i1-1,i2,i3-1))
                  bzs=0.25*(bz(i1,i2,i3)+bz(i1,i2-1,i3)
     &                      +bz(i1-1,i2,i3)+bz(i1-1,i2-1,i3))
#ifdef USE_PML
                  hxs=0.25*(hx(i1,i2,i3)+hx(i1,i2,i3-1)
     &                      +hx(i1,i2-1,i3)+hx(i1,i2-1,i3-1))
                  hys=0.25*(hy(i1,i2,i3)+hy(i1,i2,i3-1)
     &                      +hy(i1-1,i2,i3)+hy(i1-1,i2,i3-1))
                  hzs=0.25*(hz(i1,i2,i3)+hz(i1,i2-1,i3)
     &                      +hz(i1-1,i2,i3)+hz(i1-1,i2-1,i3))
                  jxexis=exs*jxis
                  jyeyis=eys*jyis
                  jzezis=ezs*jzis
!                  poyxs=eys*bzs-ezs*bys
!                  poyys=ezs*bxs-exs*bzs
!                  poyzs=exs*bys-eys*bxs
                  poyxs=eys*hzs-ezs*hys
                  poyys=ezs*hxs-exs*hzs
                  poyzs=exs*hys-eys*hxs
                  ex2s=exs*exs
                  ey2s=eys*eys
                  ez2s=ezs*ezs
                  bx2s=bxs*bxs
                  by2s=bys*bys
                  bz2s=bzs*bzs
                  hx2s=hxs*hxs
                  hy2s=hys*hys
                  hz2s=hzs*hzs
#endif
                  net(i1,i2,i3)=net(i1,i2,i3)+0.5*nes/np
                  nit(i1,i2,i3)=nit(i1,i2,i3)+0.5*nis/np
                  nnt(i1,i2,i3)=nnt(i1,i2,i3)+0.5*nns/np
                  jxit(i1,i2,i3)=jxit(i1,i2,i3)+0.5*jxis/np
                  jyit(i1,i2,i3)=jyit(i1,i2,i3)+0.5*jyis/np
                  jzit(i1,i2,i3)=jzit(i1,i2,i3)+0.5*jzis/np
                  ext(i1,i2,i3)=ext(i1,i2,i3)+0.5*exs/np
                  eyt(i1,i2,i3)=eyt(i1,i2,i3)+0.5*eys/np
                  ezt(i1,i2,i3)=ezt(i1,i2,i3)+0.5*ezs/np
                  bxt(i1,i2,i3)=bxt(i1,i2,i3)+0.5*bxs/np
                  byt(i1,i2,i3)=byt(i1,i2,i3)+0.5*bys/np
                  bzt(i1,i2,i3)=bzt(i1,i2,i3)+0.5*bzs/np
#ifdef USE_PML
                  hxt(i1,i2,i3)=hxt(i1,i2,i3)+0.5*hxs/np
                  hyt(i1,i2,i3)=hyt(i1,i2,i3)+0.5*hys/np
                  hzt(i1,i2,i3)=hzt(i1,i2,i3)+0.5*hzs/np
                  jxexit(i1,i2,i3)=jxexit(i1,i2,i3)+0.5*jxexis/np
                  jyeyit(i1,i2,i3)=jyeyit(i1,i2,i3)+0.5*jyeyis/np
                  jzezit(i1,i2,i3)=jzezit(i1,i2,i3)+0.5*jzezis/np
                  poyxt(i1,i2,i3)=poyxt(i1,i2,i3)+0.5*poyxs/np
                  poyyt(i1,i2,i3)=poyyt(i1,i2,i3)+0.5*poyys/np
                  poyzt(i1,i2,i3)=poyzt(i1,i2,i3)+0.5*poyzs/np
                  ex2t(i1,i2,i3)=ex2t(i1,i2,i3)+0.5*ex2s/np
                  ey2t(i1,i2,i3)=ey2t(i1,i2,i3)+0.5*ey2s/np
                  ez2t(i1,i2,i3)=ez2t(i1,i2,i3)+0.5*ez2s/np
                  bx2t(i1,i2,i3)=bx2t(i1,i2,i3)+0.5*bx2s/np
                  by2t(i1,i2,i3)=by2t(i1,i2,i3)+0.5*by2s/np
                  bz2t(i1,i2,i3)=bz2t(i1,i2,i3)+0.5*bz2s/np
                  hx2t(i1,i2,i3)=hx2t(i1,i2,i3)+0.5*hx2s/np
                  hy2t(i1,i2,i3)=hy2t(i1,i2,i3)+0.5*hy2s/np
                  hz2t(i1,i2,i3)=hz2t(i1,i2,i3)+0.5*hz2s/np
#endif
               enddo
            enddo
         enddo

         call SERV_labelgen(mpe,node)
         call SERV_labelgen(n,label)

         call SERV_systime(cpuc)

         open(11,file=trim(data_out)//'/'//node//'tfield'//label,
     &        access='sequential',form='unformatted')

         write(11) i1n,i1x,i2n,i2x,i3n,i3x
         write(11) np
         write(11) nmax,xnpe,ynpe
         write(11) nprf,dnprf
         write(11) dt,dx,dy,dz,shift_z
         write(11) i0,e0,b0,j0,rho0,n0
         write(11) wl,wp,ld
         write(11) alpha,beta,eta
         write(11) vos,vt,cc,qq,mm,tt,eps0
         write(11) i1mn,i1mx,i2mn,i2mx,i3mn,i3mx,shift_z
         write(11) r1n,r1x,r2n,r2x,r3n,r3x
         write(11) sngl(x0),sngl(y0),sngl(z0),sngl(rot)

         do i3=i3mn,i3mx
            do i2=i2mn,i2mx
               do i1=i1mn,i1mx

                  xr=i1*dx
                  yr=cos(rot)*(i2*dy-y0)-sin(rot)*(i3*dz-z0)+y0
                  zr=cos(rot)*(i3*dz-z0)+sin(rot)*(i2*dy-y0)+z0

                  if (((r3n-1)*dz.le.zr.and.zr.le.(r3x+1)*dz).and.
     &                ((r2n-1)*dy.le.yr.and.yr.le.(r2x+1)*dy).and.
     &                ((r1n-1)*dx.le.xr.and.xr.le.(r1x+1)*dx)) then

                     nes=net(i1,i2,i3)
                     nis=nit(i1,i2,i3)
                     nns=nnt(i1,i2,i3)
                     jxis=jxit(i1,i2,i3)
                     jyis=jyit(i1,i2,i3)
                     jzis=jzit(i1,i2,i3)
                     exs=ext(i1,i2,i3)
                     eys=eyt(i1,i2,i3)
                     ezs=ezt(i1,i2,i3)
                     bxs=bxt(i1,i2,i3)
                     bys=byt(i1,i2,i3)
                     bzs=bzt(i1,i2,i3)
                     write(11) sngl(nes)
                     write(11) sngl(nis)
                     write(11) sngl(nns)
                     write(11) sngl(jxis)
                     write(11) sngl(jyis)
                     write(11) sngl(jzis)
                     write(11) sngl(exs)
                     write(11) sngl(eys)
                     write(11) sngl(ezs)
                     write(11) sngl(bxs)
                     write(11) sngl(bys)
                     write(11) sngl(bzs)
#ifdef USE_PML
                     hxs=hxt(i1,i2,i3)
                     hys=hyt(i1,i2,i3)
                     hzs=hzt(i1,i2,i3)
                     jxexis=jxexit(i1,i2,i3)
                     jyeyis=jyeyit(i1,i2,i3)
                     jzezis=jzezit(i1,i2,i3)
                     poyxs=poyxt(i1,i2,i3)
                     poyys=poyyt(i1,i2,i3)
                     poyzs=poyzt(i1,i2,i3)
                     ex2s=ex2t(i1,i2,i3)
                     ey2s=ey2t(i1,i2,i3)
                     ez2s=ez2t(i1,i2,i3)
                     bx2s=bx2t(i1,i2,i3)
                     by2s=by2t(i1,i2,i3)
                     bz2s=bz2t(i1,i2,i3)
                     hx2s=hx2t(i1,i2,i3)
                     hy2s=hy2t(i1,i2,i3)
                     hz2s=hz2t(i1,i2,i3)
                     write(11) sngl(hxs)
                     write(11) sngl(hys)
                     write(11) sngl(hzs)
                     write(11) sngl(jxexis)
                     write(11) sngl(jyeyis)
                     write(11) sngl(jzezis)
                     write(11) sngl(poyxs)
                     write(11) sngl(poyys)
                     write(11) sngl(poyzs)
                     write(11) sngl(ex2s)
                     write(11) sngl(ey2s)
                     write(11) sngl(ez2s)
                     write(11) sngl(bx2s)
                     write(11) sngl(by2s)
                     write(11) sngl(bz2s)
                     write(11) sngl(hx2s)
                     write(11) sngl(hy2s)
                     write(11) sngl(hz2s)
#endif

                  endif

               enddo
            enddo
         enddo

         close(11)

         call SERV_systime(cpud)
         s_cpud=s_cpud+cpud-cpuc

      endif


c OUTPUT AFTER PREDEFINED TIME INTERVALS


      if (n.eq.nprf) then
         nprf=nprf+dnprf

         call SERV_labelgen(mpe,node)
         call SERV_labelgen(n,label)

         call SERV_systime(cpuc)

         open(11,file=trim(data_out)//'/'//node//'pfield'//label,
     &         access='sequential',form='unformatted')

         write(11) i1n,i1x,i2n,i2x,i3n,i3x
         write(11) np
         write(11) nmax,xnpe,ynpe
         write(11) nprf,dnprf
         write(11) dt,dx,dy,dz,shift_z
         write(11) i0,e0,b0,j0,rho0,n0
         write(11) wl,wp,ld
         write(11) alpha,beta,eta
         write(11) vos,vt,cc,qq,mm,tt,eps0
         write(11) i1mn,i1mx,i2mn,i2mx,i3mn,i3mx,shift_z
         write(11) r1n,r1x,r2n,r2x,r3n,r3x
         write(11) sngl(x0),sngl(y0),sngl(z0),sngl(rot)

         do i3=i3mn,i3mx
            do i2=i2mn,i2mx
               do i1=i1mn,i1mx

                  xr=i1*dx
                  yr=dcos(rot)*(i2*dy-y0)-sin(rot)*(i3*dz-z0)+y0
                  zr=dcos(rot)*(i3*dz-z0)+sin(rot)*(i2*dy-y0)+z0

                  if (((r3n-1)*dz.le.zr.and.zr.le.(r3x+1)*dz).and.
     &                ((r2n-1)*dy.le.yr.and.yr.le.(r2x+1)*dy).and.
     &                ((r1n-1)*dx.le.xr.and.xr.le.(r1x+1)*dx)) then

                     nes=ne(i1,i2,i3)
                     nis=ni(i1,i2,i3)
                     nns=nn(i1,i2,i3)
                     jxis=0.5*(jxi(i1,i2,i3)+jxi(i1-1,i2,i3))
                     jyis=0.5*(jyi(i1,i2,i3)+jyi(i1,i2-1,i3))
                     jzis=0.5*(jzi(i1,i2,i3)+jzi(i1,i2,i3-1))
                     exs=0.5*(ex(i1,i2,i3)+ex(i1-1,i2,i3))
                     eys=0.5*(ey(i1,i2,i3)+ey(i1,i2-1,i3))
                     ezs=0.5*(ez(i1,i2,i3)+ez(i1,i2,i3-1))
                     bxs=0.25*(bx(i1,i2,i3)+bx(i1,i2,i3-1)
     &                         +bx(i1,i2-1,i3)+bx(i1,i2-1,i3-1))
                     bys=0.25*(by(i1,i2,i3)+by(i1,i2,i3-1)
     &                         +by(i1-1,i2,i3)+by(i1-1,i2,i3-1))
                     bzs=0.25*(bz(i1,i2,i3)+bz(i1,i2-1,i3)
     &                         +bz(i1-1,i2,i3)+bz(i1-1,i2-1,i3))
                     write(11) sngl(nes)
                     write(11) sngl(nis)
                     write(11) sngl(nns)
                     write(11) sngl(jxis)
                     write(11) sngl(jyis)
                     write(11) sngl(jzis)
                     write(11) sngl(exs)
                     write(11) sngl(eys)
                     write(11) sngl(ezs)
                     write(11) sngl(bxs)
                     write(11) sngl(bys)
                     write(11) sngl(bzs)

#ifdef USE_PML
                     hxs=0.25*(hx(i1,i2,i3)+hx(i1,i2,i3-1)
     &                         +hx(i1,i2-1,i3)+hx(i1,i2-1,i3-1))
                     hys=0.25*(hy(i1,i2,i3)+hy(i1,i2,i3-1)
     &                         +hy(i1-1,i2,i3)+hy(i1-1,i2,i3-1))
                     hzs=0.25*(hz(i1,i2,i3)+hz(i1,i2-1,i3)
     &                         +hz(i1-1,i2,i3)+hz(i1-1,i2-1,i3))
                     jxexis=exs*jxis
                     jyeyis=eys*jyis
                     jzezis=ezs*jzis
!                     poyxs=eys*bzs-ezs*bys
!                     poyys=ezs*bxs-exs*bzs
!                     poyzs=exs*bys-eys*bxs
                     poyxs=eys*hzs-ezs*hys
                     poyys=ezs*hxs-exs*hzs
                     poyzs=exs*hys-eys*hxs
                     ex2s=exs*exs
                     ey2s=eys*eys
                     ez2s=ezs*ezs
                     bx2s=bxs*bxs
                     by2s=bys*bys
                     bz2s=bzs*bzs
                     hx2s=hxs*hxs
                     hy2s=hys*hys
                     hz2s=hzs*hzs
                     write(11) sngl(hxs)
                     write(11) sngl(hys)
                     write(11) sngl(hzs)
                     write(11) sngl(jxexis)
                     write(11) sngl(jyeyis)
                     write(11) sngl(jzezis)
                     write(11) sngl(poyxs)
                     write(11) sngl(poyys)
                     write(11) sngl(poyzs)
                     write(11) sngl(ex2s)
                     write(11) sngl(ey2s)
                     write(11) sngl(ez2s)
                     write(11) sngl(bx2s)
                     write(11) sngl(by2s)
                     write(11) sngl(bz2s)
                     write(11) sngl(hx2s)
                     write(11) sngl(hy2s)
                     write(11) sngl(hz2s)
#endif

                  endif
               enddo
            enddo
         enddo
              
         close(11)

         call SERV_systime(cpud)
         s_cpud=s_cpud+cpud-cpuc

      endif 


      call SERV_systime(cpuh)
      s_cpuh=s_cpuh+cpuh-cpug


      cpumess=cpumess+s_cpub
      cpuinou=cpuinou+s_cpud
      cpucomp=cpucomp+s_cpuh-s_cpud-s_cpub


      end subroutine OUT_field_1

