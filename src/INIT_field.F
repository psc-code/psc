
#include "config.h"

c THIS SUBROUTINE ALLOCATES THE ELECTROMAGNETIC FIELDS.

c     ----------------------------------------------------------------------
      subroutine ALLOC_field
c     ----------------------------------------------------------------------

      use PIC_variables
      use VLA_variables

      implicit none

      allocate(ex(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))
      allocate(ey(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))
      allocate(ez(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))
      allocate(bx(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))
      allocate(by(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))
      allocate(bz(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))
      allocate(jxi(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &             i3mn-rd3:i3mx+rd3))
      allocate(jyi(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &             i3mn-rd3:i3mx+rd3))
      allocate(jzi(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &             i3mn-rd3:i3mx+rd3))
      allocate(ne(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))
      allocate(ni(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))
      allocate(nn(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))

      allocate(ext(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(eyt(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(ezt(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(bxt(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(byt(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(bzt(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(net(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(nit(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(nnt(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(jxit(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(jyit(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(jzit(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      
      allocate(cpu(1:5))
      allocate(cpu_ary(1:5,0:npe-1))

! array allocation
! size of pml should be equal to size of em field arrays

! x

      allocate(kappax(i1mn-rd1:i1mx+rd1))
      allocate(sigmax(i1mn-rd1:i1mx+rd1))
      
      allocate(cxp(i1mn-rd1:i1mx+rd1))
      allocate(cxm(i1mn-rd1:i1mx+rd1))
     
      allocate(fbx(i1mn-rd1:i1mx+rd1))
      allocate(fcx(i1mn-rd1:i1mx+rd1))
      allocate(fdx(i1mn-rd1:i1mx+rd1))
      allocate(fex(i1mn-rd1:i1mx+rd1))

      allocate(bxp(i1mn-rd1:i1mx+rd1))
      allocate(bxm(i1mn-rd1:i1mx+rd1))
      
      allocate(gbx(i1mn-rd1:i1mx+rd1))
      allocate(gcx(i1mn-rd1:i1mx+rd1))
      allocate(gdx(i1mn-rd1:i1mx+rd1))
      allocate(gex(i1mn-rd1:i1mx+rd1))

! y
      
      allocate(kappay(i2mn-rd2:i2mx+rd2))
      allocate(sigmay(i2mn-rd2:i2mx+rd2))
      
      allocate(cyp(i2mn-rd2:i2mx+rd2))
      allocate(cym(i2mn-rd2:i2mx+rd2))
      
      allocate(fby(i2mn-rd2:i2mx+rd2))
      allocate(fcy(i2mn-rd2:i2mx+rd2))
      allocate(fdy(i2mn-rd2:i2mx+rd2))
      allocate(fey(i2mn-rd2:i2mx+rd2))

      allocate(byp(i2mn-rd2:i2mx+rd2))
      allocate(bym(i2mn-rd2:i2mx+rd2))
      
      allocate(gby(i2mn-rd2:i2mx+rd2))
      allocate(gcy(i2mn-rd2:i2mx+rd2))
      allocate(gdy(i2mn-rd2:i2mx+rd2))
      allocate(gey(i2mn-rd2:i2mx+rd2))

! z
      
      allocate(kappaz(i3mn-rd3:i3mx+rd3))
      allocate(sigmaz(i3mn-rd3:i3mx+rd3))
      
      allocate(czp(i3mn-rd3:i3mx+rd3))
      allocate(czm(i3mn-rd3:i3mx+rd3))
      
      allocate(fbz(i3mn-rd3:i3mx+rd3))
      allocate(fcz(i3mn-rd3:i3mx+rd3))
      allocate(fdz(i3mn-rd3:i3mx+rd3))
      allocate(fez(i3mn-rd3:i3mx+rd3))

      allocate(bzp(i3mn-rd3:i3mx+rd3))
      allocate(bzm(i3mn-rd3:i3mx+rd3))
      
      allocate(gbz(i3mn-rd3:i3mx+rd3))
      allocate(gcz(i3mn-rd3:i3mx+rd3))
      allocate(gdz(i3mn-rd3:i3mx+rd3))
      allocate(gez(i3mn-rd3:i3mx+rd3))

      ex=0.0d0
      ey=0.0d0
      ez=0.0d0
      bx=0.0d0
      by=0.0d0
      bz=0.0d0
      ne=0.0d0
      ni=0.0d0
      nn=0.0d0
      jxi=0.0d0
      jyi=0.0d0
      jzi=0.0d0


#ifdef USE_PML
      call PML_alloc_fields
#endif

#ifdef USE_C_INTERFACE
c     need to let C code know about the newly allocated fields
      call C_init_vars(dt, dx, dy, dz,
     &     i1mn, i2mn, i3mn, i1mx, i2mx, i3mx, rd1, rd2, rd3,
     &     i1n, i2n, i3n, i1x, i2x, i3x,
     &     ne, ni, nn, jxi, jyi, jzi, ex, ey, ez, bx, by, bz,
     &     cori, alpha, beta, eta, wl, wp,
     &     n,
     &     99)
#endif

      end subroutine ALLOC_field


c THIS SUBROUTINE INITIALIZES THE ELECTROMAGNETIC FIELDS.


      subroutine INIT_field

      use PIC_variables
      use VLA_variables

      implicit none
      integer :: l1,l2,l3

      real(kind=8) :: h1,h2,h3
      real(kind=8) :: p_pulse_z2

c setting up the fields

      do i3=i3mn-rd3,i3mx+rd3
         h3=i3*dz
         do i2=i2mn-rd2,i2mx+rd2
            h2=i2*dy
            do i1=i1mn-rd1,i1mx+rd1
               h1=i1*dx
               ex(i1,i2,i3)=0.0d0
               ey(i1,i2,i3)=p_pulse_z2(h1,h2+0.5*dy,h3,-0.5*dt)
               ey(i1,i2,i3)=0.0d0
               ez(i1,i2,i3)=0.0d0
               bx(i1,i2,i3)=p_pulse_z2(h1,h2+0.5*dy,h3+0.5*dz,0.0*dt)
               bx(i1,i2,i3)=0.0d0
               by(i1,i2,i3)=0.0d0
               bz(i1,i2,i3)=0.0d0

               jxi(i1,i2,i3)=0.0d0
               jyi(i1,i2,i3)=0.0d0
               jzi(i1,i2,i3)=0.0d0

               ne(i1,i2,i3)=0.0d0
               ni(i1,i2,i3)=0.0d0
               nn(i1,i2,i3)=0.0d0
            enddo
         enddo
      enddo

      end subroutine INIT_field

C     set up various things after the main fields have been initialized

      subroutine SETUP_field

      use PIC_variables
      use VLA_variables

      implicit none
      integer :: l1,l2,l3

      real(kind=8) :: dxi,dyi,dzi
      real(kind=8) :: pxi,pyi,pzi
      real(kind=8) :: qni,mni,cni,lni,wni
      real(kind=8) :: xi,yi,zi
      real(kind=8) :: x0,y0,z0,xr,yr,zr,rot
      real(kind=8) :: u,v,w
      real(kind=8) :: h1,h2,h3
      real(kind=8) :: gmx,g0x,g1x,gmy,g0y,g1y,gmz,g0z,g1z
      real(kind=8) :: fni,fnis
      real(kind=8) :: p_pulse_z2

      x0=5.0*1.0e-6/ld
      y0=5.0*1.0e-6/ld
      z0=5.0*1.0e-6/ld
      rot=45.0*3.141592/180.0

c setting up charge and current densities

      if (niloc.gt.0) then
         dxi=1.0/dx
         dyi=1.0/dy
         dzi=1.0/dz
         fnis=alpha*alpha*cori/eta
         do l=1,niloc

            xi=p_niloc(11*l)
            yi=p_niloc(11*l+1)
            zi=p_niloc(11*l+2)
            pxi=p_niloc(11*l+3)
            pyi=p_niloc(11*l+4)
            pzi=p_niloc(11*l+5)
            qni=p_niloc(11*l+6)
            mni=p_niloc(11*l+7)
            cni=p_niloc(11*l+8)
            lni=p_niloc(11*l+9)
            wni=p_niloc(11*l+10)

            u=xi*dxi
            v=yi*dyi
            w=zi*dzi
            l1=nint(u)
            l2=nint(v)
            l3=nint(w)
            h1=l1-u
            h2=l2-v
            h3=l3-w

            gmx=0.5*(1.5-abs(h1-1.0))*(1.5-abs(h1-1.0))
            g0x=0.75-abs(h1)*abs(h1)
            g1x=0.5*(1.5-abs(h1+1.0))*(1.5-abs(h1+1.0))
            gmy=0.5*(1.5-abs(h2-1.0))*(1.5-abs(h2-1.0))
            g0y=0.75-abs(h2)*abs(h2)
            g1y=0.5*(1.5-abs(h2+1.0))*(1.5-abs(h2+1.0))
            gmz=0.5*(1.5-abs(h3-1.0))*(1.5-abs(h3-1.0))
            g0z=0.75-abs(h3)*abs(h3)
            g1z=0.5*(1.5-abs(h3+1.0))*(1.5-abs(h3+1.0))

            if (qni.lt.0.0) then
               fni=qni*wni*fnis
               ne(l1-1,l2-1,l3-1)=ne(l1-1,l2-1,l3-1)+fni*gmx*gmy*gmz
               ne(l1,l2-1,l3-1)=ne(l1,l2-1,l3-1)+fni*g0x*gmy*gmz
               ne(l1+1,l2-1,l3-1)=ne(l1+1,l2-1,l3-1)+fni*g1x*gmy*gmz
               ne(l1-1,l2,l3-1)=ne(l1-1,l2,l3-1)+fni*gmx*g0y*gmz
               ne(l1,l2,l3-1)=ne(l1,l2,l3-1)+fni*g0x*g0y*gmz
               ne(l1+1,l2,l3-1)=ne(l1+1,l2,l3-1)+fni*g1x*g0y*gmz
               ne(l1-1,l2+1,l3-1)=ne(l1-1,l2+1,l3-1)+fni*gmx*g1y*gmz
               ne(l1,l2+1,l3-1)=ne(l1,l2+1,l3-1)+fni*g0x*g1y*gmz
               ne(l1+1,l2+1,l3-1)=ne(l1+1,l2+1,l3-1)+fni*g1x*g1y*gmz
               ne(l1-1,l2-1,l3)=ne(l1-1,l2-1,l3)+fni*gmx*gmy*g0z
               ne(l1,l2-1,l3)=ne(l1,l2-1,l3)+fni*g0x*gmy*g0z
               ne(l1+1,l2-1,l3)=ne(l1+1,l2-1,l3)+fni*g1x*gmy*g0z
               ne(l1-1,l2,l3)=ne(l1-1,l2,l3)+fni*gmx*g0y*g0z
               ne(l1,l2,l3)=ne(l1,l2,l3)+fni*g0x*g0y*g0z
               ne(l1+1,l2,l3)=ne(l1+1,l2,l3)+fni*g1x*g0y*g0z
               ne(l1-1,l2+1,l3)=ne(l1-1,l2+1,l3)+fni*gmx*g1y*g0z
               ne(l1,l2+1,l3)=ne(l1,l2+1,l3)+fni*g0x*g1y*g0z
               ne(l1+1,l2+1,l3)=ne(l1+1,l2+1,l3)+fni*g1x*g1y*g0z
               ne(l1-1,l2-1,l3+1)=ne(l1-1,l2-1,l3+1)+fni*gmx*gmy*g1z
               ne(l1,l2-1,l3+1)=ne(l1,l2-1,l3+1)+fni*g0x*gmy*g1z
               ne(l1+1,l2-1,l3+1)=ne(l1+1,l2-1,l3+1)+fni*g1x*gmy*g1z
               ne(l1-1,l2,l3+1)=ne(l1-1,l2,l3+1)+fni*gmx*g0y*g1z
               ne(l1,l2,l3+1)=ne(l1,l2,l3+1)+fni*g0x*g0y*g1z
               ne(l1+1,l2,l3+1)=ne(l1+1,l2,l3+1)+fni*g1x*g0y*g1z
               ne(l1-1,l2+1,l3+1)=ne(l1-1,l2+1,l3+1)+fni*gmx*g1y*g1z
               ne(l1,l2+1,l3+1)=ne(l1,l2+1,l3+1)+fni*g0x*g1y*g1z
               ne(l1+1,l2+1,l3+1)=ne(l1+1,l2+1,l3+1)+fni*g1x*g1y*g1z
            else if (qni.gt.0.0) then
               fni=qni*wni*fnis
               ni(l1-1,l2-1,l3-1)=ni(l1-1,l2-1,l3-1)+fni*gmx*gmy*gmz
               ni(l1,l2-1,l3-1)=ni(l1,l2-1,l3-1)+fni*g0x*gmy*gmz
               ni(l1+1,l2-1,l3-1)=ni(l1+1,l2-1,l3-1)+fni*g1x*gmy*gmz
               ni(l1-1,l2,l3-1)=ni(l1-1,l2,l3-1)+fni*gmx*g0y*gmz
               ni(l1,l2,l3-1)=ni(l1,l2,l3-1)+fni*g0x*g0y*gmz
               ni(l1+1,l2,l3-1)=ni(l1+1,l2,l3-1)+fni*g1x*g0y*gmz
               ni(l1-1,l2+1,l3-1)=ni(l1-1,l2+1,l3-1)+fni*gmx*g1y*gmz
               ni(l1,l2+1,l3-1)=ni(l1,l2+1,l3-1)+fni*g0x*g1y*gmz
               ni(l1+1,l2+1,l3-1)=ni(l1+1,l2+1,l3-1)+fni*g1x*g1y*gmz
               ni(l1-1,l2-1,l3)=ni(l1-1,l2-1,l3)+fni*gmx*gmy*g0z
               ni(l1,l2-1,l3)=ni(l1,l2-1,l3)+fni*g0x*gmy*g0z
               ni(l1+1,l2-1,l3)=ni(l1+1,l2-1,l3)+fni*g1x*gmy*g0z
               ni(l1-1,l2,l3)=ni(l1-1,l2,l3)+fni*gmx*g0y*g0z
               ni(l1,l2,l3)=ni(l1,l2,l3)+fni*g0x*g0y*g0z
               ni(l1+1,l2,l3)=ni(l1+1,l2,l3)+fni*g1x*g0y*g0z
               ni(l1-1,l2+1,l3)=ni(l1-1,l2+1,l3)+fni*gmx*g1y*g0z
               ni(l1,l2+1,l3)=ni(l1,l2+1,l3)+fni*g0x*g1y*g0z
               ni(l1+1,l2+1,l3)=ni(l1+1,l2+1,l3)+fni*g1x*g1y*g0z
               ni(l1-1,l2-1,l3+1)=ni(l1-1,l2-1,l3+1)+fni*gmx*gmy*g1z
               ni(l1,l2-1,l3+1)=ni(l1,l2-1,l3+1)+fni*g0x*gmy*g1z
               ni(l1+1,l2-1,l3+1)=ni(l1+1,l2-1,l3+1)+fni*g1x*gmy*g1z
               ni(l1-1,l2,l3+1)=ni(l1-1,l2,l3+1)+fni*gmx*g0y*g1z
               ni(l1,l2,l3+1)=ni(l1,l2,l3+1)+fni*g0x*g0y*g1z
               ni(l1+1,l2,l3+1)=ni(l1+1,l2,l3+1)+fni*g1x*g0y*g1z
               ni(l1-1,l2+1,l3+1)=ni(l1-1,l2+1,l3+1)+fni*gmx*g1y*g1z
               ni(l1,l2+1,l3+1)=ni(l1,l2+1,l3+1)+fni*g0x*g1y*g1z
               ni(l1+1,l2+1,l3+1)=ni(l1+1,l2+1,l3+1)+fni*g1x*g1y*g1z
            else if (qni.eq.0.0) then
               fni=wni*fnis
               nn(l1-1,l2-1,l3-1)=nn(l1-1,l2-1,l3-1)+fni*gmx*gmy*gmz
               nn(l1,l2-1,l3-1)=nn(l1,l2-1,l3-1)+fni*g0x*gmy*gmz
               nn(l1+1,l2-1,l3-1)=nn(l1+1,l2-1,l3-1)+fni*g1x*gmy*gmz
               nn(l1-1,l2,l3-1)=nn(l1-1,l2,l3-1)+fni*gmx*g0y*gmz
               nn(l1,l2,l3-1)=nn(l1,l2,l3-1)+fni*g0x*g0y*gmz
               nn(l1+1,l2,l3-1)=nn(l1+1,l2,l3-1)+fni*g1x*g0y*gmz
               nn(l1-1,l2+1,l3-1)=nn(l1-1,l2+1,l3-1)+fni*gmx*g1y*gmz
               nn(l1,l2+1,l3-1)=nn(l1,l2+1,l3-1)+fni*g0x*g1y*gmz
               nn(l1+1,l2+1,l3-1)=nn(l1+1,l2+1,l3-1)+fni*g1x*g1y*gmz
               nn(l1-1,l2-1,l3)=nn(l1-1,l2-1,l3)+fni*gmx*gmy*g0z
               nn(l1,l2-1,l3)=nn(l1,l2-1,l3)+fni*g0x*gmy*g0z
               nn(l1+1,l2-1,l3)=nn(l1+1,l2-1,l3)+fni*g1x*gmy*g0z
               nn(l1-1,l2,l3)=nn(l1-1,l2,l3)+fni*gmx*g0y*g0z
               nn(l1,l2,l3)=nn(l1,l2,l3)+fni*g0x*g0y*g0z
               nn(l1+1,l2,l3)=nn(l1+1,l2,l3)+fni*g1x*g0y*g0z
               nn(l1-1,l2+1,l3)=nn(l1-1,l2+1,l3)+fni*gmx*g1y*g0z
               nn(l1,l2+1,l3)=nn(l1,l2+1,l3)+fni*g0x*g1y*g0z
               nn(l1+1,l2+1,l3)=nn(l1+1,l2+1,l3)+fni*g1x*g1y*g0z
               nn(l1-1,l2-1,l3+1)=nn(l1-1,l2-1,l3+1)+fni*gmx*gmy*g1z
               nn(l1,l2-1,l3+1)=nn(l1,l2-1,l3+1)+fni*g0x*gmy*g1z
               nn(l1+1,l2-1,l3+1)=nn(l1+1,l2-1,l3+1)+fni*g1x*gmy*g1z
               nn(l1-1,l2,l3+1)=nn(l1-1,l2,l3+1)+fni*gmx*g0y*g1z
               nn(l1,l2,l3+1)=nn(l1,l2,l3+1)+fni*g0x*g0y*g1z
               nn(l1+1,l2,l3+1)=nn(l1+1,l2,l3+1)+fni*g1x*g0y*g1z
               nn(l1-1,l2+1,l3+1)=nn(l1-1,l2+1,l3+1)+fni*gmx*g1y*g1z
               nn(l1,l2+1,l3+1)=nn(l1,l2+1,l3+1)+fni*g0x*g1y*g1z
               nn(l1+1,l2+1,l3+1)=nn(l1+1,l2+1,l3+1)+fni*g1x*g1y*g1z
            endif

            xr=xi
            yr=cos(rot)*(yi-y0)-sin(rot)*(zi-z0)+y0
            zr=cos(rot)*(zi-z0)+sin(rot)*(yi-y0)+z0

            if (((r3n-1)*dz.le.zr.and.zr.le.(r3x+1)*dz).and.
     &          ((r2n-1)*dy.le.yr.and.yr.le.(r2x+1)*dy).and.
     &          ((r1n-1)*dx.le.xr.and.xr.le.(r1x+1)*dx)) then

               p_niloc(11*l+9)=0.0

            endif

         enddo
      endif


      call PIC_fax(ne)
      call PIC_fay(ne)
      call PIC_faz(ne)
      call PIC_fax(ni)
      call PIC_fay(ni)
      call PIC_faz(ni)
      call PIC_fax(nn)
      call PIC_fay(nn)
      call PIC_faz(nn)


c setting up time-averaged fields


      ext=0.0d0
      eyt=0.0d0
      ezt=0.0d0
      bxt=0.0d0
      byt=0.0d0
      bzt=0.0d0
      jxit=0.0d0
      jyit=0.0d0
      jzit=0.0d0
      net=0.0d0
      nit=0.0d0
      nnt=0.0d0

#ifdef USE_PML
      call PML_setup_fields
#endif

      cpu=0.0d0
      cpu_ary=0.0d0

      end subroutine SETUP_field

#ifdef USE_PML
c     ----------------------------------------------------------------------
      subroutine PML_alloc_field
c     ----------------------------------------------------------------------

      use PIC_variables
      use VLA_variables

      implicit none

c kay: potentials      
      allocate(vphi0(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))
      allocate(vphi1(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))
      allocate(vax0(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))
      allocate(vay0(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))
      allocate(vaz0(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))
      allocate(vax1(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))
      allocate(vay1(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))
      allocate(vaz1(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))

      ! additional d field  - added by ab

      allocate(dvx(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))
      allocate(dvy(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))
      allocate(dvz(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))


      ! additional h field - added by ab

      allocate(hx(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))
      allocate(hy(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))
      allocate(hz(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))

      ! electric and magnetic permittivity tensors - added by ab

      allocate(eps(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))
      allocate(mu(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))

      allocate(rhoi(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &             i3mn-rd3:i3mx+rd3))

      allocate(vphit(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(vaxt(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(vayt(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(vazt(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(hxt(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))    ! ab
      allocate(hyt(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))    ! ab
      allocate(hzt(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))    ! ab
      allocate(dxt(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))    ! ab
      allocate(dyt(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))    ! ab
      allocate(dzt(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))    ! ab
      allocate(ex2t(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(ey2t(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(ez2t(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(bx2t(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(by2t(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(bz2t(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(hx2t(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))   ! ab
      allocate(hy2t(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))   ! ab
      allocate(hz2t(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))   ! ab
      allocate(dx2t(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))   ! ab
      allocate(dy2t(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))   ! ab
      allocate(dz2t(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))   ! ab
      allocate(rhoit(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(jxexit(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(jyeyit(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(jzezit(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(poyxt(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(poyyt(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(poyzt(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))

      allocate(exsp(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))  ! test ab
      allocate(eysp(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))  ! test ab
      allocate(ezsp(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))  ! test ab
      allocate(jxsp(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))  ! test ab
      allocate(jysp(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))  ! test ab
      allocate(jzsp(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))  ! test ab

      vphi0=0.0d0
      vphi1=0.0d0
      vax0=0.0d0
      vay0=0.0d0
      vaz0=0.0d0
      vax1=0.0d0
      vay1=0.0d0
      vaz1=0.0d0
      dvx=0.0d0
      dvy=0.0d0
      dvz=0.0d0
      hx=0.0d0
      hy=0.0d0
      hz=0.0d0
      eps = 1.0
      mu = 1.0
      rhoi=0.0d0

      end subroutine PML_alloc_field

c     ----------------------------------------------------------------------
      subroutine PML_setup_field
c     ----------------------------------------------------------------------

      use PIC_variables
      use VLA_variables

      implicit none

      vphi0=0.0d0
      vphi1=0.0d0
      vax0=0.0d0
      vay0=0.0d0
      vaz0=0.0d0
      vax1=0.0d0
      vay1=0.0d0
      vaz1=0.0d0
      
      ! initialisation of additional fields - added by ab
      ! because of normalization these are equal to e and b fields

      dvx=ex
      dvy=ey
      dvz=ez
      hx=bx
      hy=by
      hz=bz
      
      ! initialisation of permittivity tensors
      ! due to scaling these are dimensionless
      
      eps = 1.0
      mu = 1.0
      
      rhoi=0.0d0

      vphit=0.0d0
      vaxt=0.0d0
      vayt=0.0d0
      vazt=0.0d0
      hxt=0.0d0       ! ab
      hyt=0.0d0       ! ab
      hzt=0.0d0       ! ab
      dxt=0.0d0       ! ab
      dyt=0.0d0       ! ab
      dzt=0.0d0       ! ab
      ex2t=0.0d0
      ey2t=0.0d0
      ez2t=0.0d0
      bx2t=0.0d0
      by2t=0.0d0
      bz2t=0.0d0
      hx2t=0.0d0      ! ab
      hy2t=0.0d0      ! ab
      hz2t=0.0d0      ! ab
      dx2t=0.0d0      ! ab
      dy2t=0.0d0      ! ab
      dz2t=0.0d0      ! ab
      rhoit=0.0d0
      jxexit=0.0d0
      jyeyit=0.0d0
      jzezit=0.0d0
      poyxt=0.0d0
      poyyt=0.0d0
      poyzt=0.0d0
      exsp=0.0d0      ! test ab
      eysp=0.0d0      ! test ab
      ezsp=0.0d0      ! test ab
      jxsp=0.0d0      ! test ab
      jysp=0.0d0      ! test ab
      jzsp=0.0d0      ! test ab

      fluxit=0.0d0
      fluxot=0.0d0

      enEXt=0.0d0
      enEYt=0.0d0
      enEZt=0.0d0
      enBXt=0.0d0
      enBYt=0.0d0
      enBZt=0.0d0
      enHXt=0.0d0   ! ab
      enHYt=0.0d0   ! ab
      enHZt=0.0d0   ! ab
      enDXt=0.0d0   ! ab
      enDYt=0.0d0   ! ab
      enDZt=0.0d0   ! ab

      ent=0.0d0
      poxt=0.0d0
      poyt=0.0d0
      pozt=0.0d0
      jet=0.0d0

      sum_fed=0.0d0
      sum_ped=0.0d0

! deriving pml attenuation coefficients

! not moved by half step
      
      ! x

      do i1 = i1mn-rd1, i1mx+rd1
         kappax(i1) = 1.0
         sigmax(i1) = 0.0
         cxp(i1) = 2*eps0*kappax(i1)+sigmax(i1)*dt
         cxm(i1) = 2*eps0*kappax(i1)-sigmax(i1)*dt
         fbx(i1) = 2*eps0*kappax(i1)
         fcx(i1) = cxm(i1)/cxp(i1)
         fdx(i1) = 2*eps0*dt/cxp(i1)
         fex(i1) = 1.0/cxp(i1)
      end do

      ! y
      
       do i2 = i2mn-rd2, i2mx+rd2
         kappay(i2) = 1.0
         sigmay(i2) = 0.0
         cyp(i2) = 2*eps0*kappay(i2)+sigmay(i2)*dt
         cym(i2) = 2*eps0*kappay(i2)-sigmay(i2)*dt
         fby(i2) = 2*eps0*kappay(i2)
         fcy(i2) = cym(i2)/cyp(i2)
         fdy(i2) = 2*eps0*dt/cyp(i2)
         fey(i2) = 1.0/cyp(i2)
      end do

      ! z

      do i3 = i3mn-rd3, i3mx+rd3
         kappaz(i3) = 1.0
         sigmaz(i3) = 0.0
         czp(i3) = 2*eps0*kappaz(i3)+sigmaz(i3)*dt
         czm(i3) = 2*eps0*kappaz(i3)-sigmaz(i3)*dt
         fbz(i3) = 2*eps0*kappaz(i3)
         fcz(i3) = czm(i3)/czp(i3)
         fdz(i3) = 2*eps0*dt/czp(i3)
         fez(i3) = 1.0/czp(i3)
      end do


! moved by half step

      ! x

      do i1 = i1mn-rd1, i1mx+rd1
         kappax(i1) = 1.0
         sigmax(i1) = 0.0
         bxp(i1) = 2*eps0*kappax(i1)+sigmax(i1)*dt
         bxm(i1) = 2*eps0*kappax(i1)-sigmax(i1)*dt
         gbx(i1) = 2*eps0*kappax(i1)
         gcx(i1) = bxm(i1)/bxp(i1)
         gdx(i1) = 2*eps0*dt/bxp(i1)
         gex(i1) = 1.0/bxp(i1)
      end do

      ! y

      do i2 = i2mn-rd2, i2mx+rd2
         kappay(i2) = 1.0
         sigmay(i2) = 0.0
         byp(i2) = 2*eps0*kappay(i2)+sigmay(i2)*dt
         bym(i2) = 2*eps0*kappay(i2)-sigmay(i2)*dt
         gby(i2) = 2*eps0*kappay(i2)
         gcy(i2) = bym(i2)/byp(i2)
         gdy(i2) = 2*eps0*dt/byp(i2)
         gey(i2) = 1.0/byp(i2)
      end do

      ! z

      do i3 = i3mn-rd3, i3mx+rd3
         kappaz(i3) = 1.0
         sigmaz(i3) = 0.0
         bzp(i3) = 2*eps0*kappaz(i3)+sigmaz(i3)*dt
         bzm(i3) = 2*eps0*kappaz(i3)-sigmaz(i3)*dt
         gbz(i3) = 2*eps0*kappaz(i3)
         gcz(i3) = bzm(i3)/bzp(i3)
         gdz(i3) = 2*eps0*dt/bzp(i3)
         gez(i3) = 1.0/bzp(i3)
      end do

      end subroutine PML_setup_field
#endif

