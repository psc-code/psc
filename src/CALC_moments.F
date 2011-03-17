
#include "config.h"

c     ----------------------------------------------------------------------
      subroutine CALC_densities(niloc, p_niloc, ne, ni, nn)
c     ----------------------------------------------------------------------

      use PSC_globals
      use PSC_patch

      implicit none

      integer :: niloc
      real(kind=8) :: p_niloc(0:*)
      real(kind=8),dimension(i1mn-rd1:i1mx+rd1,
     &                       i2mn-rd2:i2mx+rd2,
     &                       i3mn-rd3:i3mx+rd3) :: ne, ni, nn

      integer :: l,l1,l2,l3

      real(kind=8) :: dxi,dyi,dzi
      real(kind=8) :: pxi,pyi,pzi
      real(kind=8) :: qni,mni,cni,lni,wni
      real(kind=8) :: xi,yi,zi
      real(kind=8) :: u,v,w
      real(kind=8) :: h1,h2,h3
      real(kind=8) :: gmx,g0x,g1x,gmy,g0y,g1y,gmz,g0z,g1z
      real(kind=8) :: fni,fnis

      dxi=1.0/dx
      dyi=1.0/dy
      dzi=1.0/dz
      fnis=alpha*alpha*cori/eta

      ne = 0.
      ni = 0.
      nn = 0.

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
      enddo
      
      call PIC_fax(ne)
      call PIC_fay(ne)
      call PIC_faz(ne)
      call PIC_fax(ni)
      call PIC_fay(ni)
      call PIC_faz(ni)
      call PIC_fax(nn)
      call PIC_fay(nn)
      call PIC_faz(nn)

      end subroutine CALC_densities
