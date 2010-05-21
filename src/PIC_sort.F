
#include "config.h"

c     ======================================================================
      subroutine PIC_randomize
c     ======================================================================

      use PIC_variables
      use VLA_variables

      implicit none

      real(kind=8) :: ran
      real(kind=8) :: xh,yh,zh
      real(kind=8) :: pxh,pyh,pzh
      real(kind=8) :: qnh,mnh,cnh,lnh,wnh
      real(kind=8),allocatable,dimension(:) :: rndmv


      allocate(rndmv(1:niloc))
      call random_number(rndmv)
      
      do l=1,niloc-1
         
         ran=rndmv(l)
         j=l+int((niloc+1-l)*ran)
         
         xh=p_niloc(11*l)
         yh=p_niloc(11*l+1)
         zh=p_niloc(11*l+2)
         pxh=p_niloc(11*l+3)
         pyh=p_niloc(11*l+4)
         pzh=p_niloc(11*l+5)
         qnh=p_niloc(11*l+6)
         mnh=p_niloc(11*l+7)
         cnh=p_niloc(11*l+8)
         lnh=p_niloc(11*l+9)
         wnh=p_niloc(11*l+10)
         
         p_niloc(11*l)=p_niloc(11*j)
         p_niloc(11*l+1)=p_niloc(11*j+1)
         p_niloc(11*l+2)=p_niloc(11*j+2)
         p_niloc(11*l+3)=p_niloc(11*j+3)
         p_niloc(11*l+4)=p_niloc(11*j+4)
         p_niloc(11*l+5)=p_niloc(11*j+5)
         p_niloc(11*l+6)=p_niloc(11*j+6)
         p_niloc(11*l+7)=p_niloc(11*j+7)
         p_niloc(11*l+8)=p_niloc(11*j+8)
         p_niloc(11*l+9)=p_niloc(11*j+9)
         p_niloc(11*l+10)=p_niloc(11*j+10)
         
         p_niloc(11*j)=xh
         p_niloc(11*j+1)=yh
         p_niloc(11*j+2)=zh
         p_niloc(11*j+3)=pxh
         p_niloc(11*j+4)=pyh
         p_niloc(11*j+5)=pzh
         p_niloc(11*j+6)=qnh
         p_niloc(11*j+7)=mnh
         p_niloc(11*j+8)=cnh
         p_niloc(11*j+9)=lnh
         p_niloc(11*j+10)=wnh
         
      enddo
      
      deallocate(rndmv)

      end subroutine PIC_randomize


c     ======================================================================
      subroutine PIC_find_cell_indices
c     ======================================================================

      use PIC_variables
      use VLA_variables

      implicit none

      integer :: ll
      integer :: rds1,rds2,rds3
      real(kind=8) :: u,v,w
      real(kind=8) :: xh,yh,zh
      real(kind=8) :: dxi,dyi,dzi

      dxi=1.0/dx
      dyi=1.0/dy
      dzi=1.0/dz

      do ll=1,niloc

         xh=p_niloc(11*ll)
         yh=p_niloc(11*ll+1)
         zh=p_niloc(11*ll+2)

         u=xh*dxi
         v=yh*dyi
         w=zh*dzi
         i1=int(u)
         i2=int(v)
         i3=int(w)
         if (u<0.0) i1=i1-1
         if (v<0.0) i2=i2-1
         if (w<0.0) i3=i3-1


c  -rd <= i1 <= i1mx+rd-1
c  -rd <= i2 <= i2mx+rd-1
c  -rd <= i3 <= i3mx+rd-1


         rds1=rd1
         rds2=rd2
         rds3=rd3
         
         if (i1n==i1x) rds1=0
         if (i2n==i2x) rds2=0
         if (i3n==i3x) rds3=0
         
         p_niloc(11*ll+8)=1
     &                    +(i1-i1mn+rds1)
     &                    +(i1mx-i1mn+2*rds1+1)
     &                    *(i2-i2mn+rds2)
     &                    +(i1mx-i1mn+2*rds1+1)
     &                    *(i2mx-i2mn+2*rds2+1)
     &                    *(i3-i3mn+rds3)

      enddo

      end subroutine PIC_find_cell_indices

c     ======================================================================
      subroutine PIC_sort
c     ======================================================================
c     THIS SUBROUTINE SORTS THE PARTICLES TO CELLS


      use PIC_variables
      use VLA_variables

      implicit none

      real(kind=8) :: s_cpuh


      s_cpuh=0.0d0
      call SERV_systime(cpug)

      call PIC_randomize
      call PIC_find_cell_indices
#ifdef USE_C_INTERFACE
      call C_sort(niloc, p_niloc, 99)
#else
      call PIC_sort_1(niloc, p_niloc)
#endif

      call SERV_systime(cpuh)
      s_cpuh=s_cpuh+cpuh-cpug

      cpucomp=cpucomp+s_cpuh

      end subroutine PIC_sort
