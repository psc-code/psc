
c     ======================================================================
      subroutine PIC_find_cell_indices(niloc, p_niloc)
c     ======================================================================

      use VLA_variables
      use PSC_domain
      use PSC_patch

      implicit none

      integer      :: niloc
      real(kind=8) :: p_niloc(0:*)

      integer :: ll,i1,i2,i3
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


