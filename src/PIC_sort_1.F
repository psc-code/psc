
c     ======================================================================
      subroutine PIC_randomize(niloc, p_niloc)
c     ======================================================================

      implicit none

      integer :: niloc
      real(kind=8) :: p_niloc(0:*)

      integer :: l, j
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
      subroutine PIC_find_cell_indices(niloc, p_niloc)
c     ======================================================================

      use VLA_variables

      implicit none

      integer      :: niloc
      real(kind=8) :: p_niloc(0:*)

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
      subroutine PIC_sort_1(niloc, p_niloc)
c     ======================================================================
c     SORT PARTICLES AFTER CELL NUMBER

      implicit none

      integer :: niloc
      real(kind=8) :: p_niloc(0:*)

      integer :: ii,jj,ll,ir
      real(kind=8) :: rr0,rr1,rr2,rr3,rr4,rr5,rr6,rr7,rra,rr9,rr10

      if (niloc.lt.2) goto 30
      ll=niloc/2+1
      ir=niloc
 10   continue

      if (ll.gt.1) then
         ll=ll-1
         rr0=p_niloc(11*ll+0)
         rr1=p_niloc(11*ll+1)
         rr2=p_niloc(11*ll+2)
         rr3=p_niloc(11*ll+3)
         rr4=p_niloc(11*ll+4)
         rr5=p_niloc(11*ll+5)
         rr6=p_niloc(11*ll+6)
         rr7=p_niloc(11*ll+7)
         rra=p_niloc(11*ll+8)
         rr9=p_niloc(11*ll+9)
         rr10=p_niloc(11*ll+10)
      else
         rr0=p_niloc(11*ir+0)
         rr1=p_niloc(11*ir+1)
         rr2=p_niloc(11*ir+2)
         rr3=p_niloc(11*ir+3)
         rr4=p_niloc(11*ir+4)
         rr5=p_niloc(11*ir+5)
         rr6=p_niloc(11*ir+6)
         rr7=p_niloc(11*ir+7)
         rra=p_niloc(11*ir+8)
         rr9=p_niloc(11*ir+9)
         rr10=p_niloc(11*ir+10)
         p_niloc(11*ir+0)=p_niloc(11)
         p_niloc(11*ir+1)=p_niloc(12)
         p_niloc(11*ir+2)=p_niloc(13)
         p_niloc(11*ir+3)=p_niloc(14)
         p_niloc(11*ir+4)=p_niloc(15)
         p_niloc(11*ir+5)=p_niloc(16)
         p_niloc(11*ir+6)=p_niloc(17)
         p_niloc(11*ir+7)=p_niloc(18)
         p_niloc(11*ir+8)=p_niloc(19)
         p_niloc(11*ir+9)=p_niloc(20)
         p_niloc(11*ir+10)=p_niloc(21)
         ir=ir-1
         if (ir.eq.1) then
            p_niloc(11)=rr0
            p_niloc(12)=rr1
            p_niloc(13)=rr2
            p_niloc(14)=rr3
            p_niloc(15)=rr4
            p_niloc(16)=rr5
            p_niloc(17)=rr6
            p_niloc(18)=rr7
            p_niloc(19)=rra
            p_niloc(20)=rr9
            p_niloc(21)=rr10
            goto 30
         endif
      endif
      
      ii=ll
      jj=ll+ll
      
 20   if (jj.le.ir) then
         if (jj.lt.ir) then
            if (p_niloc(11*jj+8).lt.p_niloc(11*(jj+1)+8)) jj=jj+1
         endif
         
         if (rra.lt.p_niloc(11*jj+8)) then
            p_niloc(11*ii+0)=p_niloc(11*jj+0)
            p_niloc(11*ii+1)=p_niloc(11*jj+1)
            p_niloc(11*ii+2)=p_niloc(11*jj+2)
            p_niloc(11*ii+3)=p_niloc(11*jj+3)
            p_niloc(11*ii+4)=p_niloc(11*jj+4)
            p_niloc(11*ii+5)=p_niloc(11*jj+5)
            p_niloc(11*ii+6)=p_niloc(11*jj+6)
            p_niloc(11*ii+7)=p_niloc(11*jj+7)
            p_niloc(11*ii+8)=p_niloc(11*jj+8)
            p_niloc(11*ii+9)=p_niloc(11*jj+9)
            p_niloc(11*ii+10)=p_niloc(11*jj+10)
            ii=jj
            jj=jj+jj
         else
            jj=ir+1
         endif
         goto 20
      endif
      
      p_niloc(11*ii+0)=rr0
      p_niloc(11*ii+1)=rr1
      p_niloc(11*ii+2)=rr2
      p_niloc(11*ii+3)=rr3
      p_niloc(11*ii+4)=rr4
      p_niloc(11*ii+5)=rr5
      p_niloc(11*ii+6)=rr6
      p_niloc(11*ii+7)=rr7
      p_niloc(11*ii+8)=rra
      p_niloc(11*ii+9)=rr9
      p_niloc(11*ii+10)=rr10
      goto 10
      
 30   continue
      
      end subroutine PIC_sort_1

