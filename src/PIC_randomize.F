
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


