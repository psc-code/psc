
c     ======================================================================
      subroutine PIC_sort(niloc, p_niloc)
c     ======================================================================
c     THIS SUBROUTINE SORTS THE PARTICLES TO CELLS

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
      
      end subroutine PIC_sort
