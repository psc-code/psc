      subroutine SERV_labelgen(nstart,label)
      implicit none
      
      integer nstart
      character*5 label
      
      write(label,"(i5.5)") nstart

      return
      end
