#include "config.h"

c CHECKPOINTING ROUTINE FOR DATA WRITEOUT

      subroutine SERV_read

      use VLA_variables, only: cpu,cpu_ary,npe,nstart
      use PIC_variables, only: niloc,p_niloc
      use fields

      implicit none

      allocate(cpu(1:5))
      allocate(cpu_ary(1:5,0:npe-1))

      call SERV_read_1(nstart,niloc)

      call ALLOC_particles(niloc)
      call FIELDS_alloc

      call SERV_read_2(niloc, p_niloc,
     &     ne, ni, nn, jxi, jyi, jzi, ex, ey, ez, hx, hy, hz)

      end subroutine SERV_read



      subroutine SERV_read_1(nstart_,niloc_)

      use VLA_variables
      use PIC_variables

      implicit none

      integer nstart_,niloc_
      character*5 node

      call SERV_labelgen(mpe,node)


      open(10,file=trim(data_chk)//'/'//node//'CPNEW',
     &     access='sequential',form='unformatted')

      read(10) nstart_
      read(10) shift_c,shift_z
      read(10) nprfo,nprco
      read(10) nprpartio
      read(10) tmnvfo,tmxvfo
      read(10) tmnvpo,tmxvpo
      read(10) tmnvco,tmxvco

      if (nprf.lt.nstart_) nprf=nprfo 
      if (nprc.lt.nstart_) nprc=nprco
      if (nprparti.lt.nstart_) nprparti=nprpartio
      if (tmnvf.lt.tmnvfo) tmnvf=tmnvfo
      if (tmxvf.lt.tmxvfo) tmxvf=tmxvfo
      if (tmnvp.lt.tmnvpo) tmnvp=tmnvpo
      if (tmxvp.lt.tmxvpo) tmxvp=tmxvpo
      if (tmnvc.lt.tmnvco) tmnvc=tmnvco
      if (tmxvc.lt.tmxvco) tmxvc=tmxvco

      call OUT_poyn_serv_read

      read(10) niloc_,nialloc,cori
      read(10) i1mn,i1mx,i2mn,i2mx,i3mn,i3mx 

      end subroutine SERV_read_1



      subroutine SERV_read_2(niloc, p_niloc,
     &     ne, ni, nn, jxi, jyi, jzi, ex, ey, ez, hx, hy, hz)

      use VLA_variables, only: i1mn,i2mn,i3mn,i1mx,i2mx,i3mx,
     *     rd1,rd2,rd3,data_chk,mpe,nprf,nprc,shift_c,shift_z,
     *     tmnvc,tmxvc,tmnvf,tmnvp,tmxvf,tmxvp,cpu,cpu_ary,npe

      use PIC_variables, only: cori,nialloc,nprparti

      implicit none

      integer      :: niloc
      real(kind=8) :: p_niloc(0:*)
      real(kind=8),dimension(i1mn-rd1:i1mx+rd1,
     &                       i2mn-rd2:i2mx+rd2,
     &                       i3mn-rd3:i3mx+rd3) :: ne, ni, nn
      real(kind=8),dimension(i1mn-rd1:i1mx+rd1,
     &                       i2mn-rd2:i2mx+rd2,
     &                       i3mn-rd3:i3mx+rd3) :: jxi, jyi, jzi
      real(kind=8),dimension(i1mn-rd1:i1mx+rd1,
     &                       i2mn-rd2:i2mx+rd2,
     &                       i3mn-rd3:i3mx+rd3) :: ex, ey, ez
      real(kind=8),dimension(i1mn-rd1:i1mx+rd1,
     &                       i2mn-rd2:i2mx+rd2,
     &                       i3mn-rd3:i3mx+rd3) :: hx, hy, hz

      integer i1,i2,i3
      integer i1a,i1e,i2a,i2e,i3a,i3e

      do i1=0,11*niloc+10,100
         i1e=min(i1+99,11*niloc+10)
         read(10) (p_niloc(i1a),i1a=i1,i1e)
      enddo

      do i3=i3mn-rd3,i3mx+rd3,100
         i3e=min(i3+99,i3mx+rd3)
         do i2=i2mn-rd2,i2mx+rd2,100
            i2e=min(i2+99,i2mx+rd2)
            do i1=i1mn-rd1,i1mx+rd1,100
               i1e=min(i1+99,i1mx+rd1)
               read(10) (((ne(i1a,i2a,i3a),i1a=i1,i1e),
     &                     i2a=i2,i2e),i3a=i3,i3e)
               read(10) (((ni(i1a,i2a,i3a),i1a=i1,i1e),
     &                     i2a=i2,i2e),i3a=i3,i3e)
               read(10) (((nn(i1a,i2a,i3a),i1a=i1,i1e),
     &                     i2a=i2,i2e),i3a=i3,i3e)
               read(10) (((jxi(i1a,i2a,i3a),i1a=i1,i1e),
     &                     i2a=i2,i2e),i3a=i3,i3e)
               read(10) (((jyi(i1a,i2a,i3a),i1a=i1,i1e),
     &                     i2a=i2,i2e),i3a=i3,i3e)
               read(10) (((jzi(i1a,i2a,i3a),i1a=i1,i1e),
     &                     i2a=i2,i2e),i3a=i3,i3e)
               read(10) (((ex(i1a,i2a,i3a),i1a=i1,i1e),
     &                     i2a=i2,i2e),i3a=i3,i3e)
               read(10) (((ey(i1a,i2a,i3a),i1a=i1,i1e),
     &                     i2a=i2,i2e),i3a=i3,i3e)
               read(10) (((ez(i1a,i2a,i3a),i1a=i1,i1e),
     &                     i2a=i2,i2e),i3a=i3,i3e)
               read(10) (((hx(i1a,i2a,i3a),i1a=i1,i1e),
     &                     i2a=i2,i2e),i3a=i3,i3e)
               read(10) (((hy(i1a,i2a,i3a),i1a=i1,i1e),
     &                     i2a=i2,i2e),i3a=i3,i3e)
               read(10) (((hz(i1a,i2a,i3a),i1a=i1,i1e),
     &                     i2a=i2,i2e),i3a=i3,i3e)
            enddo
         enddo
      enddo

      call OUT_FIELDS_serv_read

      cpu=0.0d0
      cpu_ary=0.0d0

      close(10)

      end subroutine SERV_read_2
