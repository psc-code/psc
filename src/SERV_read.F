#include "config.h"

c CHECKPOINTING ROUTINE FOR DATA WRITEOUT

      subroutine SERV_read_1(nstart_,niloc_)

      use VLA_variables
      use OUT_params
      use OUT_particles
      use OUT_fields
      use PIC_bin_collision ! for I/O params
      use OUT_poyn
      use PSC_globals
      use PSC_patch

      implicit none

      integer nstart_,niloc_
      character*5 node

      call SERV_labelgen(mpe,node)


      open(10,file=trim(data_chk)//'/'//node//'CPNEW',
     &     access='sequential',form='unformatted')

      read(10) nstart_
      read(10) shift_c,shift_z
      read(10) nprpartio
      read(10) tmnvfo,tmxvfo
      read(10) tmnvpo,tmxvpo

      if (nprparti.lt.nstart_) nprparti=nprpartio
      if (tmnvf.lt.tmnvfo) tmnvf=tmnvfo
      if (tmxvf.lt.tmxvfo) tmxvf=tmxvfo
      if (tmnvp.lt.tmnvpo) tmnvp=tmnvpo
      if (tmxvp.lt.tmxvpo) tmxvp=tmxvpo

      call OUT_poyn_serv_read

      read(10) niloc_,cori
      read(10) i1mn,i1mx,i2mn,i2mx,i3mn,i3mx 

      end subroutine SERV_read_1



      subroutine SERV_read_2(niloc, p_niloc,
     &     jxi, jyi, jzi, ex, ey, ez, hx, hy, hz)

      use OUT_params
      use OUT_particles
      use OUT_fields
      use PIC_bin_collision ! for I/O params
      use PSC_domain, only: mpe
      use PSC_patch

      implicit none

      integer      :: niloc
      real(kind=8) :: p_niloc(0:*)
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

      close(10)

      end subroutine SERV_read_2
