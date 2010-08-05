#include "config.h"

c CHECKPOINTING ROUTINE FOR DATA WRITEOUT

      subroutine SERV_read

      use VLA_variables
      use PIC_variables
      use fields

      implicit none

      integer i1a,i1e,i2a,i2e,i3a,i3e
      character*5 node


      call SERV_labelgen(mpe,node)


      open(10,file=trim(data_chk)//'/'//node//'CPNEW',
     &     access='sequential',form='unformatted')

      read(10) nstart
      read(10) shift_c,shift_z
      read(10) nprfo,nprco
      read(10) nprpartio
      read(10) tmnvfo,tmxvfo
      read(10) tmnvpo,tmxvpo
      read(10) tmnvco,tmxvco

      if (nprf.lt.nstart) nprf=nprfo 
      if (nprc.lt.nstart) nprc=nprco
      if (nprparti.lt.nstart) nprparti=nprpartio
      if (tmnvf.lt.tmnvfo) tmnvf=tmnvfo
      if (tmxvf.lt.tmxvfo) tmxvf=tmxvfo
      if (tmnvp.lt.tmnvpo) tmnvp=tmnvpo
      if (tmxvp.lt.tmxvpo) tmxvp=tmxvpo
      if (tmnvc.lt.tmnvco) tmnvc=tmnvco
      if (tmxvc.lt.tmxvco) tmxvc=tmxvco

      read(10) fluxit,fluxot
      read(10) ent,poxt,poyt,pozt,jet
      read(10) enEXt,enEYt,enEZt
      read(10) enBXt,enBYt,enBZt
#ifdef USE_PML      
      read(10) enHXt,enHYt,enHZt
#endif
      read(10) sum_fed,sum_ped


      read(10) niloc,nialloc,cori
      read(10) i1mn,i1mx,i2mn,i2mx,i3mn,i3mx 


      call ALLOC_particles(niloc)

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

      call OUT_alloc_fields
#ifdef USE_PML
      call OUT_alloc_fields_pml
#endif

      allocate(cpu(1:5))
      allocate(cpu_ary(1:5,0:npe-1))


      p_niloc=0.0
      ex=0.0d0
      ey=0.0d0
      ez=0.0d0
      bx=0.0d0
      by=0.0d0
      bz=0.0d0
      jxi=0.0d0
      jyi=0.0d0
      jzi=0.0d0
      ne=0.0d0
      ni=0.0d0
      nn=0.0d0

      call OUT_setup_fields
#ifdef USE_PML
      call OUT_setup_fields_pml
#endif

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
               read(10) (((bx(i1a,i2a,i3a),i1a=i1,i1e),
     &                     i2a=i2,i2e),i3a=i3,i3e)
               read(10) (((by(i1a,i2a,i3a),i1a=i1,i1e),
     &                     i2a=i2,i2e),i3a=i3,i3e)
               read(10) (((bz(i1a,i2a,i3a),i1a=i1,i1e),
     &                     i2a=i2,i2e),i3a=i3,i3e)
            enddo
         enddo
      enddo

      call OUT_serv_read_fields

#ifdef USE_PML
      call OUT_serv_read_fields_pml
#endif

      cpu=0.0d0
      cpu_ary=0.0d0

      close(10)

      end subroutine SERV_read
