#include "config.h"

c CHECKPOINTING ROUTINE FOR DATA WRITEOUT

      subroutine SERV_read

      use VLA_variables
      use PIC_variables
      use fields
      use out_fields

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

      allocate(ext(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(eyt(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(ezt(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(bxt(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(byt(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(bzt(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(net(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(nit(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(nnt(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(jxit(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(jyit(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(jzit(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
#ifdef USE_PML
      allocate(hxt(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(hyt(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(hzt(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(ex2t(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(ey2t(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(ez2t(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(bx2t(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(by2t(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(bz2t(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(hx2t(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(hy2t(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(hz2t(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(jxexit(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(jyeyit(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(jzezit(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(poyxt(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(poyyt(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
      allocate(poyzt(i1mn:i1mx,i2mn:i2mx,i3mn:i3mx))
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

      ext=0.0d0
      eyt=0.0d0
      ezt=0.0d0
      bxt=0.0d0
      byt=0.0d0
      bzt=0.0d0
      net=0.0d0
      nit=0.0d0
      nnt=0.0d0
      jxit=0.0d0
      jyit=0.0d0
      jzit=0.0d0
#ifdef USE_PML
      hxt=0.0d0
      hyt=0.0d0
      hzt=0.0d0
      ex2t=0.0d0
      ey2t=0.0d0
      ez2t=0.0d0
      bx2t=0.0d0
      by2t=0.0d0
      bz2t=0.0d0
      hx2t=0.0d0
      hy2t=0.0d0
      hz2t=0.0d0
      jxexit=0.0d0
      jyeyit=0.0d0
      jzezit=0.0d0
      poyxt=0.0d0
      poyyt=0.0d0
      poyzt=0.0d0
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

      do i3=i3mn,i3mx,100
         i3e=min(i3+99,i3mx)
         do i2=i2mn,i2mx,100
            i2e=min(i2+99,i2mx)
            do i1=i1mn,i1mx,100
               i1e=min(i1+99,i1mx)
               read(10) (((ext(i1a,i2a,i3a),i1a=i1,i1e),
     &                     i2a=i2,i2e),i3a=i3,i3e)
               read(10) (((eyt(i1a,i2a,i3a),i1a=i1,i1e),
     &                     i2a=i2,i2e),i3a=i3,i3e)
               read(10) (((ezt(i1a,i2a,i3a),i1a=i1,i1e),
     &                     i2a=i2,i2e),i3a=i3,i3e)
               read(10) (((bxt(i1a,i2a,i3a),i1a=i1,i1e),
     &                     i2a=i2,i2e),i3a=i3,i3e)
               read(10) (((byt(i1a,i2a,i3a),i1a=i1,i1e),
     &                     i2a=i2,i2e),i3a=i3,i3e)
               read(10) (((bzt(i1a,i2a,i3a),i1a=i1,i1e),
     &                     i2a=i2,i2e),i3a=i3,i3e)
               read(10) (((net(i1a,i2a,i3a),i1a=i1,i1e),
     &                     i2a=i2,i2e),i3a=i3,i3e)
               read(10) (((nit(i1a,i2a,i3a),i1a=i1,i1e),
     &                     i2a=i2,i2e),i3a=i3,i3e)
               read(10) (((nnt(i1a,i2a,i3a),i1a=i1,i1e),
     &                     i2a=i2,i2e),i3a=i3,i3e)
               read(10) (((jxit(i1a,i2a,i3a),i1a=i1,i1e),
     &                     i2a=i2,i2e),i3a=i3,i3e)
               read(10) (((jyit(i1a,i2a,i3a),i1a=i1,i1e),
     &                     i2a=i2,i2e),i3a=i3,i3e)
               read(10) (((jzit(i1a,i2a,i3a),i1a=i1,i1e),
     &                     i2a=i2,i2e),i3a=i3,i3e)
            enddo
         enddo
      enddo

#ifdef USE_PML
      do i3=i3mn,i3mx,100
         i3e=min(i3+99,i3mx)
         do i2=i2mn,i2mx,100
            i2e=min(i2+99,i2mx)
            do i1=i1mn,i1mx,100
               i1e=min(i1+99,i1mx)
               read(10) (((hxt(i1a,i2a,i3a),i1a=i1,i1e),
     &                     i2a=i2,i2e),i3a=i3,i3e)
               read(10) (((hyt(i1a,i2a,i3a),i1a=i1,i1e),
     &                     i2a=i2,i2e),i3a=i3,i3e)
               read(10) (((hzt(i1a,i2a,i3a),i1a=i1,i1e),
     &                     i2a=i2,i2e),i3a=i3,i3e)
               read(10) (((ex2t(i1a,i2a,i3a),i1a=i1,i1e),
     &                     i2a=i2,i2e),i3a=i3,i3e)
               read(10) (((ey2t(i1a,i2a,i3a),i1a=i1,i1e),
     &                     i2a=i2,i2e),i3a=i3,i3e)
               read(10) (((ez2t(i1a,i2a,i3a),i1a=i1,i1e),
     &                     i2a=i2,i2e),i3a=i3,i3e)
               read(10) (((bx2t(i1a,i2a,i3a),i1a=i1,i1e),
     &                     i2a=i2,i2e),i3a=i3,i3e)
               read(10) (((by2t(i1a,i2a,i3a),i1a=i1,i1e),
     &                     i2a=i2,i2e),i3a=i3,i3e)
               read(10) (((bz2t(i1a,i2a,i3a),i1a=i1,i1e),
     &                     i2a=i2,i2e),i3a=i3,i3e)
               read(10) (((hx2t(i1a,i2a,i3a),i1a=i1,i1e),
     &                     i2a=i2,i2e),i3a=i3,i3e)
               read(10) (((hy2t(i1a,i2a,i3a),i1a=i1,i1e),
     &                     i2a=i2,i2e),i3a=i3,i3e)
               read(10) (((hz2t(i1a,i2a,i3a),i1a=i1,i1e),
     &                     i2a=i2,i2e),i3a=i3,i3e)
               read(10) (((jxexit(i1a,i2a,i3a),i1a=i1,i1e),
     &                     i2a=i2,i2e),i3a=i3,i3e)
               read(10) (((jyeyit(i1a,i2a,i3a),i1a=i1,i1e),
     &                     i2a=i2,i2e),i3a=i3,i3e)
               read(10) (((jzezit(i1a,i2a,i3a),i1a=i1,i1e),
     &                     i2a=i2,i2e),i3a=i3,i3e)
               read(10) (((poyxt(i1a,i2a,i3a),i1a=i1,i1e),
     &                     i2a=i2,i2e),i3a=i3,i3e)
               read(10) (((poyyt(i1a,i2a,i3a),i1a=i1,i1e),
     &                     i2a=i2,i2e),i3a=i3,i3e)
               read(10) (((poyzt(i1a,i2a,i3a),i1a=i1,i1e),
     &                     i2a=i2,i2e),i3a=i3,i3e)
            enddo
         enddo
      enddo
#endif

      cpu=0.0d0
      cpu_ary=0.0d0

      close(10)

      end subroutine SERV_read
