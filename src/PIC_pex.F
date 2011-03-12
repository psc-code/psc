
c THIS SUBROUTINE EXCHANGES PARTICLES IN X-DIRECTION.
c THE ROUTINE IS CALLED BY PIC_move_part.f PIC_pex.f, 
c PIC_pey.f, and PIC_pez.f CAN BE INTERCHANGED.

!     ------------------------------------------------------------
      subroutine PIC_pex_a(niloc, p_niloc, niloc_n)
!     ------------------------------------------------------------

      use PIC_variables, only: p_e1i, p_e3i, p_re1i, p_re3i, e1i,
     *     e3i, lpi, re1i, re3i
      use VLA_variables
      use PSC_patch

      implicit none
      include './mpif.h'

      integer      :: niloc, niloc_n
      real(kind=8) :: p_niloc(0:*)

      character*(5) :: node

      integer :: la,le
      integer :: nodei,nodej,nodek
      integer :: mtag,status(MPI_STATUS_SIZE)

      real(kind=8) :: resti,xs
      real(kind=8) :: xbn,xbx,xpn,xpx
      real(kind=8) :: xi,yi,zi
      real(kind=8) :: pxi,pyi,pzi
      real(kind=8) :: qni,mni,cni,lni,wni

      real(kind=8) :: s_cpub  


      s_cpub=0.0d0


c---------------------------------------------------------------------
c TOPOLOGY AND CONVENTIONS (example of 12 nodes)
c---------------------------------------------------------------------
c  topology:  npe=12
c
c             -------------------------
c             |  2  |  5  |  8  | 11  |
c             -------------------------      
c  x, xnpe=3  |  1  |  4  |  7  | 10  |      0.le.mpe.le.npe-1
c             -------------------------
c             |  0  |  3  |  6  |  9  |
c             -------------------------
c                     y, ynpe=4
c
c
c  transcription:     x, xnpe=4
c
c             -------------------------
c             | 31  | 32  | 33  | 34  |      nodei=seg_i1(mpe)
c             -------------------------      nodej=seg_i2(mpe)      
c  x, xnpe=3  | 21  | 22  | 23  | 24  |
c             -------------------------      1.le.nodei.le.xnpe
c             | 11  | 12  | 13  | 14  |      1.le.nodej.le.ynpe
c             -------------------------
c                     y, ynpe=4
c
C
c  memory on node 7 = node 23:
c
c                         e3              
c  i1mx+rd   -----------------------------
c            | ------------------------- |     
c            | |            (i1mx,i2mx)| |      rdi grid points in
c            | |                       | |      each spatial direction
c         e4 | |           7           | | e2   are kept in excess. e5
c            | |                       | |      is in the back and e6
c            | |                       | |      in the front.
c            | |(i1mn,i2mn)            | |
c            | ------------------------- |
c  i1mn-rd   -----------------------------
c                         e1              
c          i2mn-rd                   i2mx+rd
c
c         rd: width of additional data space
c      e1-e4: edge regions of the grid
c
c---------------------------------------------------------------------


c SETUP


      mtag=200

      nodei=seg_i1(mpe)
      nodej=seg_i2(mpe)
      nodek=seg_i3(mpe)

      resti=nodei/2.0-int(nodei/2.0)     ! resti=0.5 => nodei odd


c CALCULATION OF ARRAY SIZE FOR DATA EXCHANGE

! modification of simbox size for particles in y - added by ab
! particles shall not enter pml region


      xbn=i1mn*dx
      xbx=i1mx*dx
      xpn=i1n*dx
      xpx=i1x*dx


      if (nodei.eq.1.and.(boundary_fld_x1.eq.2.or.
     &    boundary_fld_x1.eq.3)) then
         xbn=(i1mn+size+1)*dx
         xpn=(i1n+size+1)*dx
      end if
      if (nodei.eq.xnpe.and.(boundary_fld_x2.eq.2.or.
     &    boundary_fld_x2.eq.3)) then
         xbx=(i1mx-size-1)*dx
         xpx=(i1x+size+1)*dx
      end if
  

      e1i=0
      e3i=0
      lpi=0

      if (niloc.gt.0) then
         do l=1,niloc

            xi=p_niloc(11*l)
            pxi=p_niloc(11*l+3)

            if (xi<xbn) then
               e1i=e1i+1
            else if (xi>xbx) then
               e3i=e3i+1
            else
               lpi=lpi+1
            endif

         enddo
      endif


c LOADING MESSAGE PASSING ARRAYS


      allocate(p_e1i(0:11*e1i+10))
      allocate(p_e3i(0:11*e3i+10))

      lpi=0
      e1i=0
      e3i=0

      if (niloc.gt.0) then
         do l=1,niloc

            xi=p_niloc(11*l)
            yi=p_niloc(11*l+1)
            zi=p_niloc(11*l+2)
            pxi=p_niloc(11*l+3)
            pyi=p_niloc(11*l+4)
            pzi=p_niloc(11*l+5)
            qni=p_niloc(11*l+6)
            mni=p_niloc(11*l+7)
            cni=p_niloc(11*l+8)
            lni=p_niloc(11*l+9)
            wni=p_niloc(11*l+10)

            if (xi<xbn) then
               e1i=e1i+1
               p_e1i(11*e1i)=xi
               p_e1i(11*e1i+1)=yi
               p_e1i(11*e1i+2)=zi
               p_e1i(11*e1i+3)=pxi
               p_e1i(11*e1i+4)=pyi
               p_e1i(11*e1i+5)=pzi
               p_e1i(11*e1i+6)=qni
               p_e1i(11*e1i+7)=mni
               p_e1i(11*e1i+8)=cni
               p_e1i(11*e1i+9)=lni
               p_e1i(11*e1i+10)=wni
            else if (xi>xbx) then
               e3i=e3i+1
               p_e3i(11*e3i)=xi
               p_e3i(11*e3i+1)=yi
               p_e3i(11*e3i+2)=zi
               p_e3i(11*e3i+3)=pxi
               p_e3i(11*e3i+4)=pyi
               p_e3i(11*e3i+5)=pzi
               p_e3i(11*e3i+6)=qni
               p_e3i(11*e3i+7)=mni
               p_e3i(11*e3i+8)=cni
               p_e3i(11*e3i+9)=lni
               p_e3i(11*e3i+10)=wni
            else
               lpi=lpi+1
               p_niloc(11*lpi)=xi
               p_niloc(11*lpi+1)=yi
               p_niloc(11*lpi+2)=zi
               p_niloc(11*lpi+3)=pxi
               p_niloc(11*lpi+4)=pyi
               p_niloc(11*lpi+5)=pzi
               p_niloc(11*lpi+6)=qni
               p_niloc(11*lpi+7)=mni
               p_niloc(11*lpi+8)=cni
               p_niloc(11*lpi+9)=lni
               p_niloc(11*lpi+10)=wni
            endif

         enddo
      endif


C ARRAY DIMENSION EXCHANGE


      call SERV_systime(cpua)
      if (nodei.lt.xnpe.and.resti<0.25) then               ! UPDATING e1
         pec=seg_inv(nodei+1,nodej,nodek)
         call MPI_SSEND(e3i,1,MPI_INTEGER,
     &                  pec,mtag,MPI_COMM_WORLD,info)
      endif
      if (1.lt.nodei.and.resti>0.25) then
         pec=seg_inv(nodei-1,nodej,nodek)
         call MPI_RECV(re1i,1,MPI_INTEGER,
     &                 pec,mtag,MPI_COMM_WORLD,status,info)
      endif

      if (nodei.lt.xnpe.and.resti>0.25) then
         pec=seg_inv(nodei+1,nodej,nodek)
         call MPI_SSEND(e3i,1,MPI_INTEGER,
     &                  pec,mtag,MPI_COMM_WORLD,info)
      endif
      if (1.lt.nodei.and.resti<0.25) then
         pec=seg_inv(nodei-1,nodej,nodek)
         call MPI_RECV(re1i,1,MPI_INTEGER,
     &                 pec,mtag,MPI_COMM_WORLD,status,info)
      endif


      if (boundary_part_x==1) then
         if (xnpe.gt.1) then                         !PERIODIC CONTINUATION
            if (nodei.eq.xnpe) then
               pec=seg_inv(1,nodej,nodek)
               call MPI_SSEND(e3i,1,MPI_INTEGER,
     &                        pec,mtag,MPI_COMM_WORLD,info)
            endif
            if (1.eq.nodei) then
               pec=seg_inv(xnpe,nodej,nodek)
               call MPI_RECV(re1i,1,MPI_INTEGER,
     &                       pec,mtag,MPI_COMM_WORLD,status,info)
           endif
         else
            re1i=e3i
         endif
      endif


      if (1.lt.nodei.and.resti<0.25) then               ! UPDATING e3
         pec=seg_inv(nodei-1,nodej,nodek)
         call MPI_SSEND(e1i,1,MPI_INTEGER,
     &                  pec,mtag,MPI_COMM_WORLD,info)
      endif
      if (nodei.lt.xnpe.and.resti>0.25) then
         pec=seg_inv(nodei+1,nodej,nodek)
         call MPI_RECV(re3i,1,MPI_INTEGER,
     &                 pec,mtag,MPI_COMM_WORLD,status,info)
      endif

      if (1.lt.nodei.and.resti>0.25) then
         pec=seg_inv(nodei-1,nodej,nodek)
         call MPI_SSEND(e1i,1,MPI_INTEGER,
     &                  pec,mtag,MPI_COMM_WORLD,info)
      endif
      if (nodei.lt.xnpe.and.resti<0.25) then
         pec=seg_inv(nodei+1,nodej,nodek)
         call MPI_RECV(re3i,1,MPI_INTEGER,
     &                 pec,mtag,MPI_COMM_WORLD,status,info)
      endif


      if (boundary_part_x==1) then 
         if (xnpe.gt.1) then                         !PERIODIC CONTINUATION
            if (1.eq.nodei) then
               pec=seg_inv(xnpe,nodej,nodek)
               call MPI_SSEND(e1i,1,MPI_INTEGER,
     &                        pec,mtag,MPI_COMM_WORLD,info)
            endif
            if (nodei.eq.xnpe) then
               pec=seg_inv(1,nodej,nodek)
               call MPI_RECV(re3i,1,MPI_INTEGER,
     &                       pec,mtag,MPI_COMM_WORLD,status,info)
            endif
         else
            re3i=e1i
         endif
      endif


c PARTICLE EXCHANGE


      allocate(p_re1i(0:11*re1i+10))
      allocate(p_re3i(0:11*re3i+10))


c      xs=(i1x-i1n+1)*dx
c      xs=(i1x-i1n)*dx
      xs=xpx-xpn+dx


      if (nodei.lt.xnpe.and.resti<0.25) then               ! UPDATING e1
         pec=seg_inv(nodei+1,nodej,nodek)
         call MPI_SSEND(p_e3i,11*e3i+11,MPI_DOUBLE_PRECISION,
     &                  pec,mtag,MPI_COMM_WORLD,info)
      endif
      if (1.lt.nodei.and.resti>0.25) then
         pec=seg_inv(nodei-1,nodej,nodek)
         call MPI_RECV(p_re1i,11*re1i+11,MPI_DOUBLE_PRECISION,
     &                 pec,mtag,MPI_COMM_WORLD,status,info)
      endif

      if (nodei.lt.xnpe.and.resti>0.25) then
         pec=seg_inv(nodei+1,nodej,nodek)
         call MPI_SSEND(p_e3i,11*e3i+11,MPI_DOUBLE_PRECISION,
     &                  pec,mtag,MPI_COMM_WORLD,info)
      endif
      if (1.lt.nodei.and.resti<0.25) then
         pec=seg_inv(nodei-1,nodej,nodek)
         call MPI_RECV(p_re1i,11*re1i+11,MPI_DOUBLE_PRECISION,
     &                 pec,mtag,MPI_COMM_WORLD,status,info)
      endif


      if (boundary_part_x==1) then
         if (xnpe.gt.1) then                       !PERIODIC CONTINUATION
            if (nodei.eq.xnpe) then
               pec=seg_inv(1,nodej,nodek)
               call MPI_SSEND(p_e3i,11*e3i+11,MPI_DOUBLE_PRECISION,
     &                        pec,mtag,MPI_COMM_WORLD,info)
            endif
            if (1.eq.nodei) then
               pec=seg_inv(xnpe,nodej,nodek)
               call MPI_RECV(p_re1i,11*re1i+11,MPI_DOUBLE_PRECISION,
     &                       pec,mtag,MPI_COMM_WORLD,status,info)
               do l=1,re1i
                  p_re1i(11*l)=p_re1i(11*l)-xs
               enddo
            endif
         else
            do l=1,re1i
               p_re1i(11*l)=p_e3i(11*l)-xs
               p_re1i(11*l+1)=p_e3i(11*l+1)
               p_re1i(11*l+2)=p_e3i(11*l+2)
               p_re1i(11*l+3)=p_e3i(11*l+3)
               p_re1i(11*l+4)=p_e3i(11*l+4)
               p_re1i(11*l+5)=p_e3i(11*l+5)
               p_re1i(11*l+6)=p_e3i(11*l+6)
               p_re1i(11*l+7)=p_e3i(11*l+7)
               p_re1i(11*l+8)=p_e3i(11*l+8)
               p_re1i(11*l+9)=p_e3i(11*l+9)
               p_re1i(11*l+10)=p_e3i(11*l+10)
            enddo
         endif
      endif


      if (boundary_part_x==0) then
         if (xnpe.gt.1) then                       !MIRROR REFLECTION
            if (nodei.eq.xnpe) then
               do l=1,e3i
                  lpi=lpi+1
                  p_niloc(11*lpi)=2.0*xbx-p_e3i(11*l)
                  p_niloc(11*lpi+1)=p_e3i(11*l+1)
                  p_niloc(11*lpi+2)=p_e3i(11*l+2)
                  p_niloc(11*lpi+3)=-p_e3i(11*l+3)
                  p_niloc(11*lpi+4)=p_e3i(11*l+4)
                  p_niloc(11*lpi+5)=p_e3i(11*l+5)
                  p_niloc(11*lpi+6)=p_e3i(11*l+6)
                  p_niloc(11*lpi+7)=p_e3i(11*l+7)
                  p_niloc(11*lpi+8)=p_e3i(11*l+8)
                  p_niloc(11*lpi+9)=p_e3i(11*l+9)
                  p_niloc(11*lpi+10)=p_e3i(11*l+10)
               enddo
            endif
         else
            do l=1,e3i
               lpi=lpi+1
               p_niloc(11*lpi)=2.0*xbx-p_e3i(11*l)
               p_niloc(11*lpi+1)=p_e3i(11*l+1)
               p_niloc(11*lpi+2)=p_e3i(11*l+2)
               p_niloc(11*lpi+3)=-p_e3i(11*l+3)
               p_niloc(11*lpi+4)=p_e3i(11*l+4)
               p_niloc(11*lpi+5)=p_e3i(11*l+5)
               p_niloc(11*lpi+6)=p_e3i(11*l+6)
               p_niloc(11*lpi+7)=p_e3i(11*l+7)
               p_niloc(11*lpi+8)=p_e3i(11*l+8)
               p_niloc(11*lpi+9)=p_e3i(11*l+9)
               p_niloc(11*lpi+10)=p_e3i(11*l+10)
            enddo
         endif
      endif


      if (1.lt.nodei.and.resti<0.25) then               ! UPDATING e3
         pec=seg_inv(nodei-1,nodej,nodek)
         call MPI_SSEND(p_e1i,11*e1i+11,MPI_DOUBLE_PRECISION,
     &                  pec,mtag,MPI_COMM_WORLD,info)
      endif
      if (nodei.lt.xnpe.and.resti>0.25) then
         pec=seg_inv(nodei+1,nodej,nodek)
         call MPI_RECV(p_re3i,11*re3i+11,MPI_DOUBLE_PRECISION,
     &                 pec,mtag,MPI_COMM_WORLD,status,info)
      endif

      if (1.lt.nodei.and.resti>0.25) then
         pec=seg_inv(nodei-1,nodej,nodek)
         call MPI_SSEND(p_e1i,11*e1i+11,MPI_DOUBLE_PRECISION,
     &                  pec,mtag,MPI_COMM_WORLD,info)
      endif
      if (nodei.lt.xnpe.and.resti<0.25) then
         pec=seg_inv(nodei+1,nodej,nodek)
         call MPI_RECV(p_re3i,11*re3i+11,MPI_DOUBLE_PRECISION,
     &                 pec,mtag,MPI_COMM_WORLD,status,info)
      endif


      if (boundary_part_x==1) then
         if (xnpe.gt.1) then                       !PERIODIC CONTINUATION
            if (1.eq.nodei) then
               pec=seg_inv(xnpe,nodej,nodek)
               call MPI_SSEND(p_e1i,11*e1i+11,MPI_DOUBLE_PRECISION,
     &                        pec,mtag,MPI_COMM_WORLD,info)
            endif
            if (nodei.eq.xnpe) then
               pec=seg_inv(1,nodej,nodek)
               call MPI_RECV(p_re3i,11*re3i+11,MPI_DOUBLE_PRECISION,
     &                       pec,mtag,MPI_COMM_WORLD,status,info)
               do l=1,re3i
                  p_re3i(11*l)=p_re3i(11*l)+xs
               enddo
            endif
         else
            do l=1,re3i
               p_re3i(11*l)=p_e1i(11*l)+xs
               p_re3i(11*l+1)=p_e1i(11*l+1)
               p_re3i(11*l+2)=p_e1i(11*l+2)
               p_re3i(11*l+3)=p_e1i(11*l+3)
               p_re3i(11*l+4)=p_e1i(11*l+4)
               p_re3i(11*l+5)=p_e1i(11*l+5)
               p_re3i(11*l+6)=p_e1i(11*l+6)
               p_re3i(11*l+7)=p_e1i(11*l+7)
               p_re3i(11*l+8)=p_e1i(11*l+8)
               p_re3i(11*l+9)=p_e1i(11*l+9)
               p_re3i(11*l+10)=p_e1i(11*l+10)
            enddo
         endif
      endif


      if (boundary_part_x==0) then
         if (xnpe.gt.1) then                       !MIRROR REFLECTION
            if (1.eq.nodei) then
               do l=1,e1i
                  lpi=lpi+1
                  p_niloc(11*lpi)=2.0*xbn-p_e1i(11*l)
                  p_niloc(11*lpi+1)=p_e1i(11*l+1)
                  p_niloc(11*lpi+2)=p_e1i(11*l+2)
                  p_niloc(11*lpi+3)=-p_e1i(11*l+3)
                  p_niloc(11*lpi+4)=p_e1i(11*l+4)
                  p_niloc(11*lpi+5)=p_e1i(11*l+5)
                  p_niloc(11*lpi+6)=p_e1i(11*l+6)
                  p_niloc(11*lpi+7)=p_e1i(11*l+7)
                  p_niloc(11*lpi+8)=p_e1i(11*l+8)
                  p_niloc(11*lpi+9)=p_e1i(11*l+9)
                  p_niloc(11*lpi+10)=p_e1i(11*l+10)
               enddo
            endif
         else
            do l=1,e1i
               lpi=lpi+1
               p_niloc(11*lpi)=2.0*xbn-p_e1i(11*l)
               p_niloc(11*lpi+1)=p_e1i(11*l+1)
               p_niloc(11*lpi+2)=p_e1i(11*l+2)
               p_niloc(11*lpi+3)=-p_e1i(11*l+3)
               p_niloc(11*lpi+4)=p_e1i(11*l+4)
               p_niloc(11*lpi+5)=p_e1i(11*l+5)
               p_niloc(11*lpi+6)=p_e1i(11*l+6)
               p_niloc(11*lpi+7)=p_e1i(11*l+7)
               p_niloc(11*lpi+8)=p_e1i(11*l+8)
               p_niloc(11*lpi+9)=p_e1i(11*l+9)
               p_niloc(11*lpi+10)=p_e1i(11*l+10)
            enddo
         endif
      endif


      call SERV_systime(cpub)
      s_cpub=s_cpub+cpub-cpua


c LOCAL PARTICLE ARRAY UPDATE


      deallocate(p_e1i)
      deallocate(p_e3i)

      niloc_n=lpi+re1i+re3i

      cpumess=cpumess+s_cpub

      end subroutine PIC_pex_a


!     ------------------------------------------------------------
      subroutine PIC_pex_b(niloc, p_niloc)
!     ------------------------------------------------------------

      use PIC_variables, only: p_e1i, p_e3i, p_re1i, p_re3i, e1i,
     *     e3i, lpi, re1i, re3i
      use PSC_patch

      implicit none
      include './mpif.h'

      integer      :: niloc
      real(kind=8) :: p_niloc(0:*)

      integer      :: l

      niloc=lpi
      if (re1i>0) then
         do l=1,re1i
            niloc=niloc+1
            p_niloc(11*niloc)=p_re1i(11*l)
            p_niloc(11*niloc+1)=p_re1i(11*l+1)
            p_niloc(11*niloc+2)=p_re1i(11*l+2)
            p_niloc(11*niloc+3)=p_re1i(11*l+3)
            p_niloc(11*niloc+4)=p_re1i(11*l+4)
            p_niloc(11*niloc+5)=p_re1i(11*l+5)
            p_niloc(11*niloc+6)=p_re1i(11*l+6)
            p_niloc(11*niloc+7)=p_re1i(11*l+7)
            p_niloc(11*niloc+8)=p_re1i(11*l+8)
            p_niloc(11*niloc+9)=p_re1i(11*l+9)
            p_niloc(11*niloc+10)=p_re1i(11*l+10)
         enddo
      endif

      if (re3i>0) then
         do l=1,re3i
            niloc=niloc+1
            p_niloc(11*niloc)=p_re3i(11*l)
            p_niloc(11*niloc+1)=p_re3i(11*l+1)
            p_niloc(11*niloc+2)=p_re3i(11*l+2)
            p_niloc(11*niloc+3)=p_re3i(11*l+3)
            p_niloc(11*niloc+4)=p_re3i(11*l+4)
            p_niloc(11*niloc+5)=p_re3i(11*l+5)
            p_niloc(11*niloc+6)=p_re3i(11*l+6)
            p_niloc(11*niloc+7)=p_re3i(11*l+7)
            p_niloc(11*niloc+8)=p_re3i(11*l+8)
            p_niloc(11*niloc+9)=p_re3i(11*l+9)
            p_niloc(11*niloc+10)=p_re3i(11*l+10)
         enddo
      endif


      deallocate(p_re1i)
      deallocate(p_re3i)


      end subroutine PIC_pex_b
