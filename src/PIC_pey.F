
c THIS SUBROUTINE EXCHANGES PARTICLES IN Y-DIRECTION.
c THE ROUTINE IS CALLED BY PIC_move_part.f. PIC_pex.f, 
c PIC_pey.f, and PIC_pez.f CAN BE INTERCHANGED.

!     ------------------------------------------------------------
      module PIC_pey
!     ------------------------------------------------------------

      integer :: re2i,re4i
      real(kind=8), dimension(:), allocatable :: p_re2i
      real(kind=8), dimension(:), allocatable :: p_re4i

      end module PIC_pey

!     ------------------------------------------------------------
      subroutine PIC_pey_a(niloc, p_niloc, niloc_n)
!     ------------------------------------------------------------

      use PIC_pey
      use VLA_variables
      use PSC_domain
      use PSC_patch

      implicit none
      include './mpif.h'

      integer      :: niloc, niloc_n
      real(kind=8) :: p_niloc(0:*)

      character*(5) :: node

      integer :: lpi, e2i, e4i
      integer :: nodei,nodej,nodek,pec
      integer :: mtag,status(MPI_STATUS_SIZE)

      real(kind=8), dimension(:), allocatable :: p_e2i
      real(kind=8), dimension(:), allocatable :: p_e4i
      real(kind=8) :: restj,ys
      real(kind=8) :: ybn,ybx,ypn,ypx
      real(kind=8) :: xi,yi,zi
      real(kind=8) :: pxi,pyi,pzi
      real(kind=8) :: qni,mni,cni,lni,wni


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
c            | |            (i1mx,i2mx)| |      rd grid points in
c            | |                       | |      each spatial direction
c         e4 | |           7           | | e2   are kept in excess.
c            | |                       | |
c            | |                       | |
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

      restj=nodej/2.0-int(nodej/2.0)     ! restj=0.5 => nodej odd


c CALCULATION OF ARRAY SIZE FOR DATA EXCHANGE

! modification of simbox size for particles in y - added by ab
! particles shall not enter pml region

      ybn=i2mn*dy
      ybx=i2mx*dy
      ypn=i2n*dy
      ypx=i2x*dy

      if (nodej.eq.1.and.(boundary_fld_y1.eq.2.or.
     &    boundary_fld_y1.eq.3)) then
         ybn=(i2mn+size+1)*dy
         ypn=(i2n+size+1)*dy
      end if
      if (nodej.eq.ynpe.and.(boundary_fld_y2.eq.2.or.
     &    boundary_fld_y2.eq.3)) then
         ybx=(i2mx-size-1)*dy
         ypx=(i2x-size-1)*dy
      end if


      e2i=0
      e4i=0

      if (niloc.gt.0) then
         do l=1,niloc

            yi=p_niloc(11*l+1)
            pyi=p_niloc(11*l+4)

            if (yi<ybn) then
               e4i=e4i+1
            else if (yi>ybx) then
               e2i=e2i+1
            endif

         enddo
      endif


c LOADING MESSAGE PASSING ARRAYS


      allocate(p_e2i(0:11*e2i+10))
      allocate(p_e4i(0:11*e4i+10))

      lpi=0
      e2i=0
      e4i=0

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

            if (yi<ybn) then
               e4i=e4i+1
               p_e4i(11*e4i)=xi
               p_e4i(11*e4i+1)=yi
               p_e4i(11*e4i+2)=zi
               p_e4i(11*e4i+3)=pxi
               p_e4i(11*e4i+4)=pyi
               p_e4i(11*e4i+5)=pzi
               p_e4i(11*e4i+6)=qni
               p_e4i(11*e4i+7)=mni
               p_e4i(11*e4i+8)=cni
               p_e4i(11*e4i+9)=lni
               p_e4i(11*e4i+10)=wni
            else if (yi>ybx) then
               e2i=e2i+1
               p_e2i(11*e2i)=xi
               p_e2i(11*e2i+1)=yi
               p_e2i(11*e2i+2)=zi
               p_e2i(11*e2i+3)=pxi
               p_e2i(11*e2i+4)=pyi
               p_e2i(11*e2i+5)=pzi
               p_e2i(11*e2i+6)=qni
               p_e2i(11*e2i+7)=mni
               p_e2i(11*e2i+8)=cni
               p_e2i(11*e2i+9)=lni
               p_e2i(11*e2i+10)=wni
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


      if (nodej.lt.ynpe.and.restj<0.25) then               ! UPDATING e4
         pec=seg_inv(nodei,nodej+1,nodek)
         call MPI_SSEND(e2i,1,MPI_INTEGER,
     &                  pec,mtag,MPI_COMM_WORLD,info)
      endif
      if (1.lt.nodej.and.restj>0.25) then
         pec=seg_inv(nodei,nodej-1,nodek)
         call MPI_RECV(re4i,1,MPI_INTEGER,
     &                 pec,mtag,MPI_COMM_WORLD,status,info)
      endif

      if (nodej.lt.ynpe.and.restj>0.25) then
         pec=seg_inv(nodei,nodej+1,nodek)
         call MPI_SSEND(e2i,1,MPI_INTEGER,
     &                  pec,mtag,MPI_COMM_WORLD,info)
      endif
      if (1.lt.nodej.and.restj<0.25) then
         pec=seg_inv(nodei,nodej-1,nodek)
         call MPI_RECV(re4i,1,MPI_INTEGER,
     &                 pec,mtag,MPI_COMM_WORLD,status,info)
      endif


      if (boundary_part_y==1) then
         if (ynpe.gt.1) then                         !PERIODIC CONTINUATION
            if (nodej.eq.ynpe) then
               pec=seg_inv(nodei,1,nodek)
               call MPI_SSEND(e2i,1,MPI_INTEGER,
     &                        pec,mtag,MPI_COMM_WORLD,info)
            endif
            if (1.eq.nodej) then
               pec=seg_inv(nodei,ynpe,nodek)
               call MPI_RECV(re4i,1,MPI_INTEGER,
     &                       pec,mtag,MPI_COMM_WORLD,status,info)
            endif
         else
            re4i=e2i
         endif
      endif


      if (1.lt.nodej.and.restj<0.25) then               ! UPDATING e2
         pec=seg_inv(nodei,nodej-1,nodek)
         call MPI_SSEND(e4i,1,MPI_INTEGER,
     &                  pec,mtag,MPI_COMM_WORLD,info)
      endif
      if (nodej.lt.ynpe.and.restj>0.25) then
         pec=seg_inv(nodei,nodej+1,nodek)
         call MPI_RECV(re2i,1,MPI_INTEGER,
     &                 pec,mtag,MPI_COMM_WORLD,status,info)
      endif

      if (1.lt.nodej.and.restj>0.25) then
         pec=seg_inv(nodei,nodej-1,nodek)
         call MPI_SSEND(e4i,1,MPI_INTEGER,
     &                  pec,mtag,MPI_COMM_WORLD,info)
      endif
      if (nodej.lt.ynpe.and.restj<0.25) then
         pec=seg_inv(nodei,nodej+1,nodek)
         call MPI_RECV(re2i,1,MPI_INTEGER,
     &                 pec,mtag,MPI_COMM_WORLD,status,info)
      endif


      if (boundary_part_y==1) then
         if (ynpe.gt.1) then                         !PERIODIC CONTINUATION
            if (1.eq.nodej) then
               pec=seg_inv(nodei,ynpe,nodek)
               call MPI_SSEND(e4i,1,MPI_INTEGER,
     &                        pec,mtag,MPI_COMM_WORLD,info)
            endif
            if (nodej.eq.ynpe) then
               pec=seg_inv(nodei,1,nodek)
               call MPI_RECV(re2i,1,MPI_INTEGER,
     &                       pec,mtag,MPI_COMM_WORLD,status,info)
            endif
         else
            re2i=e4i
         endif
      endif


c PARTICLE EXCHANGE


      allocate(p_re2i(0:11*re2i+10))
      allocate(p_re4i(0:11*re4i+10))


c      ys=(i2x-i2n+1)*dy
c      ys=(i2x-i2n)*dy
      ys=ypx-ypn+dy


      if (nodej.lt.ynpe.and.restj<0.25) then               ! UPDATING e4
         pec=seg_inv(nodei,nodej+1,nodek)
         call MPI_SSEND(p_e2i,11*e2i+11,MPI_DOUBLE_PRECISION,
     &                  pec,mtag,MPI_COMM_WORLD,info)
      endif
      if (1.lt.nodej.and.restj>0.25) then
         pec=seg_inv(nodei,nodej-1,nodek)
         call MPI_RECV(p_re4i,11*re4i+11,MPI_DOUBLE_PRECISION,
     &                 pec,mtag,MPI_COMM_WORLD,status,info)
      endif

      if (nodej.lt.ynpe.and.restj>0.25) then
         pec=seg_inv(nodei,nodej+1,nodek)
         call MPI_SSEND(p_e2i,11*e2i+11,MPI_DOUBLE_PRECISION,
     &                  pec,mtag,MPI_COMM_WORLD,info)
      endif
      if (1.lt.nodej.and.restj<0.25) then
         pec=seg_inv(nodei,nodej-1,nodek)
         call MPI_RECV(p_re4i,11*re4i+11,MPI_DOUBLE_PRECISION,
     &                 pec,mtag,MPI_COMM_WORLD,status,info)
      endif


      if (boundary_part_y==1) then
         if (ynpe.gt.1) then                       !PERIODIC CONTINUATION
            if (nodej.eq.ynpe) then
               pec=seg_inv(nodei,1,nodek)
               call MPI_SSEND(p_e2i,11*e2i+11,MPI_DOUBLE_PRECISION,
     &                        pec,mtag,MPI_COMM_WORLD,info)
            endif
            if (1.eq.nodej) then
               pec=seg_inv(nodei,ynpe,nodek)
               call MPI_RECV(p_re4i,11*re4i+11,MPI_DOUBLE_PRECISION,
     &                       pec,mtag,MPI_COMM_WORLD,status,info)
               do l=1,re4i
                  p_re4i(11*l+1)=p_re4i(11*l+1)-ys
               enddo
            endif
         else
            do l=1,re4i
               p_re4i(11*l)=p_e2i(11*l)
               p_re4i(11*l+1)=p_e2i(11*l+1)-ys
               p_re4i(11*l+2)=p_e2i(11*l+2)
               p_re4i(11*l+3)=p_e2i(11*l+3)
               p_re4i(11*l+4)=p_e2i(11*l+4)
               p_re4i(11*l+5)=p_e2i(11*l+5)
               p_re4i(11*l+6)=p_e2i(11*l+6)
               p_re4i(11*l+7)=p_e2i(11*l+7)
               p_re4i(11*l+8)=p_e2i(11*l+8)
               p_re4i(11*l+9)=p_e2i(11*l+9)
               p_re4i(11*l+10)=p_e2i(11*l+10)
            enddo
         endif
      endif


      if (boundary_part_y==0) then
         if (ynpe.gt.1) then                       !MIRROR REFLECTION
            if (nodej.eq.ynpe) then
               do l=1,e2i
                  lpi=lpi+1
                  p_niloc(11*lpi)=p_e2i(11*l)
                  p_niloc(11*lpi+1)=2.0*ybx-p_e2i(11*l+1)
                  p_niloc(11*lpi+2)=p_e2i(11*l+2)
                  p_niloc(11*lpi+3)=p_e2i(11*l+3)
                  p_niloc(11*lpi+4)=-p_e2i(11*l+4)
                  p_niloc(11*lpi+5)=p_e2i(11*l+5)
                  p_niloc(11*lpi+6)=p_e2i(11*l+6)
                  p_niloc(11*lpi+7)=p_e2i(11*l+7)
                  p_niloc(11*lpi+8)=p_e2i(11*l+8)
                  p_niloc(11*lpi+9)=p_e2i(11*l+9)
                  p_niloc(11*lpi+10)=p_e2i(11*l+10)
               enddo
            endif
         else
            do l=1,e2i
               lpi=lpi+1
               p_niloc(11*lpi)=p_e2i(11*l)
               p_niloc(11*lpi+1)=2.0*ybx-p_e2i(11*l+1)
               p_niloc(11*lpi+2)=p_e2i(11*l+2)
               p_niloc(11*lpi+3)=p_e2i(11*l+3)
               p_niloc(11*lpi+4)=-p_e2i(11*l+4)
               p_niloc(11*lpi+5)=p_e2i(11*l+5)
               p_niloc(11*lpi+6)=p_e2i(11*l+6)
               p_niloc(11*lpi+7)=p_e2i(11*l+7)
               p_niloc(11*lpi+8)=p_e2i(11*l+8)
               p_niloc(11*lpi+9)=p_e2i(11*l+9)
               p_niloc(11*lpi+10)=p_e2i(11*l+10)
            enddo
         endif
      endif


      if (1.lt.nodej.and.restj<0.25) then               ! UPDATING e2
         pec=seg_inv(nodei,nodej-1,nodek)
         call MPI_SSEND(p_e4i,11*e4i+11,MPI_DOUBLE_PRECISION,
     &                  pec,mtag,MPI_COMM_WORLD,info)
      endif
      if (nodej.lt.ynpe.and.restj>0.25) then
         pec=seg_inv(nodei,nodej+1,nodek)
         call MPI_RECV(p_re2i,11*re2i+11,MPI_DOUBLE_PRECISION,
     &                 pec,mtag,MPI_COMM_WORLD,status,info)
      endif

      if (1.lt.nodej.and.restj>0.25) then
         pec=seg_inv(nodei,nodej-1,nodek)
         call MPI_SSEND(p_e4i,11*e4i+11,MPI_DOUBLE_PRECISION,
     &                  pec,mtag,MPI_COMM_WORLD,info)
      endif
      if (nodej.lt.ynpe.and.restj<0.25) then
         pec=seg_inv(nodei,nodej+1,nodek)
         call MPI_RECV(p_re2i,11*re2i+11,MPI_DOUBLE_PRECISION,
     &                 pec,mtag,MPI_COMM_WORLD,status,info)
      endif


      if (boundary_part_y==1) then 
         if (ynpe.gt.1) then                       !PERIODIC CONTINUATION
            if (1.eq.nodej) then
               pec=seg_inv(nodei,ynpe,nodek)
               call MPI_SSEND(p_e4i,11*e4i+11,MPI_DOUBLE_PRECISION,
     &                        pec,mtag,MPI_COMM_WORLD,info)
            endif
            if (nodej.eq.ynpe) then
               pec=seg_inv(nodei,1,nodek)
               call MPI_RECV(p_re2i,11*re2i+11,MPI_DOUBLE_PRECISION,
     &                       pec,mtag,MPI_COMM_WORLD,status,info)
               do l=1,re2i
                  p_re2i(11*l+1)=p_re2i(11*l+1)+ys
               enddo
            endif
         else
            do l=1,re2i
               p_re2i(11*l)=p_e4i(11*l)
               p_re2i(11*l+1)=p_e4i(11*l+1)+ys
               p_re2i(11*l+2)=p_e4i(11*l+2)
               p_re2i(11*l+3)=p_e4i(11*l+3)
               p_re2i(11*l+4)=p_e4i(11*l+4)
               p_re2i(11*l+5)=p_e4i(11*l+5)
               p_re2i(11*l+6)=p_e4i(11*l+6)
               p_re2i(11*l+7)=p_e4i(11*l+7)
               p_re2i(11*l+8)=p_e4i(11*l+8)
               p_re2i(11*l+9)=p_e4i(11*l+9)
               p_re2i(11*l+10)=p_e4i(11*l+10)
            enddo
         endif
      endif


      if (boundary_part_y==0) then
         if (ynpe.gt.1) then                       !MIRROR REFLECTION
            if (1.eq.nodej) then
               do l=1,e4i
                  lpi=lpi+1
                  p_niloc(11*lpi)=p_e4i(11*l)
                  p_niloc(11*lpi+1)=2.0*ybn-p_e4i(11*l+1)
                  p_niloc(11*lpi+2)=p_e4i(11*l+2)
                  p_niloc(11*lpi+3)=p_e4i(11*l+3)
                  p_niloc(11*lpi+4)=-p_e4i(11*l+4)
                  p_niloc(11*lpi+5)=p_e4i(11*l+5)
                  p_niloc(11*lpi+6)=p_e4i(11*l+6)
                  p_niloc(11*lpi+7)=p_e4i(11*l+7)
                  p_niloc(11*lpi+8)=p_e4i(11*l+8)
                  p_niloc(11*lpi+9)=p_e4i(11*l+9)
                  p_niloc(11*lpi+10)=p_e4i(11*l+10)
               enddo
            endif
         else
            do l=1,e4i
               lpi=lpi+1
               p_niloc(11*lpi)=p_e4i(11*l)
               p_niloc(11*lpi+1)=2.0*ybn-p_e4i(11*l+1)
               p_niloc(11*lpi+2)=p_e4i(11*l+2)
               p_niloc(11*lpi+3)=p_e4i(11*l+3)
               p_niloc(11*lpi+4)=-p_e4i(11*l+4)
               p_niloc(11*lpi+5)=p_e4i(11*l+5)
               p_niloc(11*lpi+6)=p_e4i(11*l+6)
               p_niloc(11*lpi+7)=p_e4i(11*l+7)
               p_niloc(11*lpi+8)=p_e4i(11*l+8)
               p_niloc(11*lpi+9)=p_e4i(11*l+9)
               p_niloc(11*lpi+10)=p_e4i(11*l+10)
            enddo
         endif
      endif


c LOCAL PARTICLE ARRAY UPDATE


      deallocate(p_e2i)
      deallocate(p_e4i)

      niloc=lpi
      niloc_n=lpi+re2i+re4i

      end subroutine PIC_pey_a


!     ------------------------------------------------------------
      subroutine PIC_pey_b(niloc, p_niloc)
!     ------------------------------------------------------------

      use PIC_pey
      use PSC_patch

      implicit none
      include './mpif.h'

      integer      :: niloc
      real(kind=8) :: p_niloc(0:*)

      integer      :: l

      if (re2i>0) then
         do l=1,re2i
            niloc=niloc+1
            p_niloc(11*niloc)=p_re2i(11*l)
            p_niloc(11*niloc+1)=p_re2i(11*l+1)
            p_niloc(11*niloc+2)=p_re2i(11*l+2)
            p_niloc(11*niloc+3)=p_re2i(11*l+3)
            p_niloc(11*niloc+4)=p_re2i(11*l+4)
            p_niloc(11*niloc+5)=p_re2i(11*l+5)
            p_niloc(11*niloc+6)=p_re2i(11*l+6)
            p_niloc(11*niloc+7)=p_re2i(11*l+7)
            p_niloc(11*niloc+8)=p_re2i(11*l+8)
            p_niloc(11*niloc+9)=p_re2i(11*l+9)
            p_niloc(11*niloc+10)=p_re2i(11*l+10)
         enddo
      endif

      if (re4i>0) then
         do l=1,re4i
            niloc=niloc+1
            p_niloc(11*niloc)=p_re4i(11*l)
            p_niloc(11*niloc+1)=p_re4i(11*l+1)
            p_niloc(11*niloc+2)=p_re4i(11*l+2)
            p_niloc(11*niloc+3)=p_re4i(11*l+3)
            p_niloc(11*niloc+4)=p_re4i(11*l+4)
            p_niloc(11*niloc+5)=p_re4i(11*l+5)
            p_niloc(11*niloc+6)=p_re4i(11*l+6)
            p_niloc(11*niloc+7)=p_re4i(11*l+7)
            p_niloc(11*niloc+8)=p_re4i(11*l+8)
            p_niloc(11*niloc+9)=p_re4i(11*l+9)
            p_niloc(11*niloc+10)=p_re4i(11*l+10)
         enddo
      endif


      deallocate(p_re2i)
      deallocate(p_re4i)

      end subroutine PIC_pey_b
