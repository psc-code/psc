c THIS SUBROUTINE IS CALLED BY INIT_field.f AND
c PIC_move_part.f. IT CORRECTS FIELDS DERIVED FROM
c PARTICLES. PIC_fax.f, PIC_fay.f, and PIC_faz.f
c CAN BE FREELY INTERCHANGED.


      subroutine PIC_fay(fd)

      use PIC_variables
      use VLA_variables

      implicit none
      include './mpif.h'

      integer rysize
      integer nodei,nodej,nodek
      integer mtag,status(MPI_STATUS_SIZE)
      
      real(kind=8) :: restj
      real(kind=8) :: fd(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &                   i3mn-rd3:i3mx+rd3)

      real(kind=8),allocatable,dimension(:,:,:) :: rimy


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
c
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


c initialization


      mtag=300

      nodei=seg_i1(mpe)
      nodej=seg_i2(mpe)
      nodek=seg_i3(mpe)

      restj=nodej/2.0-int(nodej/2.0)                       ! restj=0.5 => nodej odd

      rysize=(i1mx-i1mn+2*rd1+1)*rd2*(i3mx-i3mn+2*rd2+1)

      allocate(rimy(i1mn-rd1:i1mx+rd1,1:rd2,i3mn-rd3:i3mx+rd3))


c*************************************************************************
c e2 addition
c*************************************************************************


      if (nodej.lt.ynpe.and.restj<0.25) then
         pec=seg_inv(nodei,nodej+1,nodek)
         do i3=i3mn-rd3,i3mx+rd3
            do i2=1,rd2 
               do i1=i1mn-rd1,i1mx+rd1 
                  rimy(i1,i2,i3)=fd(i1,i2mx+i2,i3)
               enddo
            enddo
         enddo
         call MPI_SSEND(rimy,rysize,MPI_DOUBLE_PRECISION,
     &                  pec,mtag,MPI_COMM_WORLD,info)
      endif
      if (1.lt.nodej.and.restj>0.25) then
         pec=seg_inv(nodei,nodej-1,nodek)
         call MPI_RECV(rimy,rysize,MPI_DOUBLE_PRECISION,
     &                 pec,mtag,MPI_COMM_WORLD,status,info)
         do i3=i3mn-rd3,i3mx+rd3
            do i2=1,rd2 
               do i1=i1mn-rd1,i1mx+rd1
                  fd(i1,i2mn+i2-1,i3)=fd(i1,i2mn+i2-1,i3)
     &                                +rimy(i1,i2,i3)
               enddo
            enddo
         enddo
      endif

      if (nodej.lt.ynpe.and.restj>0.25) then
         pec=seg_inv(nodei,nodej+1,nodek)
         do i3=i3mn-rd3,i3mx+rd3
            do i2=1,rd2 
               do i1=i1mn-rd1,i1mx+rd1 
                  rimy(i1,i2,i3)=fd(i1,i2mx+i2,i3)
               enddo
            enddo
         enddo
         call MPI_SSEND(rimy,rysize,MPI_DOUBLE_PRECISION,
     &                  pec,mtag,MPI_COMM_WORLD,info)
      endif
      if (1.lt.nodej.and.restj<0.25) then
         pec=seg_inv(nodei,nodej-1,nodek)
         call MPI_RECV(rimy,rysize,MPI_DOUBLE_PRECISION,
     &                 pec,mtag,MPI_COMM_WORLD,status,info)
         do i3=i3mn-rd3,i3mx+rd3
            do i2=1,rd2
               do i1=i1mn-rd1,i1mx+rd1
                  fd(i1,i2mn+i2-1,i3)=fd(i1,i2mn+i2-1,i3)
     &                                +rimy(i1,i2,i3)
               enddo
            enddo
         enddo
      endif

c*************************************************************************
c e2 boundary addition
c*************************************************************************

      if (ynpe.gt.1) then
         if (nodej.eq.ynpe) then
            pec=seg_inv(nodei,1,nodek)
            do i3=i3mn-rd3,i3mx+rd3
               do i2=1,rd2
                  do i1=i1mn-rd1,i1mx+rd1 
                     rimy(i1,i2,i3)=fd(i1,i2mx+i2,i3)
                  enddo
               enddo
            enddo
            call MPI_SSEND(rimy,rysize,MPI_DOUBLE_PRECISION,
     &                     pec,mtag,MPI_COMM_WORLD,info)
         endif
         if (1.eq.nodej) then
            pec=seg_inv(nodei,ynpe,nodek)
            call MPI_RECV(rimy,rysize,MPI_DOUBLE_PRECISION,
     &                    pec,mtag,MPI_COMM_WORLD,status,info)
            do i3=i3mn-rd3,i3mx+rd3
               do i2=1,rd2 
                  do i1=i1mn-rd1,i1mx+rd1
                     fd(i1,i2mn+i2-1,i3)=fd(i1,i2mn+i2-1,i3)
     &                                   +rimy(i1,i2,i3)
                  enddo
               enddo
            enddo
         endif
      else
         if (i2n.ne.i2x) then
            do i3=i3mn-rd3,i3mx+rd3
               do i2=1,rd2
                  do i1=i1mn-rd1,i1mx+rd1
                     rimy(i1,i2,i3)=fd(i1,i2mx+i2,i3)
                     fd(i1,i2mn+i2-1,i3)=fd(i1,i2mn+i2-1,i3)
     &                                   +rimy(i1,i2,i3)
                  enddo
               enddo
            enddo
         else
            do i3=i3mn-rd3,i3mx+rd3
               do i2=1,rd2
                  do i1=i1mn-rd1,i1mx+rd1
                     rimy(i1,i2,i3)=fd(i1,i2mx+i2,i3)
                     fd(i1,i2mn,i3)=fd(i1,i2mn,i3)+rimy(i1,i2,i3)
                  enddo
               enddo
            enddo
         endif
      endif

c*************************************************************************
c e4 addition
c*************************************************************************

      if (1.lt.nodej.and.restj<0.25) then
         pec=seg_inv(nodei,nodej-1,nodek)
         do i3=i3mn-rd3,i3mx+rd3
            do i2=1,rd2
               do i1=i1mn-rd1,i1mx+rd1 
                  rimy(i1,i2,i3)=fd(i1,i2mn-i2,i3)
               enddo
            enddo
         enddo
         call MPI_SSEND(rimy,rysize,MPI_DOUBLE_PRECISION,
     &                  pec,mtag,MPI_COMM_WORLD,info)
      endif
      if (nodej.lt.ynpe.and.restj>0.25) then
         pec=seg_inv(nodei,nodej+1,nodek)
         call MPI_RECV(rimy,rysize,MPI_DOUBLE_PRECISION,
     &                 pec,mtag,MPI_COMM_WORLD,status,info)
         do i3=i3mn-rd3,i3mx+rd3
            do i2=1,rd2
               do i1=i1mn-rd1,i1mx+rd1
                  fd(i1,i2mx-i2+1,i3)=fd(i1,i2mx-i2+1,i3)
     &                                +rimy(i1,i2,i3)
               enddo
            enddo
         enddo
      endif

      if (1.lt.nodej.and.restj>0.25) then
         pec=seg_inv(nodei,nodej-1,nodek)
         do i3=i3mn-rd3,i3mx+rd3
            do i2=1,rd2
               do i1=i1mn-rd1,i1mx+rd1
                  rimy(i1,i2,i3)=fd(i1,i2mn-i2,i3)
               enddo
            enddo
         enddo
         call MPI_SSEND(rimy,rysize,MPI_DOUBLE_PRECISION,
     &                  pec,mtag,MPI_COMM_WORLD,info)
      endif
      if (nodej.lt.ynpe.and.restj<0.25) then
         pec=seg_inv(nodei,nodej+1,nodek)
         call MPI_RECV(rimy,rysize,MPI_DOUBLE_PRECISION,
     &                 pec,mtag,MPI_COMM_WORLD,status,info)
         do i3=i3mn-rd3,i3mx+rd3
            do i2=1,rd2
               do i1=i1mn-rd1,i1mx+rd1
                  fd(i1,i2mx-i2+1,i3)=fd(i1,i2mx-i2+1,i3)
     &                                +rimy(i1,i2,i3)
               enddo
            enddo
         enddo
      endif

c*************************************************************************
c e4 boundary addition
c*************************************************************************

      if (ynpe.gt.1) then
         if (1.eq.nodej) then
            pec=seg_inv(nodei,ynpe,nodek)
            do i3=i3mn-rd3,i3mx+rd3
               do i2=1,rd2
                  do i1=i1mn-rd1,i1mx+rd1 
                     rimy(i1,i2,i3)=fd(i1,i2mn-i2,i3)
                  enddo
               enddo
            enddo
            call MPI_SSEND(rimy,rysize,MPI_DOUBLE_PRECISION,
     &                     pec,mtag,MPI_COMM_WORLD,info)
         endif
         if (nodej.eq.ynpe) then
            pec=seg_inv(nodei,1,nodek)
            call MPI_RECV(rimy,rysize,MPI_DOUBLE_PRECISION,
     &                    pec,mtag,MPI_COMM_WORLD,status,info)
            do i3=i3mn-rd3,i3mx+rd3
               do i2=1,rd2
                  do i1=i1mn-rd1,i1mx+rd1
                     fd(i1,i2mx-i2+1,i3)=fd(i1,i2mx-i2+1,i3)
     &                                   +rimy(i1,i2,i3)
                  enddo
               enddo
            enddo
         endif
      else
         if (i2n.ne.i2x) then
            do i3=i3mn-rd3,i3mx+rd3
               do i2=1,rd2
                  do i1=i1mn-rd1,i1mx+rd1 
                     rimy(i1,i2,i3)=fd(i1,i2mn-i2,i3)
                     fd(i1,i2mx-i2+1,i3)=fd(i1,i2mx-i2+1,i3)
     &                                   +rimy(i1,i2,i3)
                  enddo
               enddo
            enddo
         else
            do i3=i3mn-rd3,i3mx+rd3
               do i2=1,rd2
                  do i1=i1mn-rd1,i1mx+rd1 
                     rimy(i1,i2,i3)=fd(i1,i2mn-i2,i3)
                     fd(i1,i2mx,i3)=fd(i1,i2mx,i3)+rimy(i1,i2,i3)
                  enddo
               enddo
            enddo
         endif
      endif
 

      deallocate(rimy)

 
      end subroutine PIC_fay
