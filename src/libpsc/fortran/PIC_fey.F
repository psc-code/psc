
c This subroutine exchanges the Maxwell fields between
c computation domains for the parallel Maxwell solver.
c Is called in PIC_msa.f and PIC_msb.f. 


      subroutine PIC_fey(fd)
      
      use VLA_variables, only: n
      use PSC_domain
      use PSC_patch
      use coeff_pml

      implicit none
      include './mpif.h'

      integer nodei,nodej,nodek,pec,i1,i2,i3,ierr
      integer mtag,status(MPI_STATUS_SIZE)
      integer rysize
      
      real(kind=8) restj
      real(kind=8) fd(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &                i3mn-rd3:i3mx+rd3)

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


c INITIALIZATION


      mtag=300

      nodei=seg_i1(mpe)
      nodej=seg_i2(mpe)
      nodek=seg_i3(mpe)

      restj=nodej/2.0-int(nodej/2.0)                              ! restj=0.5 => nodej odd

      rysize=(i1mx-i1mn+2*rd1+1)*rd2*(i3mx-i3mn+2*rd3+1)

      allocate(rimy(i1mn-rd1:i1mx+rd1,1:rd2,i3mn-rd3:i3mx+rd3))


c UPDATING LOCAL e4


      if (nodej.lt.ynpe.and.restj<0.25) then
         pec=seg_inv(nodei,nodej+1,nodek)
         do i3=i3mn-rd3,i3mx+rd3
            do i2=1,rd2
               do i1=i1mn-rd1,i1mx+rd1
                  rimy(i1,i2,i3)=fd(i1,i2mx-i2+1,i3)
               enddo
            enddo
         enddo
         call MPI_SSEND(rimy,rysize,MPI_DOUBLE_PRECISION,
     &                  pec,mtag,MPI_COMM_WORLD,ierr)
      endif
      if (1.lt.nodej.and.restj>0.25) then
         pec=seg_inv(nodei,nodej-1,nodek)
         call MPI_RECV(rimy,rysize,MPI_DOUBLE_PRECISION,
     &                 pec,mtag,MPI_COMM_WORLD,status,ierr)
         do i3=i3mn-rd3,i3mx+rd3
            do i2=1,rd2
               do i1=i1mn-rd1,i1mx+rd1
                  fd(i1,i2mn-i2,i3)=rimy(i1,i2,i3)
               enddo
            enddo
         enddo
      endif

      if (nodej.lt.ynpe.and.restj>0.25) then
         pec=seg_inv(nodei,nodej+1,nodek)
         do i3=i3mn-rd3,i3mx+rd3
            do i2=1,rd2
               do i1=i1mn-rd1,i1mx+rd1
                  rimy(i1,i2,i3)=fd(i1,i2mx-i2+1,i3)
               enddo
            enddo
         enddo
         call MPI_SSEND(rimy,rysize,MPI_DOUBLE_PRECISION,
     &                  pec,mtag,MPI_COMM_WORLD,ierr)
      endif
      if (1.lt.nodej.and.restj<0.25) then
         pec=seg_inv(nodei,nodej-1,nodek)
         call MPI_RECV(rimy,rysize,MPI_DOUBLE_PRECISION,
     &                 pec,mtag,MPI_COMM_WORLD,status,ierr)
         do i3=i3mn-rd3,i3mx+rd3
            do i2=1,rd2
               do i1=i1mn-rd1,i1mx+rd1
                  fd(i1,i2mn-i2,i3)=rimy(i1,i2,i3)
               enddo
            enddo
         enddo
      endif


c UPDATING LOCAL BOUNDARY e4 
! periodic continuation at y1 (data from y2)


      if (boundary_fld_y2==1) then       
      if (ynpe.gt.1) then
         if (nodej.eq.ynpe) then
            pec=seg_inv(nodei,1,nodek)
            do i3=i3mn-rd3,i3mx+rd3
               do i2=1,rd2
                  do i1=i1mn-rd1,i1mx+rd1
                     rimy(i1,i2,i3)=fd(i1,i2mx-i2+1,i3)
                  enddo
               enddo
            enddo
            call MPI_SSEND(rimy,rysize,MPI_DOUBLE_PRECISION,
     &                     pec,mtag,MPI_COMM_WORLD,ierr)
         endif
         if (1.eq.nodej) then
            pec=seg_inv(nodei,ynpe,nodek)
            call MPI_RECV(rimy,rysize,MPI_DOUBLE_PRECISION,
     &                    pec,mtag,MPI_COMM_WORLD,status,ierr)
            do i3=i3mn-rd3,i3mx+rd3
               do i2=1,rd2
                  do i1=i1mn-rd1,i1mx+rd1
                     fd(i1,i2mn-i2,i3)=rimy(i1,i2,i3)
                  enddo
               enddo
            enddo
         endif
      else
         do i3=i3mn-rd3,i3mx+rd3
            do i2=1,rd2
               do i1=i1mn-rd1,i1mx+rd1
                  rimy(i1,i2,i3)=fd(i1,i2mx-i2+1,i3)
                  fd(i1,i2mn-i2,i3)=rimy(i1,i2,i3)
               enddo
            enddo
         enddo
      endif
      endif

! pml continuation at y1

      if (nodej.eq.1) then
         if (boundary_fld_y1.eq.2.or.
     &      (boundary_fld_y1.eq.3.and.
     &       pos_y1.ne.0.0.and.pos_y1.lt.n*dt)) then
            do i2 = i2mn-rd2,i2mn+thick
               kappay(i2) = 1.0 + (kappay_max-1.0)
     &              *((thick+1-(i2-i2mn))*dy/deltay)**pml
               sigmay(i2) = sigmay_max*
     &              ((thick+1-(i2-i2mn))*dy/deltay)**pml
               cyp(i2) = 2*eps0*kappay(i2)+sigmay(i2)*dt
               cym(i2) = 2*eps0*kappay(i2)-sigmay(i2)*dt
               fby(i2) = 2*eps0*kappay(i2)
               fcy(i2) = cym(i2)/cyp(i2)
               fdy(i2) = 2*eps0*dt/cyp(i2)
               fey(i2) = 1.0/cyp(i2)
            end do

            do i2 = i2mn-rd2,i2mn+thick
               kappay(i2) = 1.0 + (kappay_max-1.0)
     &              *((thick+1-(i2-i2mn+0.5))*dy/deltay)**pml
               sigmay(i2) = sigmay_max*
     &              ((thick+1-(i2-i2mn+0.5))*dy/deltay)**pml
               byp(i2) = 2*eps0*kappay(i2)+sigmay(i2)*dt
               bym(i2) = 2*eps0*kappay(i2)-sigmay(i2)*dt
               gby(i2) = 2*eps0*kappay(i2)
               gcy(i2) = bym(i2)/byp(i2)
               gdy(i2) = 2*eps0*dt/byp(i2)
               gey(i2) = 1.0/byp(i2)
            end do
c            boundary_pml_y1 = 'done'
         endif
      endif

c UPDATING LOCAL e2


      if (1.lt.nodej.and.restj<0.25) then
         pec=seg_inv(nodei,nodej-1,nodek)
         do i3=i3mn-rd3,i3mx+rd3
            do i2=1,rd2
               do i1=i1mn-rd1,i1mx+rd1
                  rimy(i1,i2,i3)=fd(i1,i2mn+i2-1,i3)
               enddo
            enddo
         enddo
         call MPI_SSEND(rimy,rysize,MPI_DOUBLE_PRECISION,
     &                  pec,mtag,MPI_COMM_WORLD,ierr)
      endif
      if (nodej.lt.ynpe.and.restj>0.25) then
         pec=seg_inv(nodei,nodej+1,nodek)
         call MPI_RECV(rimy,rysize,MPI_DOUBLE_PRECISION,
     &                 pec,mtag,MPI_COMM_WORLD,status,ierr)
         do i3=i3mn-rd3,i3mx+rd3
            do i2=1,rd2
               do i1=i1mn-rd1,i1mx+rd1
                  fd(i1,i2mx+i2,i3)=rimy(i1,i2,i3)
               enddo
            enddo
         enddo
      endif

      if (1.lt.nodej.and.restj>0.25) then
         pec=seg_inv(nodei,nodej-1,nodek)
         do i3=i3mn-rd3,i3mx+rd3
            do i2=1,rd2
               do i1=i1mn-rd1,i1mx+rd1
                  rimy(i1,i2,i3)=fd(i1,i2mn+i2-1,i3)
               enddo
            enddo
         enddo
         call MPI_SSEND(rimy,rysize,MPI_DOUBLE_PRECISION,
     &                  pec,mtag,MPI_COMM_WORLD,ierr)
      endif
      if (nodej.lt.ynpe.and.restj<0.25) then
         pec=seg_inv(nodei,nodej+1,nodek)
         call MPI_RECV(rimy,rysize,MPI_DOUBLE_PRECISION,
     &                 pec,mtag,MPI_COMM_WORLD,status,ierr)
         do i3=i3mn-rd3,i3mx+rd3
            do i2=1,rd2
               do i1=i1mn-rd1,i1mx+rd1
                  fd(i1,i2mx+i2,i3)=rimy(i1,i2,i3)
               enddo
            enddo
         enddo
      endif


c UPDATING LOCAL BOUNDARY e2 
! periodic continuation at y2 (data from y1)


      if (boundary_fld_y1==1) then       
      if (ynpe.gt.1) then
         if (1.eq.nodej) then
            pec=seg_inv(nodei,ynpe,nodek)
            do i3=i3mn-rd3,i3mx+rd3
               do i2=1,rd2
                  do i1=i1mn-rd1,i1mx+rd1 
                     rimy(i1,i2,i3)=fd(i1,i2mn+i2-1,i3)
                  enddo
               enddo
            enddo
            call MPI_SSEND(rimy,rysize,MPI_DOUBLE_PRECISION,
     &                     pec,mtag,MPI_COMM_WORLD,ierr)
         endif
         if (nodej.eq.ynpe) then
            pec=seg_inv(nodei,1,nodek)
            call MPI_RECV(rimy,rysize,MPI_DOUBLE_PRECISION,
     &                    pec,mtag,MPI_COMM_WORLD,status,ierr)
            do i3=i3mn-rd3,i3mx+rd3
               do i2=1,rd2
                  do i1=i1mn-rd1,i1mx+rd1
                     fd(i1,i2mx+i2,i3)=rimy(i1,i2,i3)
                  enddo
               enddo
            enddo
         endif
      else
         do i3=i3mn-rd3,i3mx+rd3
            do i2=1,rd2
               do i1=i1mn-rd1,i1mx+rd1 
                  rimy(i1,i2,i3)=fd(i1,i2mn+i2-1,i3)
                  fd(i1,i2mx+i2,i3)=rimy(i1,i2,i3)
               enddo
            enddo
         enddo
      endif
      endif


! pml continuation at y2

      if (nodej.eq.ynpe) then
         if (boundary_fld_y2.eq.2.or.
     &      (boundary_fld_y2.eq.3.and.
     &       pos_y2.ne.0.0.and.pos_y2.lt.n*dt)) then         
            do i2 = i2mx-thick,i2mx+rd2
               kappay(i2) = 1.0 + (kappay_max-1.0)
     &              *((i2-i2mx+thick)*dy/deltay)**pml
               sigmay(i2) = sigmay_max*
     &              ((i2-i2mx+thick)*dy/deltay)**pml
               cyp(i2) = 2*eps0*kappay(i2)+sigmay(i2)*dt
               cym(i2) = 2*eps0*kappay(i2)-sigmay(i2)*dt
               fby(i2) = 2*eps0*kappay(i2)
               fcy(i2) = cym(i2)/cyp(i2)
               fdy(i2) = 2*eps0*dt/cyp(i2)
               fey(i2) = 1.0/cyp(i2)
            end do

            do i2 = i2mx-thick,i2mx+rd2
               kappay(i2) = 1.0 + (kappay_max-1.0)
     &              *(((i2+0.5)-i2mx+thick)*dy/deltay)**pml
               sigmay(i2) = sigmay_max*
     &              (((i2+0.5)-i2mx+thick)*dy/deltay)**pml
               byp(i2) = 2*eps0*kappay(i2)+sigmay(i2)*dt
               bym(i2) = 2*eps0*kappay(i2)-sigmay(i2)*dt
               gby(i2) = 2*eps0*kappay(i2)
               gcy(i2) = bym(i2)/byp(i2)
               gdy(i2) = 2*eps0*dt/byp(i2)
               gey(i2) = 1.0/byp(i2)
            end do
c            boundary_pml_y2 = 'done'
         endif
      endif

      deallocate(rimy)


      end subroutine PIC_fey
