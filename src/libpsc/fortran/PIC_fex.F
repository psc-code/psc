
c This subroutine exchanges the Maxwell fields between
c computation domains for the parallel Maxwell solver.
c Is called in PIC_msa.f and PIC_msb.f.


      subroutine PIC_fex(fd)
      
      use VLA_variables, only: n
      use PSC_patch
      use PSC_domain
      use coeff_pml

      implicit none
      include './mpif.h'

      integer nodei,nodej,nodek,pec,i1,i2,i3,ierr
      integer mtag,status(MPI_STATUS_SIZE)
      integer rxsize
      
      real(kind=8) resti
      real(kind=8) fd(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &                i3mn-rd3:i3mx+rd3)

      real(kind=8),allocatable,dimension(:,:,:) :: rimx


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


c INITIALIZATION


      mtag=300

      nodei=seg_i1(mpe)
      nodej=seg_i2(mpe)
      nodek=seg_i3(mpe)

      resti=nodei/2.0-int(nodei/2.0)                        ! resti=0.5 => nodei odd
      rxsize=rd1*(i2mx-i2mn+2*rd2+1)*(i3mx-i3mn+2*rd3+1)

      allocate(rimx(1:rd1,i2mn-rd2:i2mx+rd2,i3mn-rd3:i3mx+rd3))


c UPDATING LOCAL e1


      if (nodei.lt.xnpe.and.resti<0.25) then
         pec=seg_inv(nodei+1,nodej,nodek)
         do i3=i3mn-rd3,i3mx+rd3
            do i2=i2mn-rd2,i2mx+rd2
               do i1=1,rd1
                  rimx(i1,i2,i3)=fd(i1mx-i1+1,i2,i3)
               enddo
            enddo
         enddo
         call MPI_SSEND(rimx,rxsize,MPI_DOUBLE_PRECISION,
     &                  pec,mtag,MPI_COMM_WORLD,ierr)
      endif
      if (1.lt.nodei.and.resti>0.25) then
         pec=seg_inv(nodei-1,nodej,nodek)
         call MPI_RECV(rimx,rxsize,MPI_DOUBLE_PRECISION,
     &                 pec,mtag,MPI_COMM_WORLD,status,ierr)
         do i3=i3mn-rd3,i3mx+rd3
            do i2=i2mn-rd2,i2mx+rd2
               do i1=1,rd1
                  fd(i1mn-i1,i2,i3)=rimx(i1,i2,i3)
               enddo
            enddo
         enddo
      endif

      if (nodei.lt.xnpe.and.resti>0.25) then
         pec=seg_inv(nodei+1,nodej,nodek)
         do i3=i3mn-rd3,i3mx+rd3
            do i2=i2mn-rd2,i2mx+rd2
               do i1=1,rd1
                  rimx(i1,i2,i3)=fd(i1mx-i1+1,i2,i3)
               enddo
            enddo
         enddo
         call MPI_SSEND(rimx,rxsize,MPI_DOUBLE_PRECISION,
     &                  pec,mtag,MPI_COMM_WORLD,ierr)
      endif
      if (1.lt.nodei.and.resti<0.25) then
         pec=seg_inv(nodei-1,nodej,nodek)
         call MPI_RECV(rimx,rxsize,MPI_DOUBLE_PRECISION,
     &                 pec,mtag,MPI_COMM_WORLD,status,ierr)
         do i3=i3mn-rd3,i3mx+rd3
            do i2=i2mn-rd2,i2mx+rd2
               do i1=1,rd1
                  fd(i1mn-i1,i2,i3)=rimx(i1,i2,i3)
               enddo
            enddo
         enddo
      endif


c UPDATING LOCAL BOUNDARY e1
! periodic continuation at x1 (data from x2)

      if (boundary_fld_x2==1) then       
      if (xnpe.gt.1) then
         if (nodei.eq.xnpe) then
            pec=seg_inv(1,nodej,nodek)
            do i3=i3mn-rd3,i3mx+rd3
               do i2=i2mn-rd2,i2mx+rd2
                  do i1=1,rd1
                     rimx(i1,i2,i3)=fd(i1mx-i1+1,i2,i3)
                  enddo
               enddo
            enddo
            call MPI_SSEND(rimx,rxsize,MPI_DOUBLE_PRECISION,
     &                     pec,mtag,MPI_COMM_WORLD,ierr)
         endif
         if (1.eq.nodei) then
            pec=seg_inv(xnpe,nodej,nodek)
            call MPI_RECV(rimx,rxsize,MPI_DOUBLE_PRECISION,
     &                    pec,mtag,MPI_COMM_WORLD,status,ierr)
            do i3=i3mn-rd3,i3mx+rd3
               do i2=i2mn-rd2,i2mx+rd2
                  do i1=1,rd1
                     fd(i1mn-i1,i2,i3)=rimx(i1,i2,i3)
                  enddo
               enddo
            enddo
         endif
      else
         do i3=i3mn-rd3,i3mx+rd3
            do i2=i2mn-rd2,i2mx+rd2
               do i1=1,rd1
                  rimx(i1,i2,i3)=fd(i1mx-i1+1,i2,i3)
                  fd(i1mn-i1,i2,i3)=rimx(i1,i2,i3)
               enddo
            enddo
         enddo
      endif
      endif


! pml continuation at x1

      if (nodei.eq.1) then
         if (boundary_fld_x1.eq.2.or.
     &      (boundary_fld_x1.eq.3.and.
     &       pos_x1.ne.0.0.and.pos_x1.lt.n*dt)) then
            do i1 = i1mn-rd1,i1mn+thick
               kappax(i1) = 1.0 + (kappax_max-1.0)
     &              *((thick+1-(i1-i1mn))*dx/deltax)**pml
               sigmax(i1) = sigmax_max*
     &              ((thick+1-(i1-i1mn))*dx/deltax)**pml
               cxp(i1) = 2*eps0*kappax(i1)+sigmax(i1)*dt
               cxm(i1) = 2*eps0*kappax(i1)-sigmax(i1)*dt
               fbx(i1) = 2*eps0*kappax(i1)
               fcx(i1) = cxm(i1)/cxp(i1)
               fdx(i1) = 2*eps0*dt/cxp(i1)
               fex(i1) = 1.0/cxp(i1)
            end do

            do i1 = i1mn-rd1,i1mn+thick
               kappax(i1) = 1.0 + (kappax_max-1.0)
     &              *((thick+1-(i1-i1mn+0.5))*dx/deltax)**pml
               sigmax(i1) = sigmax_max*
     &              ((thick+1-(i1-i1mn+0.5))*dx/deltax)**pml
               bxp(i1) = 2*eps0*kappax(i1)+sigmax(i1)*dt
               bxm(i1) = 2*eps0*kappax(i1)-sigmax(i1)*dt
               gbx(i1) = 2*eps0*kappax(i1)
               gcx(i1) = bxm(i1)/bxp(i1)
               gdx(i1) = 2*eps0*dt/bxp(i1)
               gex(i1) = 1.0/bxp(i1)
            end do 
c            boundary_pml_x1 = 'done'
         endif
      endif

c UPDATING LOCAL e3


      if (1.lt.nodei.and.resti<0.25) then
         pec=seg_inv(nodei-1,nodej,nodek)
         do i3=i3mn-rd3,i3mx+rd3
            do i2=i2mn-rd2,i2mx+rd2
               do i1=1,rd1
                  rimx(i1,i2,i3)=fd(i1mn+i1-1,i2,i3)
               enddo
            enddo
         enddo
         call MPI_SSEND(rimx,rxsize,MPI_DOUBLE_PRECISION,
     &                  pec,mtag,MPI_COMM_WORLD,ierr)
      endif
      if (nodei.lt.xnpe.and.resti>0.25) then
         pec=seg_inv(nodei+1,nodej,nodek)
         call MPI_RECV(rimx,rxsize,MPI_DOUBLE_PRECISION,
     &                 pec,mtag,MPI_COMM_WORLD,status,ierr)
         do i3=i3mn-rd3,i3mx+rd3
            do i2=i2mn-rd2,i2mx+rd2
               do i1=1,rd1
                  fd(i1mx+i1,i2,i3)=rimx(i1,i2,i3)
               enddo
            enddo
         enddo
      endif

      if (1.lt.nodei.and.resti>0.25) then
         pec=seg_inv(nodei-1,nodej,nodek)
         do i3=i3mn-rd3,i3mx+rd3
            do i2=i2mn-rd2,i2mx+rd2
               do i1=1,rd1
                  rimx(i1,i2,i3)=fd(i1mn+i1-1,i2,i3)
               enddo
            enddo
         enddo
         call MPI_SSEND(rimx,rxsize,MPI_DOUBLE_PRECISION,
     &                  pec,mtag,MPI_COMM_WORLD,ierr)
      endif
      if (nodei.lt.xnpe.and.resti<0.25) then
         pec=seg_inv(nodei+1,nodej,nodek)
         call MPI_RECV(rimx,rxsize,MPI_DOUBLE_PRECISION,
     &                 pec,mtag,MPI_COMM_WORLD,status,ierr)
         do i3=i3mn-rd3,i3mx+rd3
            do i2=i2mn-rd2,i2mx+rd2
               do i1=1,rd1
                  fd(i1mx+i1,i2,i3)=rimx(i1,i2,i3)
               enddo
            enddo
         enddo
      endif


c UPDATING LOCAL BOUNDARY e3 
! periodic continuation at x2 (data from x1)

      if (boundary_fld_x1==1) then       
      if (xnpe.gt.1) then
         if (1.eq.nodei) then
            pec=seg_inv(xnpe,nodej,nodek)
            do i3=i3mn-rd3,i3mx+rd3
               do i2=i2mn-rd2,i2mx+rd2
                  do i1=1,rd1 
                     rimx(i1,i2,i3)=fd(i1mn+i1-1,i2,i3)
                  enddo
               enddo
            enddo
            call MPI_SSEND(rimx,rxsize,MPI_DOUBLE_PRECISION,
     &                     pec,mtag,MPI_COMM_WORLD,ierr)
         endif
         if (nodei.eq.xnpe) then
            pec=seg_inv(1,nodej,nodek)
            call MPI_RECV(rimx,rxsize,MPI_DOUBLE_PRECISION,
     &                    pec,mtag,MPI_COMM_WORLD,status,ierr)
            do i3=i3mn-rd3,i3mx+rd3
               do i2=i2mn-rd2,i2mx+rd2
                  do i1=1,rd1
                     fd(i1mx+i1,i2,i3)=rimx(i1,i2,i3)
                  enddo
               enddo
            enddo
         endif
      else
         do i3=i3mn-rd3,i3mx+rd3
            do i2=i2mn-rd2,i2mx+rd2
               do i1=1,rd1 
                  rimx(i1,i2,i3)=fd(i1mn+i1-1,i2,i3)
                  fd(i1mx+i1,i2,i3)=rimx(i1,i2,i3)
               enddo
            enddo
         enddo
      endif
      endif


! pml continuation at x2

      if (nodei.eq.xnpe) then
         if (boundary_fld_x2.eq.2.or.
     &      (boundary_fld_x2.eq.3.and.
     &       pos_x2.ne.0.0.and.pos_x2.lt.n*dt)) then
            do i1 = i1mx-thick,i1mx+rd1
               kappax(i1) = 1.0 + (kappax_max-1.0)
     &              *((i1-i1mx+thick)*dx/deltax)**pml
               sigmax(i1) = sigmax_max*
     &              ((i1-i1mx+thick)*dx/deltax)**pml
               cxp(i1) = 2*eps0*kappax(i1)+sigmax(i1)*dt
               cxm(i1) = 2*eps0*kappax(i1)-sigmax(i1)*dt
               fbx(i1) = 2*eps0*kappax(i1)
               fcx(i1) = cxm(i1)/cxp(i1)
               fdx(i1) = 2*eps0*dt/cxp(i1)
               fex(i1) = 1.0/cxp(i1)
            end do

            do i1 = i1mx-thick,i1mx+rd1
               kappax(i1) = 1.0 + (kappax_max-1.0)
     &              *(((i1+0.5)-i1mx+thick)*dx/deltax)**pml
               sigmax(i1) = sigmax_max*
     &              (((i1+0.5)-i1mx+thick)*dx/deltax)**pml
               bxp(i1) = 2*eps0*kappax(i1)+sigmax(i1)*dt
               bxm(i1) = 2*eps0*kappax(i1)-sigmax(i1)*dt
               gbx(i1) = 2*eps0*kappax(i1)
               gcx(i1) = bxm(i1)/bxp(i1)
               gdx(i1) = 2*eps0*dt/bxp(i1)
               gex(i1) = 1.0/bxp(i1)
            end do
c            boundary_pml_x2 = 'done'
         endif
      endif

      deallocate(rimx)


      end subroutine PIC_fex
