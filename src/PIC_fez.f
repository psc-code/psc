c This subroutine exchanges the Maxwell fields between
c computation domains for the parallel Maxwell solver.
c Is called in PIC_msa.f and PIC_msb.f. 


      subroutine PIC_fez(fd)
      
      use PIC_variables
      use VLA_variables

      implicit none
      include './mpif.h'

      integer nodei,nodej,nodek
      integer mtag,status(MPI_STATUS_SIZE)
      integer rzsize
      
      real(kind=8) restk
      real(kind=8) fd(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &                i3mn-rd3:i3mx+rd3)

      real(kind=8),allocatable,dimension(:,:,:) :: rimz


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

      restk=nodek/2.0-int(nodek/2.0)                        ! restk=0.5 => nodek odd

      rzsize=(i1mx-i1mn+2*rd1+1)*(i2mx-i2mn+2*rd2+1)*rd3

      allocate(rimz(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,1:rd3))


c UPDATING LOCAL e6


      if (nodek.lt.znpe.and.restk<0.25) then
         pec=seg_inv(nodei,nodej,nodek+1)
         do i3=1,rd3
            do i2=i2mn-rd2,i2mx+rd2
               do i1=i1mn-rd1,i1mx+rd1
                  rimz(i1,i2,i3)=fd(i1,i2,i3mx-i3+1)
               enddo
            enddo
         enddo
         call MPI_SSEND(rimz,rzsize,MPI_DOUBLE_PRECISION,
     &                  pec,mtag,MPI_COMM_WORLD,info)
      endif
      if (1.lt.nodek.and.restk>0.25) then
         pec=seg_inv(nodei,nodej,nodek-1)
         call MPI_RECV(rimz,rzsize,MPI_DOUBLE_PRECISION,
     &                 pec,mtag,MPI_COMM_WORLD,status,info)
         do i3=1,rd3
            do i2=i2mn-rd2,i2mx+rd2
               do i1=i1mn-rd1,i1mx+rd1
                  fd(i1,i2,i3mn-i3)=rimz(i1,i2,i3)
               enddo
            enddo
         enddo
      endif

      if (nodek.lt.znpe.and.restk>0.25) then
         pec=seg_inv(nodei,nodej,nodek+1)
         do i3=1,rd3
            do i2=i2mn-rd2,i2mx+rd2
               do i1=i1mn-rd1,i1mx+rd1
                  rimz(i1,i2,i3)=fd(i1,i2,i3mx-i3+1)
               enddo
            enddo
         enddo
         call MPI_SSEND(rimz,rzsize,MPI_DOUBLE_PRECISION,
     &                  pec,mtag,MPI_COMM_WORLD,info)
      endif
      if (1.lt.nodek.and.restk<0.25) then
         pec=seg_inv(nodei,nodej,nodek-1)
         call MPI_RECV(rimz,rzsize,MPI_DOUBLE_PRECISION,
     &                 pec,mtag,MPI_COMM_WORLD,status,info)
         do i3=1,rd3
            do i2=i2mn-rd2,i2mx+rd2
               do i1=i1mn-rd1,i1mx+rd1
                  fd(i1,i2,i3mn-i3)=rimz(i1,i2,i3)
               enddo
            enddo
         enddo
      endif


c UPDATING LOCAL BOUNDARY e6 (periodic continuation in z)


      if (boundary_field_z==1) then       
      if (znpe.gt.1) then
         if (nodek.eq.znpe) then
            pec=seg_inv(nodei,nodej,1)
            do i3=1,rd3
               do i2=i2mn-rd2,i2mx+rd2
                  do i1=i1mn-rd1,i1mx+rd1
                     rimz(i1,i2,i3)=fd(i1,i2,i3mx-i3+1)
                  enddo
               enddo
            enddo
            call MPI_SSEND(rimz,rzsize,MPI_DOUBLE_PRECISION,
     &                     pec,mtag,MPI_COMM_WORLD,info)
         endif
         if (1.eq.nodek) then
            pec=seg_inv(nodei,nodej,znpe)
            call MPI_RECV(rimz,rzsize,MPI_DOUBLE_PRECISION,
     &                    pec,mtag,MPI_COMM_WORLD,status,info)
            do i3=1,rd3
               do i2=i2mn-rd2,i2mx+rd2
                  do i1=i1mn-rd1,i1mx+rd1
                     fd(i1,i2,i3mn-i3)=rimz(i1,i2,i3)
                  enddo
               enddo
            enddo
         endif
      else
         do i3=1,rd3
            do i2=i2mn-rd2,i2mx+rd2
               do i1=i1mn-rd1,i1mx+rd1
                  rimz(i1,i2,i3)=fd(i1,i2,i3mx-i3+1)
                  fd(i1,i2,i3mn-i3)=rimz(i1,i2,i3)
               enddo
            enddo
         enddo
      endif
      endif


c UPDATING LOCAL BOUNDARY e6 (pml continuation in z) - ab

      if (nodek.eq.znpe) then
         if (boundary_pml_z2.eq.'true'.or.
     &        (boundary_pml_z2.eq.'time'.and.
     &        pos_z2.ne.0.0.and.pos_z2.lt.n*dt)) then
 
c Oberer Knoten in z-Richtung
            do i3 = i3mx-thick, i3mx+rd3
               kappaz(i3) = 1.0 + (kappaz_max-1.0)
     &              *((i3-i3mx+thick)*dz/deltaz)**pml
               sigmaz(i3) = sigmaz_max*
     &              ((i3-i3mx+thick)*dz/deltaz)**pml
               czp(i3) = 2*eps0*kappaz(i3)+sigmaz(i3)*dt
               czm(i3) = 2*eps0*kappaz(i3)-sigmaz(i3)*dt
               fbz(i3) = 2*eps0*kappaz(i3)
               fcz(i3) = czm(i3)/czp(i3)
               fdz(i3) = 2*eps0*dt/czp(i3)
               fez(i3) = 1.0/czp(i3)
            end do

            do i3 = i3mx-thick, i3mx+rd3
               kappaz(i3) = 1.0 + (kappaz_max-1.0)
     &              *(((i3+0.5)-i3mx+thick)*dz/deltaz)**pml
               sigmaz(i3) = sigmaz_max*
     &              (((i3+0.5)-i3mx+thick)*dz/deltaz)**pml
               bzp(i3) = 2*eps0*kappaz(i3)+sigmaz(i3)*dt
               bzm(i3) = 2*eps0*kappaz(i3)-sigmaz(i3)*dt
               gbz(i3) = 2*eps0*kappaz(i3)
               gcz(i3) = bzm(i3)/bzp(i3)
               gdz(i3) = 2*eps0*dt/bzp(i3)
               gez(i3) = 1.0/bzp(i3)
            end do

            boundary_pml_z2 = 'done'
c            write(6,*) 'OBEN'
         endif
      endif


c UPDATING LOCAL e5


      if (1.lt.nodek.and.restk<0.25) then
         pec=seg_inv(nodei,nodej,nodek-1)
         do i3=1,rd3
            do i2=i2mn-rd2,i2mx+rd2
               do i1=i1mn-rd1,i1mx+rd1
                  rimz(i1,i2,i3)=fd(i1,i2,i3mn+i3-1)
               enddo
            enddo
         enddo
         call MPI_SSEND(rimz,rzsize,MPI_DOUBLE_PRECISION,
     &                  pec,mtag,MPI_COMM_WORLD,info)
      endif
      if (nodek.lt.znpe.and.restk>0.25) then
         pec=seg_inv(nodei,nodej,nodek+1)
         call MPI_RECV(rimz,rzsize,MPI_DOUBLE_PRECISION,
     &                 pec,mtag,MPI_COMM_WORLD,status,info)
         do i3=1,rd3
            do i2=i2mn-rd2,i2mx+rd2
               do i1=i1mn-rd1,i1mx+rd1
                  fd(i1,i2,i3mx+i3)=rimz(i1,i2,i3)
               enddo
            enddo
         enddo
      endif

      if (1.lt.nodek.and.restk>0.25) then
         pec=seg_inv(nodei,nodej,nodek-1)
         do i3=1,rd3
            do i2=i2mn-rd2,i2mx+rd2
               do i1=i1mn-rd1,i1mx+rd1
                  rimz(i1,i2,i3)=fd(i1,i2,i3mn+i3-1)
               enddo
            enddo
         enddo
         call MPI_SSEND(rimz,rzsize,MPI_DOUBLE_PRECISION,
     &                  pec,mtag,MPI_COMM_WORLD,info)
      endif
      if (nodek.lt.znpe.and.restk<0.25) then
         pec=seg_inv(nodei,nodej,nodek+1)
         call MPI_RECV(rimz,rzsize,MPI_DOUBLE_PRECISION,
     &                 pec,mtag,MPI_COMM_WORLD,status,info)
         do i3=1,rd3
            do i2=i2mn-rd2,i2mx+rd2
               do i1=i1mn-rd1,i1mx+rd1
                  fd(i1,i2,i3mx+i3)=rimz(i1,i2,i3)
               enddo
            enddo
         enddo
      endif


c UPDATING LOCAL BOUNDARY e5 (periodic continuation in z)


      if (boundary_field_z==1) then       
      if (znpe.gt.1) then
         if (1.eq.nodek) then
            pec=seg_inv(nodei,nodej,znpe)
            do i3=1,rd3
               do i2=i2mn-rd2,i2mx+rd2
                  do i1=i1mn-rd1,i1mx+rd1 
                     rimz(i1,i2,i3)=fd(i1,i2,i3mn+i3-1)
                  enddo
               enddo
            enddo
            call MPI_SSEND(rimz,rzsize,MPI_DOUBLE_PRECISION,
     &                     pec,mtag,MPI_COMM_WORLD,info)
         endif
         if (nodek.eq.znpe) then
            pec=seg_inv(nodei,nodej,1)
            call MPI_RECV(rimz,rzsize,MPI_DOUBLE_PRECISION,
     &                    pec,mtag,MPI_COMM_WORLD,status,info)
            do i3=1,rd3
               do i2=i2mn-rd2,i2mx+rd2
                  do i1=i1mn-rd1,i1mx+rd1
                     fd(i1,i2,i3mx+i3)=rimz(i1,i2,i3)
                  enddo
               enddo
            enddo
         endif
      else
         do i3=1,rd3
            do i2=i2mn-rd2,i2mx+rd2
               do i1=i1mn-rd1,i1mx+rd1 
                  rimz(i1,i2,i3)=fd(i1,i2,i3mn+i3-1)
                  fd(i1,i2,i3mx+i3)=rimz(i1,i2,i3)
               enddo
            enddo
         enddo
      endif
      endif


c UPDATING LOCAL BOUNDARY e5 (pml continuation in z)


      if (nodek.eq.1) then
         if (boundary_pml_z1.eq.'true'.or.
     &        (boundary_pml_z1.eq.'time'.and.
     &        pos_z1.ne.0.0.and.pos_z1.lt.n*dt)) then

            open (13, file="pml.txt",status="replace")
            write(13,*) n
            write(13,*) pos_z1
            close(13)

c     Unterer Knoten in z-Richtung
            do i3 = i3mn-rd3, thick
               kappaz(i3) = 1.0 + (kappaz_max-1.0)
     &              *((thick+1-i3)*dz/deltaz)**pml
               sigmaz(i3) = sigmaz_max*
     &              ((thick+1-i3)*dz/deltaz)**pml
               czp(i3) = 2*eps0*kappaz(i3)+sigmaz(i3)*dt
               czm(i3) = 2*eps0*kappaz(i3)-sigmaz(i3)*dt
               fbz(i3) = 2*eps0*kappaz(i3)
               fcz(i3) = czm(i3)/czp(i3)
               fdz(i3) = 2*eps0*dt/czp(i3)
               fez(i3) = 1.0/czp(i3)
            end do
            
            do i3 = i3mn-rd3, thick
               kappaz(i3) = 1.0 + (kappaz_max-1.0)
     &              *((thick+1-(i3+0.5))*dz/deltaz)**pml
               sigmaz(i3) = sigmaz_max*
     &              ((thick+1-(i3+0.5))*dz/deltaz)**pml
               bzp(i3) = 2*eps0*kappaz(i3)+sigmaz(i3)*dt
               bzm(i3) = 2*eps0*kappaz(i3)-sigmaz(i3)*dt
               gbz(i3) = 2*eps0*kappaz(i3)
               gcz(i3) = bzm(i3)/bzp(i3)
               gdz(i3) = 2*eps0*dt/bzp(i3)
               gez(i3) = 1.0/bzp(i3)
            end do
            boundary_pml_z1 = 'done'

c     write(6,*) 'UNTEN'
         endif
      endif


      deallocate(rimz)


      end subroutine PIC_fez
