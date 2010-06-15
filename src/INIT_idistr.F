
c THIS SUBROUTINE INITIALIZES THE PARTICLE DISTRIBUTION
c AND LOAD BALANCE.

c***************************************************************
c Keine adaptiven Gewichte
c***************************************************************


      subroutine INIT_idistr(part_label_off,
     *     rd1n,rd1x,rd2n,rd2x,rd3n,rd3x)

      use PIC_variables
      use VLA_variables

      implicit none
      include './mpif.h'

      integer :: part_label_off
      integer :: rd1n,rd1x,rd2n,rd2x,rd3n,rd3x

      integer :: cell,ncel
      integer :: rds1,rds2,rds3

      integer :: xmin,xmax     ! ab
      integer :: ymin,ymax     ! ab
      integer :: zmin,zmax     ! ab

      real(kind=8) :: x,y,z
      real(kind=8) :: px,py,pz
      real(kind=8) :: ran1,ran2,ran3,ran4,ran5,ran6
      real(kind=8) :: qni,mni,cni,lni,wni,tnxi,tnyi,tnzi
      real(kind=8) :: dens,INIT_den

      real(kind=8),allocatable,dimension(:) :: rndmv

      allocate(rndmv(1:6*npe))

      xmin=i1n       ! ab
      xmax=i1x       ! ab
      ymin=i2n       ! ab
      ymax=i2x       ! ab
      zmin=i3n       ! ab
      zmax=i3x       ! ab

      if (boundary_pml_x1.ne.'false') xmin=i1n+size+1       ! ab
      if (boundary_pml_x2.ne.'false') xmax=i1x-size-1       ! ab
      if (boundary_pml_y1.ne.'false') ymin=i2n+size+1       ! ab
      if (boundary_pml_y2.ne.'false') ymax=i2x-size-1       ! ab
      if (boundary_pml_z1.ne.'false') zmin=i3n+size+1       ! ab
      if (boundary_pml_z2.ne.'false') zmax=i3x-size-1       ! ab

c SETTING UP PARTICLE PHASE SPACE

      p_niloc=0.0d0
      niloc=0
      do i3=i3mn+rd3n,i3mx+rd3x
         z=i3*dz
         do i2=i2mn+rd2n,i2mx+rd2x
            y=i2*dy
            do i1=i1mn+rd1n,i1mx+rd1x
               x=i1*dx

               cell=(i1-i1mn+rds1+1)
     &              +(i1mx-i1mn+2*rds1+1)
     &              *(i2-i2mn+rds2)
     &              +(i1mx-i1mn+2*rds1+1)
     &              *(i2mx-i2mn+2*rds2+1)
     &              *(i3-i3mn+rds3)

               dens=INIT_den(x,y,z)*pmlcheck(i1,i2,i3)    ! ab
               if (dens.gt.0.001) then
               ncel=nint(dens/cori)
               do l=1,ncel

                  qni=-1.0               ! e-
                  mni=+1.0
                  cni=+cell
                  lni=+part_label_off+niloc
                  wni=+dens
                  tnxi=0.0
                  tnyi=0.0
                  tnzi=0.05

c MAXWELLIAN DISTRIBUTION

                  call random_number(rndmv)
                  ran1=min(0.99999999999999999d0,rndmv(6*mpe+1))
                  ran2=rndmv(6*mpe+2)
                  ran3=min(0.99999999999999999d0,rndmv(6*mpe+3))
                  ran4=rndmv(6*mpe+4)
                  ran5=min(0.99999999999999999d0,rndmv(6*mpe+5))
                  ran6=rndmv(6*mpe+6)

                  px=sqrt(-tnxi*beta**2*log(1.0-ran1)/mni)
     &               *cos(6.2831853*ran2)
                  py=sqrt(-tnyi*beta**2*log(1.0-ran3)/mni)
     &               *cos(6.2831853*ran4)
                  pz=sqrt(-tnzi*beta**2*log(1.0-ran5)/mni)
     &               *cos(6.2831853*ran6)

                  niloc=niloc+1
                  p_niloc(11*niloc+0)=x
                  p_niloc(11*niloc+1)=y
                  p_niloc(11*niloc+2)=z
                  p_niloc(11*niloc+3)=px
                  p_niloc(11*niloc+4)=py
                  p_niloc(11*niloc+5)=pz
                  p_niloc(11*niloc+6)=qni
                  p_niloc(11*niloc+7)=mni
                  p_niloc(11*niloc+8)=cni
                  p_niloc(11*niloc+9)=lni
                  p_niloc(11*niloc+10)=wni

               enddo
               endif

            enddo
         enddo
      enddo


      do i3=i3mn+rd3n,i3mx+rd3x
         z=i3*dz
         do i2=i2mn+rd2n,i2mx+rd2x
            y=i2*dy
            do i1=i1mn+rd1n,i1mx+rd1x
               x=i1*dx

               cell=(i1-i1mn+rds1+1)
     &              +(i1mx-i1mn+2*rds1+1)
     &              *(i2-i2mn+rds2)
     &              +(i1mx-i1mn+2*rds1+1)
     &              *(i2mx-i2mn+2*rds2+1)
     &              *(i3-i3mn+rds3)

               dens=INIT_den(x,y,z)*pmlcheck(i1,i2,i3)    ! ab
               if (dens.gt.0.001) then
               ncel=nint(dens/cori)
               do l=1,ncel

                  qni=+1.0               ! p+
                  mni=+1836.0
                  cni=+cell
                  lni=+part_label_off+niloc
                  wni=+dens
                  tnxi=+0.0
                  tnyi=+0.0
                  tnzi=+0.0

c MAXWELLIAN DISTRIBUTION

                  call random_number(rndmv)
                  ran1=min(0.99999999999999999d0,rndmv(6*mpe+1))
                  ran2=rndmv(6*mpe+2)
                  ran3=min(0.99999999999999999d0,rndmv(6*mpe+3))
                  ran4=rndmv(6*mpe+4)
                  ran5=min(0.99999999999999999d0,rndmv(6*mpe+5))
                  ran6=rndmv(6*mpe+6)

                  px=sqrt(-tnxi*beta**2*log(1.0-ran1)/mni)
     &               *cos(6.2831853*ran2)
                  py=sqrt(-tnyi*beta**2*log(1.0-ran3)/mni)
     &               *cos(6.2831853*ran4)
                  pz=sqrt(-tnzi*beta**2*log(1.0-ran5)/mni)
     &               *cos(6.2831853*ran6)

                  niloc=niloc+1
                  p_niloc(11*niloc+0)=x
                  p_niloc(11*niloc+1)=y
                  p_niloc(11*niloc+2)=z
                  p_niloc(11*niloc+3)=px
                  p_niloc(11*niloc+4)=py
                  p_niloc(11*niloc+5)=pz
                  p_niloc(11*niloc+6)=qni
                  p_niloc(11*niloc+7)=mni
                  p_niloc(11*niloc+8)=cni
                  p_niloc(11*niloc+9)=lni
                  p_niloc(11*niloc+10)=wni

               enddo
               endif

            enddo
         enddo
      enddo


c CHECK FOR ZERO PARTICLE WEIGHTS


      if (niloc.gt.0) then
         do l=1,niloc
            wni=p_niloc(11*l+10)
            if (wni.eq.0.0) then
               if (mpe.eq.0) then
                  write(6,*) 'Zero particle weights not permissible!'
               endif
               call MPI_finalize(info)
               stop
            endif
         enddo
      endif

      deallocate(rndmv)

      contains                         ! ab
      function pmlcheck(x,y,z)
      implicit none
      integer :: x,y,z
      real(kind=8) :: pmlcheck
      
      pmlcheck = 1.0
      if (x.lt.xmin.or.x.gt.xmax
     &     .or.y.lt.ymin.or.y.gt.ymax
     &     .or.z.lt.zmin.or.z.gt.zmax) then
         pmlcheck = 0.0
      end if
      end function

      end subroutine INIT_idistr
