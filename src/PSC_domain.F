
      module PSC_domain
      implicit none

      integer :: i1tot,i2tot,i3tot
      integer :: i1n,i1x,i2n,i2x,i3n,i3x

c lengthx: x extension of simulation box in m
c lengthy: y extension of simulation box in m
c lengthz: z extension of simulation box in m

      real(kind=8) lengthx,lengthy,lengthz

c boundary_part_x: particle boundary conditions in x
c boundary_part_y: particle boundary conditions in y
c boundary_part_z: particle boundary conditions in z

      integer boundary_part_x
      integer boundary_part_y
      integer boundary_part_z

! boundary_fld_x1, boundary_fld_x2: boundary condition for bnd in x-front,back
! boundary_fld_y1, boundary_fld_y2: boundary condition for bnd in y-left,right
! boundary_fld_z1, boundary_fld_z2: boundary condition for bnd in z-up,down
!
! condition = 0 : radiating boundaries
! condition = 1 : periodic boundaries
! condition = 2 : upmls
! condition = 3 : radiating boundaries with upmls at a given time

      integer boundary_fld_x1
      integer boundary_fld_x2
      integer boundary_fld_y1
      integer boundary_fld_y2
      integer boundary_fld_z1
      integer boundary_fld_z2

! whether to PML is being used

      logical use_pml

! thick: thickness of pml in gridpoints
! size: thickness of pml and buffer region in gridpoints

      integer thick, size

c npe: total number of PEs in the partition
c mpe: number of the PE the code is running on

      integer npe,mpe

c xnpe: number of PEs in x-direction
c ynpe: number of PEs in y-direction
c znpe: number of PEs in z-direction

      integer :: xnpe,ynpe,znpe

c seg_i1: x-label of distributed data segment
c seg_i2: y-label of distributed data segment
c seg_i3: z-label of distributed data segment
c seg_inv: returns node number after x-, y-label, and z-label input

      integer,allocatable, dimension(:) :: seg_i1,seg_i2,seg_i3
      integer,allocatable, dimension(:,:,:) :: seg_inv

      end module PSC_domain

      subroutine PSC_set_domain(length, itot, in, ix,
     &     bnd_fld_lo, bnd_fld_hi, bnd_part, nproc, rd, use_pml_)

      use PSC_domain
      use PSC_patch

      implicit none

      real(kind=8) :: length(3)
      integer :: itot(3), in(3), ix(3)
      integer :: bnd_fld_lo(3), bnd_fld_hi(3), bnd_part(3)
      integer :: nproc(3), rd(3)
      logical :: use_pml_
      integer :: i1, i2, i3, pec


      lengthx=length(1)
      lengthy=length(2)
      lengthz=length(3)

      i1tot=itot(1)
      i2tot=itot(2)
      i3tot=itot(3)

      i1n=in(1)
      i2n=in(2)
      i3n=in(3)

      i1x=ix(1)
      i2x=ix(2)
      i3x=ix(3)

      xnpe=nproc(1)
      ynpe=nproc(2)
      znpe=nproc(3)

      boundary_fld_x1=bnd_fld_lo(1)
      boundary_fld_y1=bnd_fld_lo(2)
      boundary_fld_z1=bnd_fld_lo(3)
      boundary_fld_x2=bnd_fld_hi(1)
      boundary_fld_y2=bnd_fld_hi(2)
      boundary_fld_z2=bnd_fld_hi(3)

      boundary_part_x=bnd_part(1)
      boundary_part_y=bnd_part(2)
      boundary_part_z=bnd_part(3)

      rd1=rd(1)
      rd2=rd(2)
      rd3=rd(3)

      use_pml = use_pml_

      if (allocated(seg_i1)) then
         return
      endif

      allocate(seg_i1(0:xnpe*ynpe*znpe-1))
      allocate(seg_i2(0:xnpe*ynpe*znpe-1))
      allocate(seg_i3(0:xnpe*ynpe*znpe-1))
      allocate(seg_inv(1:xnpe,1:ynpe,1:znpe))

      do i3=1,znpe
         do i2=1,ynpe
            do i1=1,xnpe
               pec=(i3-1)*xnpe*ynpe+(i2-1)*xnpe+(i1-1)
               seg_i1(pec)=i1
               seg_i2(pec)=i2
               seg_i3(pec)=i3
               seg_inv(i1,i2,i3)=pec
            enddo
         enddo
      enddo

      end subroutine PSC_set_domain

