
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

! thick: thickness of pml in gridpoints
! size: thickness of pml and buffer region in gridpoints

      integer thick, size

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
