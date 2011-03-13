
      module PSC_domain
      implicit none

      integer :: i1tot,i2tot,i3tot
      integer :: i1n,i1x,i2n,i2x,i3n,i3x

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

      end module PSC_domain
