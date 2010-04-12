c THIS SUBROUTINE IS THE PARTICLE MOVER.


      subroutine PIC_move_part

      use PIC_variables
      use VLA_variables

      if (i1x-i1n.gt.0.and.
     &    i2x-i2n.gt.0.and.
     &    i3x-i3n.gt.0) then 
         call PIC_move_part_xyz
      elseif (i1x-i1n.gt.0.and.
     &    i2x-i2n.gt.0) then
         call PIC_move_part_xy
      elseif (i1x-i1n.gt.0.and.
     &    i3x-i3n.gt.0) then
         call PIC_move_part_xz
      elseif (i2x-i2n.gt.0.and.
     &    i3x-i3n.gt.0) then
         call PIC_move_part_yz
      elseif (i1x-i1n.gt.0) then
         call PIC_move_part_x
      elseif (i2x-i2n.gt.0) then
         call PIC_move_part_y
      elseif (i3x-i3n.gt.0) then
         call PIC_move_part_z
      endif

      end subroutine PIC_move_part
