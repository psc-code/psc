
c     ======================================================================
      module PSC_patch
c     ======================================================================
      implicit none

      integer :: i1mn,i2mn,i3mn ! patch: lower bound
      integer :: i1mx,i2mx,i3mx ! patch: upper bound
      integer :: rd1,rd2,rd3    ! patch: # ghost cells
      real(kind=8) :: dt        ! timestep (dimensionless)
      real(kind=8) :: dx        ! grid spacing (dimensionless)
      real(kind=8) :: dy        ! grid spacing (dimensionless)
      real(kind=8) :: dz        ! grid spacing (dimensionless)
      
      end module PSC_patch

c     ======================================================================
      subroutine PSC_set_patch(imn, imx, rd, dt_, dxyz)
c     ======================================================================
      use PSC_patch
      implicit none

      integer :: imn(3), imx(3), rd(3)
      real(kind=8) :: dt_, dxyz(3)

      i1mn = imn(1)
      i2mn = imn(2)
      i3mn = imn(3)
      i1mx = imx(1)
      i2mx = imx(2)
      i3mx = imx(3)
      rd1 = rd(1)
      rd2 = rd(2)
      rd3 = rd(3)

      dt = dt_
      dx = dxyz(1)
      dy = dxyz(2)
      dz = dxyz(3)
	  
      end subroutine PSC_set_patch

