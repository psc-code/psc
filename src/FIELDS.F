
#include "config.h"

c     ======================================================================
      module FIELDS
c     ======================================================================
      implicit none

c ex-hz: em field arrays 
c ne-nn: charge densities
c jxi-jzi: current densities

      real(kind=8),allocatable,dimension(:,:,:) :: ex,ey,ez
      real(kind=8),allocatable,dimension(:,:,:) :: hx,hy,hz
      real(kind=8),allocatable,dimension(:,:,:) :: jxi,jyi,jzi
      real(kind=8),allocatable,dimension(:,:,:) :: ne,ni,nn

! for the pml the declaration of d and h field is required
! dvx-hz: additional em field arrays
! dvx-dvz is chosen because of the existence of dx-dz for spatial range

      real(kind=8),allocatable,dimension(:,:,:) :: dvx,dvy,dvz   ! added by ab
      real(kind=8),allocatable,dimension(:,:,:) :: bx,by,bz   ! added by ab

! electric and magnetic permittivity tensors

      real(kind=8),allocatable,dimension(:,:,:) :: eps, mu  ! added by ab

      end module FIELDS

c     ----------------------------------------------------------------------
      subroutine FIELDS_alloc
c     ----------------------------------------------------------------------

      use VLA_variables
      use FIELDS
      use PSC_patch

      implicit none

      allocate(ex(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))
      allocate(ey(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))
      allocate(ez(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))
      allocate(hx(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))
      allocate(hy(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))
      allocate(hz(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))
      allocate(jxi(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &             i3mn-rd3:i3mx+rd3))
      allocate(jyi(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &             i3mn-rd3:i3mx+rd3))
      allocate(jzi(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &             i3mn-rd3:i3mx+rd3))
      allocate(ne(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))
      allocate(ni(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))
      allocate(nn(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))

      ex=0.0d0
      ey=0.0d0
      ez=0.0d0
      hx=0.0d0
      hy=0.0d0
      hz=0.0d0
      ne=0.0d0
      ni=0.0d0
      nn=0.0d0
      jxi=0.0d0
      jyi=0.0d0
      jzi=0.0d0

      if (use_pml) call FIELDS_alloc_pml

c     need to let C code know about the newly allocated fields
      call C_fields_alloc_cb(ne, ni, nn, jxi, jyi, jzi,
     &     ex, ey, ez, hx, hy, hz, dvx, dvy, dvz, bx, by, bz, eps, mu)

      end subroutine FIELDS_alloc

c     ----------------------------------------------------------------------
      subroutine FIELDS_free
c     ----------------------------------------------------------------------

      use VLA_variables
      use FIELDS

      implicit none

      deallocate(ex)
      deallocate(ey)
      deallocate(ez)
      deallocate(hx)
      deallocate(hy)
      deallocate(hz)
      deallocate(jxi)
      deallocate(jyi)
      deallocate(jzi)
      deallocate(ne)
      deallocate(ni)
      deallocate(nn)

      if (use_pml) call FIELDS_free_pml

      end subroutine FIELDS_free

c     ======================================================================
c     PML specific extensions

c     ----------------------------------------------------------------------
      subroutine FIELDS_alloc_pml
c     ----------------------------------------------------------------------

      use VLA_variables
      use FIELDS
      use PSC_patch

      implicit none

      allocate(dvx(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))
      allocate(dvy(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))
      allocate(dvz(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))
      allocate(bx(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))
      allocate(by(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))
      allocate(bz(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))
      allocate(eps(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))
      allocate(mu(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))

      dvx=0.0d0
      dvy=0.0d0
      dvz=0.0d0
      bx=0.0d0
      by=0.0d0
      bz=0.0d0
      eps = 1.0
      mu = 1.0

      end subroutine FIELDS_alloc_pml

c     ----------------------------------------------------------------------
      subroutine FIELDS_free_pml
c     ----------------------------------------------------------------------

      use FIELDS

      implicit none

      deallocate(dvx)
      deallocate(dvy)
      deallocate(dvz)
      deallocate(bx)
      deallocate(by)
      deallocate(bz)
      deallocate(eps)
      deallocate(mu)

      end subroutine FIELDS_free_pml

