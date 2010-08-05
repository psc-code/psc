
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

      end module FIELDS

c     ----------------------------------------------------------------------
      subroutine FIELDS_alloc
c     ----------------------------------------------------------------------

      use PIC_variables
      use VLA_variables
      use fields

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


#ifdef USE_C_INTERFACE
c     need to let C code know about the newly allocated fields
      call C_fields_alloc_cb(ne, ni, nn, jxi, jyi, jzi,
     &     ex, ey, ez, hx, hy, hz)
#endif

      end subroutine FIELDS_alloc

c     ----------------------------------------------------------------------
      subroutine FIELDS_free
c     ----------------------------------------------------------------------

      use PIC_variables
      use VLA_variables
      use fields

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

      end subroutine FIELDS_free

