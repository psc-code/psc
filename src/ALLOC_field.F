
#include "config.h"

c THIS SUBROUTINE ALLOCATES THE ELECTROMAGNETIC FIELDS.

c     ----------------------------------------------------------------------
      subroutine ALLOC_field
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
      allocate(bx(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))
      allocate(by(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))
      allocate(bz(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
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
      bx=0.0d0
      by=0.0d0
      bz=0.0d0
      ne=0.0d0
      ni=0.0d0
      nn=0.0d0
      jxi=0.0d0
      jyi=0.0d0
      jzi=0.0d0


#ifdef USE_C_INTERFACE
c     need to let C code know about the newly allocated fields
      call C_ALLOC_field_cb(ne, ni, nn, jxi, jyi, jzi,
     &     ex, ey, ez, bx, by, bz)
#endif

      end subroutine ALLOC_field

c     ----------------------------------------------------------------------
      subroutine FREE_field
c     ----------------------------------------------------------------------

      use PIC_variables
      use VLA_variables
      use fields

      implicit none

      deallocate(ex)
      deallocate(ey)
      deallocate(ez)
      deallocate(bx)
      deallocate(by)
      deallocate(bz)
      deallocate(jxi)
      deallocate(jyi)
      deallocate(jzi)
      deallocate(ne)
      deallocate(ni)
      deallocate(nn)

      end subroutine FREE_field

