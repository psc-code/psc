
#include "config.h"

c     ----------------------------------------------------------------------
      subroutine SET_subdomain(i1mn_, i1mx_, i2mn_, i2mx_, i3mn_, i3mx_)
c     ----------------------------------------------------------------------

      use VLA_variables

      implicit none

      integer :: i1mn_, i1mx_, i2mn_, i2mx_, i3mn_, i3mx_

      i1mn = i1mn_
      i2mn = i2mn_
      i3mn = i3mn_
      i1mx = i1mx_
      i2mx = i2mx_
      i3mx = i3mx_

      end subroutine SET_subdomain

c     ----------------------------------------------------------------------
      subroutine GET_subdomain(i1mn_, i1mx_, i2mn_, i2mx_, i3mn_, i3mx_)
c     ----------------------------------------------------------------------

      use VLA_variables

      implicit none

      integer :: i1mn_, i1mx_, i2mn_, i2mx_, i3mn_, i3mx_

      i1mn_ = i1mn
      i2mn_ = i2mn
      i3mn_ = i3mn
      i1mx_ = i1mx
      i2mx_ = i2mx
      i3mx_ = i3mx

      end subroutine GET_subdomain
