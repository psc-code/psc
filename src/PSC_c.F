
#include "config.h"

c     wrapper for C to be able to set these parameters

      subroutine SET_param_psc(qq_, mm_, tt_, cc_, eps0_,
     *     nmax_, lw_, i0_, n0_, e0_, b0_, j0_, rho0_, phi0_,
     *     a0_)

      use PSC_params
      use VLA_variables

      implicit none

      real(kind=8) :: qq_, mm_, tt_, cc_, eps0_
      real(kind=8) :: lw_
      real(kind=8) :: i0_, n0_, e0_, b0_, j0_, rho0_, phi0_, a0_
      integer :: nmax_

      qq = qq_
      mm = mm_
      tt = tt_
      cc = cc_
      eps0 = eps0_

      nmax = nmax_
      lw = lw_
      i0 = i0_
      n0 = n0_
      e0 = e0_
      b0 = b0_
      j0 = j0_
      rho0 = rho0_
      phi0 = phi0_
      a0 = a0_

      end subroutine SET_param_psc

c     wrapper for C to be able to set these coefficients

      subroutine SET_param_coeff(beta_,
     *     wl_, ld_, vos_, vt_, wp_, np_, nnp_)

      use PSC_params
      use VLA_variables

      implicit none

      real(kind=8) :: alpha_, beta_, eta_
      real(kind=8) :: wl_, ld_, vos_, vt_, wp_
      integer :: np_, nnp_

      beta = beta_
      wl = wl_
      ld = ld_
      vos = vos_
      vt = vt_
      wp = wp_

      np = np_
      nnp = nnp_

      end subroutine SET_param_coeff

