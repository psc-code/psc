
      module PSC_params

c i0: irradiance in Wum**2/m**2
c n0: peak atom density in m**3
c lw: laser wavelength
c wl: laser frequency
c qq: reference charge in As
c mm: reference mass in kg
c tt: atom temperature in J
c cc: light velocity in m/s
c eps0: eps0 in As/Vm
c wp: plasma frequency
c vt: thermal velocity
c ld: plasma wavelength

c nmax: maximum number of time steps 

      real(kind=8) i0,n0,lw,wl
      real(kind=8) qq,mm,tt,cc,eps0
      real(kind=8) ld,wp,vos,vt
      real(kind=8) e0,b0,j0,rho0, phi0, a0

      real(kind=8) beta

      integer :: nmax

      end module PSC_params

      subroutine PSC_set_params(qq_, mm_, tt_, cc_, eps0_,
     *     nmax_, lw_, i0_, n0_, e0_, b0_, j0_, rho0_, phi0_,
     *     a0_)

      use PSC_params

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

      lw = lw_
      i0 = i0_
      n0 = n0_
      e0 = e0_
      b0 = b0_
      j0 = j0_
      rho0 = rho0_
      phi0 = phi0_
      a0 = a0_

      end subroutine PSC_set_params

      subroutine PSC_set_coeff(beta_, wl_, ld_, vos_, vt_, wp_)

      use PSC_params

      implicit none

      real(kind=8) :: beta_
      real(kind=8) :: wl_, ld_, vos_, vt_, wp_

      beta = beta_
      wl = wl_
      ld = ld_
      vos = vos_
      vt = vt_
      wp = wp_

      end subroutine PSC_set_coeff

