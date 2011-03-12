
#include "config.h"

c     wrapper for C to be able to access these parameters

      subroutine GET_param_domain(length, itot, in, ix,
     &     bnd_fld_lo, bnd_fld_hi, bnd_part, nproc, rd, use_pml_)

      use VLA_variables
      use PIC_variables
      use PSC_patch

      implicit none

      real(kind=8) :: length(3)
      integer :: itot(3), in(3), ix(3)
      integer :: bnd_fld_lo(3), bnd_fld_hi(3), bnd_part(3)
      integer :: nproc(3), rd(3)
      logical :: use_pml_

      length(1)=lengthx
      length(2)=lengthy
      length(3)=lengthz

      itot(1)=i1tot
      itot(2)=i2tot
      itot(3)=i3tot

      in(1)=i1n
      in(2)=i2n
      in(3)=i3n

      ix(1)=i1x
      ix(2)=i2x
      ix(3)=i3x

      nproc(1)=xnpe
      nproc(2)=ynpe
      nproc(3)=znpe

      bnd_fld_lo(1)=boundary_fld_x1
      bnd_fld_lo(2)=boundary_fld_y1
      bnd_fld_lo(3)=boundary_fld_z1
      bnd_fld_hi(1)=boundary_fld_x2
      bnd_fld_hi(2)=boundary_fld_y2
      bnd_fld_hi(3)=boundary_fld_z2

      bnd_part(1)=boundary_part_x
      bnd_part(2)=boundary_part_y
      bnd_part(3)=boundary_part_z

      rd(1)=rd1
      rd(2)=rd2
      rd(3)=rd3

      use_pml_ = use_pml

      end subroutine GET_param_domain

c     wrapper for C to be able to set these parameters

      subroutine SET_param_domain(length, itot, in, ix,
     &     bnd_fld_lo, bnd_fld_hi, bnd_part, nproc, rd, use_pml_)

      use VLA_variables
      use PIC_variables
      use PSC_patch

      implicit none

      real(kind=8) :: length(3)
      integer :: itot(3), in(3), ix(3)
      integer :: bnd_fld_lo(3), bnd_fld_hi(3), bnd_part(3)
      integer :: nproc(3), rd(3)
      logical :: use_pml_

      lengthx=length(1)
      lengthy=length(2)
      lengthz=length(3)

      i1tot=itot(1)
      i2tot=itot(2)
      i3tot=itot(3)

      i1n=in(1)
      i2n=in(2)
      i3n=in(3)

      i1x=ix(1)
      i2x=ix(2)
      i3x=ix(3)

      xnpe=nproc(1)
      ynpe=nproc(2)
      znpe=nproc(3)

      boundary_fld_x1=bnd_fld_lo(1)
      boundary_fld_y1=bnd_fld_lo(2)
      boundary_fld_z1=bnd_fld_lo(3)
      boundary_fld_x2=bnd_fld_hi(1)
      boundary_fld_y2=bnd_fld_hi(2)
      boundary_fld_z2=bnd_fld_hi(3)

      boundary_part_x=bnd_part(1)
      boundary_part_y=bnd_part(2)
      boundary_part_z=bnd_part(3)

      rd1=rd(1)
      rd2=rd(2)
      rd3=rd(3)

      use_pml = use_pml_

      end subroutine SET_param_domain

c     wrapper for C to be able to access these parameters

      subroutine GET_param_psc(qq_, mm_, tt_, cc_, eps0_,
     *     nmax_, cpum_, lw_, i0_, n0_, e0_, b0_, j0_, rho0_, phi0_,
     *     a0_, nicell_)

      use VLA_variables
      use PIC_variables

      implicit none

      real(kind=8) :: qq_, mm_, tt_, cc_, eps0_
      real(kind=8) :: cpum_, lw_
      real(kind=8) :: i0_, n0_, e0_, b0_, j0_, rho0_, phi0_, a0_
      integer :: nmax_, nicell_

      qq_ = qq
      mm_ = mm
      tt_ = tt
      cc_ = cc
      eps0_ = eps0

      nmax_ = nmax
      cpum_ = cpum
      lw_ = lw
      i0_ = i0
      n0_ = n0
      e0_ = e0
      b0_ = b0
      j0_ = j0
      rho0_ = rho0
      phi0_ = phi0
      a0_ = a0
      nicell_ = nicell

      end subroutine GET_param_psc

c     wrapper for C to be able to set these parameters

      subroutine SET_param_psc(qq_, mm_, tt_, cc_, eps0_,
     *     nmax_, cpum_, lw_, i0_, n0_, e0_, b0_, j0_, rho0_, phi0_,
     *     a0_, nicell_)

      use VLA_variables
      use PIC_variables

      implicit none

      real(kind=8) :: qq_, mm_, tt_, cc_, eps0_
      real(kind=8) :: cpum_, lw_
      real(kind=8) :: i0_, n0_, e0_, b0_, j0_, rho0_, phi0_, a0_
      integer :: nmax_, nicell_

      qq = qq_
      mm = mm_
      tt = tt_
      cc = cc_
      eps0 = eps0_

      nmax = nmax_
      cpum = cpum_
      lw = lw_
      i0 = i0_
      n0 = n0_
      e0 = e0_
      b0 = b0_
      j0 = j0_
      rho0 = rho0_
      phi0 = phi0_
      a0 = a0_
      nicell = nicell_

      end subroutine SET_param_psc

c     wrapper for C to be able to set these coefficients

      subroutine SET_param_coeff(beta_,
     *     wl_, ld_, vos_, vt_, wp_, dx_, dt_, np_, nnp_)

      use VLA_variables
      use PIC_variables

      implicit none

      real(kind=8) :: alpha_, beta_, eta_
      real(kind=8) :: wl_, ld_, vos_, vt_, wp_, dt_
      real(kind=8) :: dx_(3)
      integer :: np_, nnp_

      beta = beta_
      wl = wl_
      ld = ld_
      vos = vos_
      vt = vt_
      wp = wp_

      dx = dx_(1)
      dy = dx_(2)
      dz = dx_(3)
      dt = dt_

      np = np_
      nnp = nnp_

      end subroutine SET_param_coeff

