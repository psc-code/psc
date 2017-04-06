
#ifndef PDE_MHD_PUSHEJ_C
#define PDE_MHD_PUSHEJ_C

#include "pde/pde_mhd_primbb.c"
#include "pde/pde_mhd_calc_current.c"

// ----------------------------------------------------------------------
// patch_push_ej_c

static void
patch_push_ej_c(fld3d_t p_Unext, mrc_fld_data_t dt, fld3d_t p_Ucurr,
		fld3d_t p_W, fld3d_t p_zmask)
{
  static fld3d_t p_Jec, p_Bcc;
  fld3d_setup_tmp_compat(&p_Jec, 3, _BX);
  fld3d_setup_tmp_compat(&p_Bcc, 3, _TMP1);
  fld3d_t p_Bcurr = fld3d_make_view(p_Ucurr, BX);

  patch_calc_current_ec(p_Jec, p_Bcurr);
  patch_primbb(p_Bcc, p_Ucurr);
	
  fld3d_foreach(i,j,k, 0, 0) {
    mrc_fld_data_t s = dt * F3S(p_zmask, 0, i,j,k);
    mrc_fld_data_t cx = EC_TO_CC(p_Jec, 0, i,j,k);
    mrc_fld_data_t cy = EC_TO_CC(p_Jec, 1, i,j,k);
    mrc_fld_data_t cz = EC_TO_CC(p_Jec, 2, i,j,k);
    mrc_fld_data_t ffx = s * (cy * F3S(p_Bcc, 2, i,j,k) - cz * F3S(p_Bcc, 1, i,j,k));
    mrc_fld_data_t ffy = s * (cz * F3S(p_Bcc, 0, i,j,k) - cx * F3S(p_Bcc, 2, i,j,k));
    mrc_fld_data_t ffz = s * (cx * F3S(p_Bcc, 1, i,j,k) - cy * F3S(p_Bcc, 0, i,j,k));
    mrc_fld_data_t duu = (ffx * F3S(p_W, VX, i,j,k) +
			  ffy * F3S(p_W, VY, i,j,k) +
			  ffz * F3S(p_W, VZ, i,j,k));

    F3S(p_Unext, RVX, i,j,k) += ffx;
    F3S(p_Unext, RVY, i,j,k) += ffy;
    F3S(p_Unext, RVZ, i,j,k) += ffz;
    F3S(p_Unext, UU , i,j,k) += duu;
  } fld3d_foreach_end;
}

// ----------------------------------------------------------------------
// patch_push_ej_fortran

#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)

#define push_ej_F77 F77_FUNC(push_ej,PUSH_EJ)

void push_ej_F77(real *b1x, real *b1y, real *b1z,
		 real *rv1x, real *rv1y, real *rv1z, real *uu1,
		 real *zmask, real *vx, real *vy, real *vz,
		 real *dt);

static void
patch_push_ej_fortran(fld3d_t p_Unext, mrc_fld_data_t dt, fld3d_t p_Ucurr,
		      fld3d_t p_W, fld3d_t p_zmask)
{
  push_ej_F77(F(p_Ucurr, BX), F(p_Ucurr, BY), F(p_Ucurr, BZ),
	      F(p_Unext, RVX), F(p_Unext, RVY), F(p_Unext, RVZ), F(p_Unext, UU), 
	      F(p_zmask, 0), F(p_W, VX), F(p_W, VY), F(p_W, VZ), &dt);
}

#endif

// ----------------------------------------------------------------------
// patch_push_ej

static void _mrc_unused
patch_push_ej(fld3d_t p_Unext, mrc_fld_data_t dt, fld3d_t p_Ucurr,
	      fld3d_t p_W, fld3d_t p_zmask)
{
  if (s_opt_mhd_push_ej == OPT_MHD_C) {
    patch_push_ej_c(p_Unext, dt, p_Ucurr, p_W, p_zmask);
#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)
  } else if (s_opt_mhd_push_ej == OPT_MHD_FORTRAN) {
    patch_push_ej_fortran(p_Unext, dt, p_Ucurr, p_W, p_zmask);
#endif
  } else {
    assert(0);
  }
}

#endif
