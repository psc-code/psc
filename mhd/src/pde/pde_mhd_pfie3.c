
#ifndef PDE_MHD_PFIE3_C
#define PDE_MHD_PFIE3_C

#include "pde/pde_mhd_calce.c"
#include "pde/pde_mhd_bpush.c"

// ----------------------------------------------------------------------
// patch_pfie3_c

static void
patch_pfie3_c(fld3d_t p_Unext, mrc_fld_data_t dt, fld3d_t p_Uprev,
	      fld3d_t p_Ucurr, fld3d_t p_W, fld3d_t p_zmask, fld3d_t p_rmask,
	      fld3d_t p_resis, fld3d_t p_Jcc, fld3d_t p_E)
{
  patch_calce(p_E, dt, p_Ucurr, p_W, p_zmask, p_rmask, p_resis, p_Jcc);
  patch_bpush1(p_Unext, dt, p_Uprev, p_E);
}

// ----------------------------------------------------------------------
// patch_pfie3_fortran

#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)

#define pfie3_F77 F77_FUNC(pfie3,PFIE3)

void pfie3_F77(real *b1x, real *b1y, real *b1z,
	       real *b2x, real *b2y, real *b2z,
	       real *b3x, real *b3y, real *b3z,
	       real *zmask, real *rmask, real *resis,
	       real *flx, real *fly, real *flz,
	       real *vx, real *vy, real *vz,
	       real *currx, real *curry, real *currz,
	       real *dt, real *time);

static void
patch_pfie3_fortran(fld3d_t p_Unext, mrc_fld_data_t dt, fld3d_t p_Uprev,
		    fld3d_t p_Ucurr, fld3d_t p_W, fld3d_t p_zmask, fld3d_t p_rmask,
		    fld3d_t p_resis, fld3d_t p_Jcc, fld3d_t p_E)
{
  pfie3_F77(F(p_Uprev, BX), F(p_Uprev, BY), F(p_Uprev, BZ),
  	    F(p_Ucurr, BX), F(p_Ucurr, BY), F(p_Ucurr, BZ),
  	    F(p_Unext, BX), F(p_Unext, BY), F(p_Unext, BZ),
  	    F(p_zmask, 0), F(p_rmask, 0), F(p_resis, 0),
  	    F(p_E, 0), F(p_E, 1), F(p_E, 2),
  	    F(p_W, VX), F(p_W, VY), F(p_W, VZ),
  	    F(p_Jcc, 0), F(p_Jcc, 1), F(p_Jcc, 2),
  	    &dt, &s_mhd_time);
}

#endif

// ----------------------------------------------------------------------
// patch_pfie3
//
// this dispatch can go away if we get rid of the do_legacy part
// in the Fortran version

static void
patch_pfie3(fld3d_t p_Unext, mrc_fld_data_t dt, fld3d_t p_Uprev,
	    fld3d_t p_Ucurr, fld3d_t p_W, fld3d_t p_zmask, fld3d_t p_rmask,
	    fld3d_t p_resis, fld3d_t p_Jcc, fld3d_t p_E)
{
  if (s_opt_mhd_pfie3 == OPT_MHD_C) {
    patch_pfie3_c(p_Unext, dt, p_Uprev, p_Ucurr, p_W, p_zmask, p_rmask, p_resis,
		  p_Jcc, p_E);
#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)
  } else if (s_opt_mhd_pfie3 == OPT_MHD_FORTRAN) {
    patch_pfie3_fortran(p_Unext, dt, p_Uprev, p_Ucurr, p_W, p_zmask, p_rmask, p_resis,
			p_Jcc, p_E);
#endif
  } else {
    assert(0);
  }
}

#endif
