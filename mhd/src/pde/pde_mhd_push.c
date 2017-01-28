
#include "pde/pde_mhd_rmaskn.c"
#include "pde/pde_mhd_pushfluid.c"
#include "pde/pde_mhd_pushfield.c"
#include "pde/pde_mhd_badval_checks.c"

// ----------------------------------------------------------------------
// patch_push_c

static void
patch_push_c(fld3d_t p_Unext, fld3d_t p_Uprev, fld3d_t p_Ucurr,
	     fld3d_t p_W, fld3d_t p_cmsv,
	     fld3d_t p_ymask, fld3d_t p_zmask, fld3d_t p_E,
	     mrc_fld_data_t dt, int stage)
{
  static fld3d_t p_rmask, p_resis, p_Jcc;
  fld3d_setup_tmp_compat(&p_rmask, 1, _RMASK);
  fld3d_setup_tmp_compat(&p_resis, 1, _RESIS);
  fld3d_setup_tmp_compat(&p_Jcc, 3, _CURRX);

  patch_rmaskn(p_rmask, p_zmask);
  patch_pushfluid(p_Unext, dt, p_Uprev, p_Ucurr, p_W,
		  p_cmsv, p_ymask, p_zmask, stage);
  patch_pushfield(p_Unext, dt, p_Uprev, p_Ucurr, p_W,
		  p_zmask, p_rmask, p_resis, p_Jcc, p_E, stage);
}

// ----------------------------------------------------------------------
// patch_push_fortran

#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)

#include "pde/pde_fortran.h"

#define pushpred_F77 F77_FUNC(pushpred,PUSHPRED)
#define pushcorr_F77 F77_FUNC(pushcorr,PUSHCORR)

void pushpred_F77(real *rr1, real *rv1x, real *rv1y, real *rv1z, real *uu1,
		  real *b1x, real *b1y, real *b1z,
		  real *rr2, real *rv2x, real *rv2y, real *rv2z, real *uu2,
		  real *b2x, real *b2y, real *b2z,
		  real *rr, real *vx, real *vy, real *vz, real *pp,
		  real *cmsv, real *ymask, real *zmask, real *rmask,
		  real *flx, real *fly, real *flz,
		  real *tmp1, real *tmp2, real *tmp3, real *resis,
		  real *dt, real *time);
void pushcorr_F77(real *rr1, real *rv1x, real *rv1y, real *rv1z, real *uu1,
		  real *b1x, real *b1y, real *b1z,
		  real *rr2, real *rv2x, real *rv2y, real *rv2z, real *uu2,
		  real *b2x, real *b2y, real *b2z,
		  real *rr, real *vx, real *vy, real *vz, real *pp,
		  real *cmsv, real *ymask, real *zmask, real *rmask,
		  real *flx, real *fly, real *flz,
		  real *tmp1, real *tmp2, real *tmp3, real *resis,
		  real *dt, real *time);

static void
patch_pushpred_fortran(mrc_fld_data_t dt)
{
  pushpred_F77(F(s_p_f, _RR1), F(s_p_f, _RV1X), F(s_p_f, _RV1Y), F(s_p_f, _RV1Z), F(s_p_f, _UU1),
	       F(s_p_f, _B1X), F(s_p_f, _B1Y), F(s_p_f, _B1Z),
	       F(s_p_f, _RR2), F(s_p_f, _RV2X), F(s_p_f, _RV2Y), F(s_p_f, _RV2Z), F(s_p_f, _UU2),
	       F(s_p_f, _B2X), F(s_p_f, _B2Y), F(s_p_f, _B2Z),
	       F(s_p_f, _RR), F(s_p_f, _VX), F(s_p_f, _VY), F(s_p_f, _VZ), F(s_p_f, _PP),
	       F(s_p_f, _CMSV), F(s_p_f, _YMASK), F(s_p_f, _ZMASK), F(s_p_f, _RMASK),
	       F(s_p_f, _FLX), F(s_p_f, _FLY), F(s_p_f, _FLZ),
	       F(s_p_f, _TMP1), F(s_p_f, _TMP2), F(s_p_f, _TMP3), F(s_p_f, _RESIS),
	       &dt, &s_mhd_time);
}

static void
patch_pushcorr_fortran(mrc_fld_data_t dt)
{
  pushcorr_F77(F(s_p_f, _RR1), F(s_p_f, _RV1X), F(s_p_f, _RV1Y), F(s_p_f, _RV1Z), F(s_p_f, _UU1),
	       F(s_p_f, _B1X), F(s_p_f, _B1Y), F(s_p_f, _B1Z),
	       F(s_p_f, _RR2), F(s_p_f, _RV2X), F(s_p_f, _RV2Y), F(s_p_f, _RV2Z), F(s_p_f, _UU2),
	       F(s_p_f, _B2X), F(s_p_f, _B2Y), F(s_p_f, _B2Z),
	       F(s_p_f, _RR), F(s_p_f, _VX), F(s_p_f, _VY), F(s_p_f, _VZ), F(s_p_f, _PP),
	       F(s_p_f, _CMSV), F(s_p_f, _YMASK), F(s_p_f, _ZMASK), F(s_p_f, _RMASK),
	       F(s_p_f, _FLX), F(s_p_f, _FLY), F(s_p_f, _FLZ),
	       F(s_p_f, _TMP1), F(s_p_f, _TMP2), F(s_p_f, _TMP3), F(s_p_f, _RESIS),
	       &dt, &s_mhd_time);
}

static void
patch_push_fortran(mrc_fld_data_t dt, int stage)
{
  if (stage == 0) {
    // fortran expects full timestep, but predictor got pass .5f * mhd->dt
    patch_pushpred_fortran(2.f * dt);
  } else {
    patch_pushcorr_fortran(dt);
  }
}

#endif

// ----------------------------------------------------------------------
// patch_push

static void
patch_push(fld3d_t p_Unext, fld3d_t p_Uprev, fld3d_t p_Ucurr,
	   fld3d_t p_W, fld3d_t p_cmsv,
	   fld3d_t p_ymask, fld3d_t p_zmask,
	   fld3d_t p_E, mrc_fld_data_t dt, int stage)
{
  int opt_mhd_push = stage ? s_opt_mhd_pushpred : s_opt_mhd_pushcorr;

  if (opt_mhd_push == OPT_MHD_C) {
    patch_push_c(p_Unext, p_Uprev, p_Ucurr, p_W, p_cmsv,
		 p_ymask, p_zmask, p_E, dt, stage);
#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)
  } else if (opt_mhd_push == OPT_MHD_FORTRAN) {
    patch_push_fortran(dt, stage);
#endif
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// patch_poison_bnd

static void _mrc_unused
patch_poison_bnd(fld3d_t p_U)
{
#if OPT_STAGGER == OPT_STAGGER_GGCM
  const int S = 1;
#else
  const int S = 0;
#endif

  fld3d_foreach(i,j,k, 2, 2) {
    if (i >= 0 && i < s_ldims[0] &&
	j >= 0 && j < s_ldims[1] &&
	k >= 0 && k < s_ldims[2]) {
      continue;
    }
    for (int m = 0; m < 5; m++) {
      F3S(p_U, m, i,j,k) = 999999.;
    }
    if (i+S < 0 || i+S > s_ldims[0]) {
      F3S(p_U, BX, i,j,k) = 999999.;
    }
    if (j+S < 0 || j+S > s_ldims[1]) {
      F3S(p_U, BY, i,j,k) = 999999.;
    }
    if (k+S < 0 || k+S > s_ldims[2]) {
      F3S(p_U, BZ, i,j,k) = 999999.;
    }
  } fld3d_foreach_end;
}

// ----------------------------------------------------------------------
// patch_pushstage

static void
patch_pushstage(fld3d_t p_f, mrc_fld_data_t dt, int stage)
{
  static fld3d_t p_W, p_cmsv, p_E;
  fld3d_setup_tmp_compat(&p_W   , 5, _RR);
  fld3d_setup_tmp_compat(&p_cmsv, 1, _CMSV);
  fld3d_setup_tmp_compat(&p_E   , 3, _FLX);

  fld3d_t p_Unext, p_Uprev, p_Ucurr;
  if (stage == 0) {
    fld3d_setup_view(&p_Unext, p_f, _RR2);
    fld3d_setup_view(&p_Uprev, p_f, _RR1);
    fld3d_setup_view(&p_Ucurr, p_f, _RR1);
  } else {
    fld3d_setup_view(&p_Unext, p_f, _RR1);
    fld3d_setup_view(&p_Uprev, p_f, _RR1);
    fld3d_setup_view(&p_Ucurr, p_f, _RR2);
  }
  fld3d_t p_ymask = fld3d_make_view(p_f, _YMASK);
  fld3d_t p_zmask = fld3d_make_view(p_f, _ZMASK);

  patch_primvar(p_W, p_Ucurr, p_cmsv);
  if (stage == 1) {
    patch_badval_checks_sc(p_Ucurr, p_W);
  }

  if (stage == 0) {
    static fld3d_t p_bcc;
    fld3d_setup_tmp_compat(&p_bcc, 3, _BX);
    patch_primbb(p_bcc, p_Ucurr);
    patch_zmaskn(p_zmask, p_W, p_bcc, p_ymask);
  }

  patch_push(p_Unext, p_Uprev, p_Ucurr, p_W, p_cmsv,
	     p_ymask, p_zmask, p_E, dt, stage);

  /* if (stage == 1) { */
  /*   patch_poison_bnd(p_Unext); */
  /* } */
}

// ----------------------------------------------------------------------
// pde_mhd_pushstage

static void
pde_mhd_pushstage(struct mrc_fld *x, mrc_fld_data_t dt, int stage)
{
  fld3d_t p_f;
  fld3d_setup(&p_f, x);

  pde_for_each_patch(p) {
    fld3d_get(&p_f, p);
    patch_pushstage(p_f, dt, stage);
    fld3d_put(&p_f, p);
  }
}

