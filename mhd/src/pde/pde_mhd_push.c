
#include "pde/pde_mhd_rmaskn.c"
#include "pde/pde_mhd_pushfluid.c"
#include "pde/pde_mhd_pushfield.c"
#include "pde/pde_mhd_badval_checks.c"

// ----------------------------------------------------------------------
// patch_push_c

static void
patch_push_c(fld3d_t p_f, mrc_fld_data_t dt, int stage)
{
  fld3d_t p_Unext, p_Uprev, p_Ucurr;
  fld3d_t p_W, p_cmsv, p_ymask, p_zmask, p_rmask;
  fld3d_t p_resis, p_Jcc;
  if (stage == 0) {
    fld3d_setup_view(&p_Unext, p_f, _RR2);
    fld3d_setup_view(&p_Uprev, p_f, _RR1);
    fld3d_setup_view(&p_Ucurr, p_f, _RR1);
  } else {
    fld3d_setup_view(&p_Unext, p_f, _RR1);
    fld3d_setup_view(&p_Uprev, p_f, _RR1);
    fld3d_setup_view(&p_Ucurr, p_f, _RR2);
  }
  fld3d_setup_view(&p_W    , p_f, _RR);
  fld3d_setup_view(&p_cmsv , p_f, _CMSV);
  fld3d_setup_view(&p_ymask, p_f, _YMASK);
  fld3d_setup_view(&p_zmask, p_f, _ZMASK);
  fld3d_setup_view(&p_rmask, p_f, _RMASK);
  fld3d_setup_view(&p_resis, p_f, _RESIS);
  fld3d_setup_view(&p_Jcc  , p_f, _CURRX);

  if (stage == 0) {
    dt *= .5f;
  }

  patch_rmaskn(p_rmask, p_zmask);
  patch_pushfluid(p_Unext, dt, p_Uprev, p_Ucurr, p_W,
		  p_cmsv, p_ymask, p_zmask, p_f, stage);
  patch_pushfield(p_Unext, dt, p_Uprev, p_Ucurr, p_W,
		  p_zmask, p_rmask, p_resis, p_Jcc, p_f, stage);
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
patch_pushpred_fortran(fld3d_t p_f, mrc_fld_data_t dt)
{
  pushpred_F77(F(p_f, _RR1), F(p_f, _RV1X), F(p_f, _RV1Y), F(p_f, _RV1Z), F(p_f, _UU1),
	       F(p_f, _B1X), F(p_f, _B1Y), F(p_f, _B1Z),
	       F(p_f, _RR2), F(p_f, _RV2X), F(p_f, _RV2Y), F(p_f, _RV2Z), F(p_f, _UU2),
	       F(p_f, _B2X), F(p_f, _B2Y), F(p_f, _B2Z),
	       F(p_f, _RR), F(p_f, _VX), F(p_f, _VY), F(p_f, _VZ), F(p_f, _PP),
	       F(p_f, _CMSV), F(p_f, _YMASK), F(p_f, _ZMASK), F(p_f, _RMASK),
	       F(p_f, _FLX), F(p_f, _FLY), F(p_f, _FLZ),
	       F(p_f, _TMP1), F(p_f, _TMP2), F(p_f, _TMP3), F(p_f, _RESIS),
	       &dt, &s_mhd_time);
}

static void
patch_pushcorr_fortran(fld3d_t p_f, mrc_fld_data_t dt)
{
  pushcorr_F77(F(p_f, _RR1), F(p_f, _RV1X), F(p_f, _RV1Y), F(p_f, _RV1Z), F(p_f, _UU1),
	       F(p_f, _B1X), F(p_f, _B1Y), F(p_f, _B1Z),
	       F(p_f, _RR2), F(p_f, _RV2X), F(p_f, _RV2Y), F(p_f, _RV2Z), F(p_f, _UU2),
	       F(p_f, _B2X), F(p_f, _B2Y), F(p_f, _B2Z),
	       F(p_f, _RR), F(p_f, _VX), F(p_f, _VY), F(p_f, _VZ), F(p_f, _PP),
	       F(p_f, _CMSV), F(p_f, _YMASK), F(p_f, _ZMASK), F(p_f, _RMASK),
	       F(p_f, _FLX), F(p_f, _FLY), F(p_f, _FLZ),
	       F(p_f, _TMP1), F(p_f, _TMP2), F(p_f, _TMP3), F(p_f, _RESIS),
	       &dt, &s_mhd_time);
}

static void
patch_push_fortran(fld3d_t p_f, mrc_fld_data_t dt, int stage)
{
  if (stage == 0) {
    patch_pushpred_fortran(p_f, dt);
  } else {
    patch_pushcorr_fortran(p_f, dt);
  }
}

#endif

// ----------------------------------------------------------------------
// patch_push

static void
patch_push(fld3d_t p_f, mrc_fld_data_t dt, int stage)
{
  int opt_mhd_push = stage ? s_opt_mhd_pushpred : s_opt_mhd_pushcorr;

  if (opt_mhd_push == OPT_MHD_C) {
    patch_push_c(p_f, dt, stage);
#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)
  } else if (opt_mhd_push == OPT_MHD_FORTRAN) {
    patch_push_fortran(p_f, dt, stage);
#endif
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// patch_pushstage

static void
patch_pushstage(fld3d_t p_f, mrc_fld_data_t dt, int stage)
{
  fld3d_t p_W, p_Ucurr, p_cmsv, p_bcc, p_ymask, p_zmask;
  fld3d_setup_view(&p_W    , p_f, _RR);
  if (stage == 0) {
    fld3d_setup_view(&p_Ucurr, p_f, _RR1);
  } else {
    fld3d_setup_view(&p_Ucurr, p_f, _RR2);
  }
  fld3d_setup_view(&p_cmsv , p_f, _CMSV);

  patch_primvar(p_W, p_Ucurr, p_cmsv);

  if (stage == 0) {
    fld3d_setup_view(&p_bcc  , p_f, _BX);
    fld3d_setup_view(&p_ymask, p_f, _YMASK);
    fld3d_setup_view(&p_zmask, p_f, _ZMASK);
    patch_primbb(p_bcc, p_Ucurr);
    patch_zmaskn(p_zmask, p_W, p_bcc, p_ymask);
  }

  patch_push(p_f, dt, stage);

  if (stage == 1) {
    patch_badval_checks_sc(p_Ucurr, p_W);
  }
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

