
// ======================================================================
// patch_pushpred

// ----------------------------------------------------------------------
// patch_pushpred_fortran

#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)

#include "pde/pde_fortran.h"

#define pushpred_F77 F77_FUNC(pushpred,PUSHPRED)
#define rmaskn_F77 F77_FUNC(rmaskn,RMASKN)
#define pushfluid1_F77 F77_FUNC(pushfluid1,PUSHFLUID1)
#define pushfield1_F77 F77_FUNC(pushfield1,PUSHFIELD1)

void pushpred_F77(real *rr1, real *rv1x, real *rv1y, real *rv1z, real *uu1,
		  real *b1x, real *b1y, real *b1z,
		  real *rr2, real *rv2x, real *rv2y, real *rv2z, real *uu2,
		  real *b2x, real *b2y, real *b2z,
		  real *rr, real *vx, real *vy, real *vz, real *pp,
		  real *cmsv, real *ymask, real *zmask, real *rmask,
		  real *flx, real *fly, real *flz,
		  real *tmp1, real *tmp2, real *tmp3, real *resis,
		  real *dt, real *time);
void rmaskn_F77(real *zmask, real *rmask);
void pushfluid1_F77(real *rr1, real *rv1x, real *rv1y, real *rv1z, real *uu1,
		    real *rr2, real *rv2x, real *rv2y, real *rv2z, real *uu2,
		    real *rr, real *vx, real *vy, real *vz, real *pp,
		    real *ymask, real *zmask, real *cmsv,
		    real *dth);
void pushfield1_F77(real *rr1, real *rv1x, real *rv1y, real *rv1z, real *uu1,
		    real *b1x, real *b1y, real *b1z,
		    real *rr2, real *rv2x, real *rv2y, real *rv2z, real *uu2,
		    real *b2x, real *b2y, real *b2z,
		    real *rr, real *vx, real *vy, real *vz, real *pp,
		    real *cmsv, real *ymask, real *zmask, real *rmask,
		    real *flx, real *fly, real *flz,
		    real *tmp1, real *tmp2, real *tmp3, real *resis,
		    real *dth, real *time);

static void
patch_rmaskn_fortran(fld3d_t p_f)
{
  rmaskn_F77(F(p_f, _ZMASK), F(p_f, _RMASK));
}

static void
patch_pushfluid1_fortran(fld3d_t p_f, mrc_fld_data_t dth)
{
  pushfluid1_F77(F(p_f, _RR1), F(p_f, _RV1X), F(p_f, _RV1Y), F(p_f, _RV1Z), F(p_f, _UU1),
		 F(p_f, _RR2), F(p_f, _RV2X), F(p_f, _RV2Y), F(p_f, _RV2Z), F(p_f, _UU2),
		 F(p_f, _RR), F(p_f, _VX), F(p_f, _VY), F(p_f, _VZ), F(p_f, _PP),
		 F(p_f, _YMASK), F(p_f, _ZMASK), F(p_f, _CMSV),
		 &dth);
}

static void
patch_pushfield1_fortran(fld3d_t p_f, mrc_fld_data_t dth)
{
  pushfield1_F77(F(p_f, _RR1), F(p_f, _RV1X), F(p_f, _RV1Y), F(p_f, _RV1Z), F(p_f, _UU1),
		 F(p_f, _B1X), F(p_f, _B1Y), F(p_f, _B1Z),
		 F(p_f, _RR2), F(p_f, _RV2X), F(p_f, _RV2Y), F(p_f, _RV2Z), F(p_f, _UU2),
		 F(p_f, _B2X), F(p_f, _B2Y), F(p_f, _B2Z),
		 F(p_f, _RR), F(p_f, _VX), F(p_f, _VY), F(p_f, _VZ), F(p_f, _PP),
		 F(p_f, _CMSV), F(p_f, _YMASK), F(p_f, _ZMASK), F(p_f, _RMASK),
		 F(p_f, _FLX), F(p_f, _FLY), F(p_f, _FLZ),
		 F(p_f, _TMP1), F(p_f, _TMP2), F(p_f, _TMP3), F(p_f, _RESIS),
		 &dth, &s_mhd_time);
}

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

#endif

// ======================================================================

// ----------------------------------------------------------------------
// patch_pushpred_c

static void
patch_pushpred_c(fld3d_t p_f, mrc_fld_data_t dt)
{
  mrc_fld_data_t dth = .5f * dt;

  patch_rmaskn_fortran(p_f);
  patch_pushfluid1_fortran(p_f, dth);
  patch_pushfield1_fortran(p_f, dth);
}

// ----------------------------------------------------------------------
// patch_pushpred

static void
patch_pushpred(fld3d_t p_f, mrc_fld_data_t dt)
{
  if (s_opt_mhd_pushpred == OPT_MHD_C) {
    patch_pushpred_c(p_f, dt);
#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)
  } else if (s_opt_mhd_pushpred == OPT_MHD_FORTRAN) {
    patch_pushpred_fortran(p_f, dt);
#endif
  } else {
    assert(0);
  }
}

// ======================================================================
// patch_pushcorr

// ----------------------------------------------------------------------
// patch_pushcorr_fortran

#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)

#include "pde/pde_fortran.h"

#define pushcorr_F77 F77_FUNC(pushcorr,PUSHCORR)

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

#endif

// ----------------------------------------------------------------------
// patch_pushcorr

static void
patch_pushcorr(fld3d_t p_f, mrc_fld_data_t dt)
{
  if (s_opt_mhd_pushcorr == OPT_MHD_C) {
    assert(0);
#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)
  } else if (s_opt_mhd_pushcorr == OPT_MHD_FORTRAN) {
    patch_pushcorr_fortran(p_f, dt);
#endif
  } else {
    assert(0);
  }
}

