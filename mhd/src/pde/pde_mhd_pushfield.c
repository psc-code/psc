
#ifndef PDE_MHD_PUSHFIELD_C
#define PDE_MHD_PUSHFIELD_C

#include "pde/pde_mhd_push_ej.c"
#include "pde/pde_mhd_calc_resis.c"
#include "pde/pde_mhd_pfie3.c"

// ----------------------------------------------------------------------
// vgr0

static void
vgr0(fld3d_t p_f, int m)
{
  fld3d_foreach(i,j,k, 2, 2) {
    F3S(p_f, m, i,j,k) = 0.;
  } fld3d_foreach_end;
}

// ----------------------------------------------------------------------
// patch_pushfield_c

static void
patch_pushfield_c(fld3d_t p_Unext, mrc_fld_data_t dt, fld3d_t p_Uprev,
		  fld3d_t p_Ucurr, fld3d_t p_W, fld3d_t p_zmask, fld3d_t p_rmask,
		  fld3d_t p_resis, fld3d_t p_Jcc, fld3d_t p_E, int stage)
{
  if (s_magdiffu == MAGDIFFU_NL1) {
    patch_calc_resis_nl1(p_resis);
    // FIXME, should be unnec
    vgr0(p_Jcc, 0);
    vgr0(p_Jcc, 1);
    vgr0(p_Jcc, 2);
  } else if (s_magdiffu == MAGDIFFU_RES1) {
    assert(0);
    //    calc_resis_res1(bxB,byB,bzB,currx,curry,currz,tmp1,tmp2,tmp3,flx,fly,flz,zmask,rr,pp,resis);
  } else if (s_magdiffu == MAGDIFFU_CONST) {
    patch_calc_resis_const(p_resis, p_Jcc, p_Ucurr, p_zmask);
  }

  patch_push_ej(p_Unext, dt, p_Ucurr, p_W, p_zmask);
  patch_pfie3(p_Unext, dt, p_Uprev, p_Ucurr, p_W, p_zmask, p_rmask,
	      p_resis, p_Jcc, p_E);
}

// ----------------------------------------------------------------------
// patch_pushfield_fortran

#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)

#define pushfield1_F77 F77_FUNC(pushfield1,PUSHFIELD1)
#define pushfield2_F77 F77_FUNC(pushfield2,PUSHFIELD2)

void pushfield1_F77(real *rr1, real *rv1x, real *rv1y, real *rv1z, real *uu1,
		    real *b1x, real *b1y, real *b1z,
		    real *rr2, real *rv2x, real *rv2y, real *rv2z, real *uu2,
		    real *b2x, real *b2y, real *b2z,
		    real *rr, real *vx, real *vy, real *vz, real *pp,
		    real *cmsv, real *ymask, real *zmask, real *rmask,
		    real *flx, real *fly, real *flz,
		    real *tmp1, real *tmp2, real *tmp3, real *resis,
		    real *dth, real *time);
void pushfield2_F77(real *rr1, real *rv1x, real *rv1y, real *rv1z, real *uu1,
		    real *b1x, real *b1y, real *b1z,
		    real *rr2, real *rv2x, real *rv2y, real *rv2z, real *uu2,
		    real *b2x, real *b2y, real *b2z,
		    real *rr, real *vx, real *vy, real *vz, real *pp,
		    real *cmsv, real *ymask, real *zmask, real *rmask,
		    real *flx, real *fly, real *flz,
		    real *tmp1, real *tmp2, real *tmp3, real *resis,
		    real *dth, real *time);

static void
patch_pushfield1_fortran(mrc_fld_data_t dt)
{
  pushfield1_F77(F(s_p_f, _RR1), F(s_p_f, _RV1X), F(s_p_f, _RV1Y), F(s_p_f, _RV1Z), F(s_p_f, _UU1),
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
patch_pushfield2_fortran(mrc_fld_data_t dt)
{
  pushfield2_F77(F(s_p_f, _RR1), F(s_p_f, _RV1X), F(s_p_f, _RV1Y), F(s_p_f, _RV1Z), F(s_p_f, _UU1),
		 F(s_p_f, _B1X), F(s_p_f, _B1Y), F(s_p_f, _B1Z),
		 F(s_p_f, _RR2), F(s_p_f, _RV2X), F(s_p_f, _RV2Y), F(s_p_f, _RV2Z), F(s_p_f, _UU2),
		 F(s_p_f, _B2X), F(s_p_f, _B2Y), F(s_p_f, _B2Z),
		 F(s_p_f, _RR), F(s_p_f, _VX), F(s_p_f, _VY), F(s_p_f, _VZ), F(s_p_f, _PP),
		 F(s_p_f, _CMSV), F(s_p_f, _YMASK), F(s_p_f, _ZMASK), F(s_p_f, _RMASK),
		 F(s_p_f, _FLX), F(s_p_f, _FLY), F(s_p_f, _FLZ),
		 F(s_p_f, _TMP1), F(s_p_f, _TMP2), F(s_p_f, _TMP3), F(s_p_f, _RESIS),
		 &dt, &s_mhd_time);
}

// ----------------------------------------------------------------------
// patch_pushfield_fortran

static void _mrc_unused
patch_pushfield_fortran(mrc_fld_data_t dt, int stage)
{
  if (stage == 0) {
    patch_pushfield1_fortran(dt);
  } else {
    patch_pushfield2_fortran(dt);
  }
}

#endif

// ----------------------------------------------------------------------
// patch_pushfield

static void _mrc_unused
patch_pushfield(fld3d_t p_Unext, mrc_fld_data_t dt, fld3d_t p_Uprev,
		fld3d_t p_Ucurr, fld3d_t p_W, fld3d_t p_zmask, fld3d_t p_rmask,
		fld3d_t p_resis, fld3d_t p_Jcc, fld3d_t p_E, int stage)
{
  int opt_mhd_pushfield = stage ? s_opt_mhd_pushfield2 : s_opt_mhd_pushfield1;

  if (opt_mhd_pushfield == OPT_MHD_C) {
    patch_pushfield_c(p_Unext, dt, p_Uprev, p_Ucurr, p_W,
		      p_zmask, p_rmask, p_resis, p_Jcc, p_E, stage);
#if defined(HAVE_OPENGGCM_FORTRAN) && defined(MRC_FLD_AS_FLOAT_H)
  } else if (opt_mhd_pushfield == OPT_MHD_FORTRAN) {
    patch_pushfield_fortran(dt, stage);
#endif
  } else {
    assert(0);
  }
}

#endif
