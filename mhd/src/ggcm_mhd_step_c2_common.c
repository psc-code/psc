
#include "ggcm_mhd_step_private.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_defs_extra.h"
#include "ggcm_mhd_crds.h"
#include "ggcm_mhd_diag_private.h"

#include <mrc_domain.h>
#include <mrc_profile.h>

#include <math.h>
#include <string.h>

// FIXME: major ugliness
// The fortran fields do primitive vars in the order _RR,_PP,_VX,_VY,_VZ
// but in C, we stick with the corresponding conservative var order, ie.,
// RR,VX,VY,VZ,PP
// The below hackily switches the order around in C, so that it matches fortran

#define PP 1
#define VX 2
#define VY 3
#define VZ 4

#include "pde/pde_defs.h"

// mhd options

#define OPT_EQN OPT_EQN_MHD_SCONS
#define OPT_BACKGROUND false

#include "pde/pde_mhd_compat.c"
#include "pde/pde_mhd_get_dt.c"
#include "pde/pde_mhd_pushfluid.c"
#include "pde/pde_mhd_rmaskn.c"
#include "pde/pde_mhd_calc_current.c"
#include "pde/pde_mhd_push_ej.c"
#include "pde/pde_mhd_calc_resis.c"

// FIXME, don't even know why I have to do this
#undef PP
#undef VX
#undef VY
#undef VZ

#include "mhd_sc.c"

// TODO:
// - handle various resistivity models
// - handle limit2, limit3
// - handle lowmask

#define REPS (1.e-10f)

// ======================================================================
// ggcm_mhd_step subclass "c2"
//
// this class will do full predictor / corrector steps,
// ie., including primvar() etc.

struct ggcm_mhd_step_c2 {
  struct mhd_options opt;
};

#define ggcm_mhd_step_c2(step) mrc_to_subobj(step, struct ggcm_mhd_step_c2)

static inline float
bcthy3f(mrc_fld_data_t s1, mrc_fld_data_t s2)
{
  if (s1 > 0.f && fabsf(s2) > REPS) {
/* .if(calce_aspect_low) then */
/* .call lowmask(I, 0, 0,tl1) */
/* .call lowmask( 0,J, 0,tl2) */
/* .call lowmask( 0, 0,K,tl3) */
/* .call lowmask(I,J,K,tl4) */
/*       tt=tt*(1.0-max(tl1,tl2,tl3,tl4)) */
    return s1 / s2;
  }
  return 0.f;
}

static inline void
calc_avg_dz_By(struct ggcm_mhd *mhd, struct mrc_fld *f, int m_curr, int XX, int YY, int ZZ,
	       int JX1, int JY1, int JZ1, int JX2, int JY2, int JZ2)
{
  float *bd1x = ggcm_mhd_crds_get_crd(mhd->crds, 0, BD1);
  float *bd1y = ggcm_mhd_crds_get_crd(mhd->crds, 1, BD1);
  float *bd1z = ggcm_mhd_crds_get_crd(mhd->crds, 2, BD1);

  // d_z B_y, d_y B_z on x edges
  mrc_fld_foreach(f, i,j,k, 1, 2) {
    mrc_fld_data_t bd1[3] = { bd1x[i-1], bd1y[j-1], bd1z[k-1] };

    F3(f, _TMP1, i,j,k) = bd1[ZZ] * 
      (F3(f, m_curr + _B1X + YY, i,j,k) - F3(f, m_curr + _B1X + YY, i-JX2,j-JY2,k-JZ2));
    F3(f, _TMP2, i,j,k) = bd1[YY] * 
      (F3(f, m_curr + _B1X + ZZ, i,j,k) - F3(f, m_curr + _B1X + ZZ, i-JX1,j-JY1,k-JZ1));
  } mrc_fld_foreach_end;

  // .5 * harmonic average if same sign
  mrc_fld_foreach(f, i,j,k, 1, 1) {
    mrc_fld_data_t s1, s2;
    // dz_By on y face
    s1 = F3(f, _TMP1, i+JX2,j+JY2,k+JZ2) * F3(f, _TMP1, i,j,k);
    s2 = F3(f, _TMP1, i+JX2,j+JY2,k+JZ2) + F3(f, _TMP1, i,j,k);
    F3(f, _TMP3, i,j,k) = bcthy3f(s1, s2);
    // dy_Bz on z face
    s1 = F3(f, _TMP2, i+JX1,j+JY1,k+JZ1) * F3(f, _TMP2, i,j,k);
    s2 = F3(f, _TMP2, i+JX1,j+JY1,k+JZ1) + F3(f, _TMP2, i,j,k);
    F3(f, _TMP4, i,j,k) = bcthy3f(s1, s2);
  } mrc_fld_foreach_end;
}

#define CC_TO_EC(f, m, i,j,k, I,J,K) \
  (.25f * (F3(f, m, i-I,j-J,k-K) +  \
	   F3(f, m, i-I,j   ,k   ) +  \
	   F3(f, m, i   ,j-J,k   ) +  \
	   F3(f, m, i   ,j   ,k-K)))

static inline void
calc_v_x_B(mrc_fld_data_t ttmp[2], struct mrc_fld *f, int m_curr, int i, int j, int k,
	   int XX, int YY, int ZZ, int I, int J, int K,
	   int JX1, int JY1, int JZ1, int JX2, int JY2, int JZ2,
	   float *bd2x, float *bd2y, float *bd2z, mrc_fld_data_t dt)
{
    mrc_fld_data_t bd2m[3] = { bd2x[i-1], bd2y[j-1], bd2z[k-1] };
    mrc_fld_data_t bd2[3] = { bd2x[i], bd2y[j], bd2z[k] };
    mrc_fld_data_t vbZZ;
    // edge centered velocity
    mrc_fld_data_t vvYY = CC_TO_EC(f, _VX + YY, i,j,k, I,J,K) /* - d_i * vcurrYY */;
    if (vvYY > 0.f) {
      vbZZ = F3(f, m_curr + _B1X + ZZ, i-JX1,j-JY1,k-JZ1) +
	F3(f, _TMP4, i-JX1,j-JY1,k-JZ1) * (bd2m[YY] - dt*vvYY);
    } else {
      vbZZ = F3(f, m_curr + _B1X + ZZ, i,j,k) -
	F3(f, _TMP4, i,j,k) * (bd2[YY] + dt*vvYY);
    }
    ttmp[0] = vbZZ * vvYY;

    mrc_fld_data_t vbYY;
    // edge centered velocity
    mrc_fld_data_t vvZZ = CC_TO_EC(f, _VX + ZZ, i,j,k, I,J,K) /* - d_i * vcurrZZ */;
    if (vvZZ > 0.f) {
      vbYY = F3(f, m_curr + _B1X + YY, i-JX2,j-JY2,k-JZ2) +
	F3(f, _TMP3, i-JX2,j-JY2,k-JZ2) * (bd2m[ZZ] - dt*vvZZ);
    } else {
      vbYY = F3(f, m_curr + _B1X + YY, i,j,k) -
	F3(f, _TMP3, i,j,k) * (bd2[ZZ] + dt*vvZZ);
    }
    ttmp[1] = vbYY * vvZZ;

}

static void
bcthy3z_NL1(struct ggcm_mhd *mhd, int XX, int YY, int ZZ, int I, int J, int K,
	    int JX1, int JY1, int JZ1, int JX2, int JY2, int JZ2,
	    mrc_fld_data_t dt, int m_curr)
{
  struct mrc_fld *f = mhd->fld;

  float *bd2x = ggcm_mhd_crds_get_crd(mhd->crds, 0, BD2);
  float *bd2y = ggcm_mhd_crds_get_crd(mhd->crds, 1, BD2);
  float *bd2z = ggcm_mhd_crds_get_crd(mhd->crds, 2, BD2);

  calc_avg_dz_By(mhd, f, m_curr, XX, YY, ZZ, JX1, JY1, JZ1, JX2, JY2, JZ2);

  mrc_fld_data_t diffmul=1.0;
  if (mhd->time < mhd->par.diff_timelo) { // no anomalous res at startup
    diffmul = 0.f;
  }

  // edge centered E = - v x B (+ dissipation)
  mrc_fld_foreach(f, i,j,k, 0, 1) {
    mrc_fld_data_t ttmp[2];
    calc_v_x_B(ttmp, f, m_curr, i, j, k, XX, YY, ZZ, I, J, K,
	       JX1, JY1, JZ1, JX2, JY2, JZ2, bd2x, bd2y, bd2z, dt);
    
    mrc_fld_data_t t1m = F3(f, m_curr + _B1X + ZZ, i+JX1,j+JY1,k+JZ1) - F3(f, m_curr + _B1X + ZZ, i,j,k);
    mrc_fld_data_t t1p = fabsf(F3(f, m_curr + _B1X + ZZ, i+JX1,j+JY1,k+JZ1)) + fabsf(F3(f, m_curr + _B1X + ZZ, i,j,k));
    mrc_fld_data_t t2m = F3(f, m_curr + _B1X + YY, i+JX2,j+JY2,k+JZ2) - F3(f, m_curr + _B1X + YY, i,j,k);
    mrc_fld_data_t t2p = fabsf(F3(f, m_curr + _B1X + YY, i+JX2,j+JY2,k+JZ2)) + fabsf(F3(f, m_curr + _B1X + YY, i,j,k));
    mrc_fld_data_t tp = t1p + t2p + REPS;
    mrc_fld_data_t tpi = diffmul / tp;
    mrc_fld_data_t d1 = sqr(t1m * tpi);
    mrc_fld_data_t d2 = sqr(t2m * tpi);
    if (d1 < mhd->par.diffth) d1 = 0.;
    if (d2 < mhd->par.diffth) d2 = 0.;
    ttmp[0] -= d1 * t1m * F3(f, _RMASK, i,j,k);
    ttmp[1] -= d2 * t2m * F3(f, _RMASK, i,j,k);
    F3(f, _RESIS, i,j,k) += fabsf(d1+d2) * F3(f, _ZMASK, i,j,k);
    F3(f, _FLX + XX, i,j,k) = ttmp[0] - ttmp[1];
  } mrc_fld_foreach_end;
}

static void
bcthy3z_const(struct ggcm_mhd *mhd, int XX, int YY, int ZZ, int I, int J, int K,
	      int JX1, int JY1, int JZ1, int JX2, int JY2, int JZ2, mrc_fld_data_t dt, int m_curr)
{
  struct mrc_fld *f = mhd->fld;

  float *bd2x = ggcm_mhd_crds_get_crd(mhd->crds, 0, BD2);
  float *bd2y = ggcm_mhd_crds_get_crd(mhd->crds, 1, BD2);
  float *bd2z = ggcm_mhd_crds_get_crd(mhd->crds, 2, BD2);

  calc_avg_dz_By(mhd, f, m_curr, XX, YY, ZZ, JX1, JY1, JZ1, JX2, JY2, JZ2);

  // edge centered E = - v x B (+ dissipation)
  mrc_fld_foreach(f, i,j,k, 0, 1) {
    mrc_fld_data_t ttmp[2];
    calc_v_x_B(ttmp, f, m_curr, i, j, k, XX, YY, ZZ, I, J, K,
	       JX1, JY1, JZ1, JX2, JY2, JZ2, bd2x, bd2y, bd2z, dt);

    mrc_fld_data_t vcurrXX = CC_TO_EC(f, _CURRX + XX, i,j,k, I,J,K);
    mrc_fld_data_t vresis = CC_TO_EC(f, _RESIS, i,j,k, I,J,K);
    F3(f, _FLX + XX, i,j,k) = ttmp[0] - ttmp[1] - vresis * vcurrXX;
  } mrc_fld_foreach_end;
}

static void
calce_nl1_c(struct ggcm_mhd *mhd, mrc_fld_data_t dt, int m_curr)
{
  bcthy3z_NL1(mhd, 0,1,2, 0,1,1, 0,1,0, 0,0,1, dt, m_curr);
  bcthy3z_NL1(mhd, 1,2,0, 1,0,1, 0,0,1, 1,0,0, dt, m_curr);
  bcthy3z_NL1(mhd, 2,0,1, 1,1,0, 1,0,0, 0,1,0, dt, m_curr);
}

static void
calce_const_c(struct ggcm_mhd *mhd, mrc_fld_data_t dt, int m_curr)
{
  bcthy3z_const(mhd, 0,1,2, 0,1,1, 0,1,0, 0,0,1, dt, m_curr);
  bcthy3z_const(mhd, 1,2,0, 1,0,1, 0,0,1, 1,0,0, dt, m_curr);
  bcthy3z_const(mhd, 2,0,1, 1,1,0, 1,0,0, 0,1,0, dt, m_curr);
}

static void
calce_c(struct ggcm_mhd *mhd, mrc_fld_data_t dt, int m_curr)
{
  switch (mhd->par.magdiffu) {
  case MAGDIFFU_NL1:
    return calce_nl1_c(mhd, dt, m_curr);
  case MAGDIFFU_CONST:
    return calce_const_c(mhd, dt, m_curr);
  default:
    assert(0);
  }
}

static void
bpush_c(struct ggcm_mhd *mhd, mrc_fld_data_t dt, int m_prev, int m_next)
{
  struct mrc_fld *f = mhd->fld;
  float *bd3x = ggcm_mhd_crds_get_crd(mhd->crds, 0, BD3);
  float *bd3y = ggcm_mhd_crds_get_crd(mhd->crds, 1, BD3);
  float *bd3z = ggcm_mhd_crds_get_crd(mhd->crds, 2, BD3);

  mrc_fld_foreach(f, i,j,k, 0, 0) {
    F3(f, m_next + _B1X, i,j,k) = F3(f, m_prev + _B1X, i,j,k) +
      dt * (bd3y[j] * (F3(f,_FLZ, i,j+1,k) - F3(f,_FLZ, i,j,k)) -
	    bd3z[k] * (F3(f,_FLY, i,j,k+1) - F3(f,_FLY, i,j,k)));
    F3(f, m_next + _B1Y, i,j,k) = F3(f, m_prev + _B1Y, i,j,k) +
      dt * (bd3z[k] * (F3(f,_FLX, i,j,k+1) - F3(f,_FLX, i,j,k)) -
	    bd3x[i] * (F3(f,_FLZ, i+1,j,k) - F3(f,_FLZ, i,j,k)));
    F3(f, m_next + _B1Z, i,j,k) = F3(f, m_prev + _B1Z, i,j,k) +
      dt * (bd3x[i] * (F3(f,_FLY, i+1,j,k) - F3(f,_FLY, i,j,k)) -
	    bd3y[j] * (F3(f,_FLX, i,j+1,k) - F3(f,_FLX, i,j,k)));
  } mrc_fld_foreach_end;
}

static void
pushstage_c(struct ggcm_mhd *mhd, mrc_fld_data_t dt, int m_prev, int m_curr, int m_next,
	    bool limit)
{
  struct mrc_fld *f = mhd->fld;
  fld3d_t p_f;
  fld3d_setup(&p_f, f);
  pde_patch_set(0);
  fld3d_get(&p_f, 0);

  int stage = m_curr == _RR1 ? 0 : 1;

  fld3d_t p_rmask = fld3d_make_view(p_f, _RMASK);
  fld3d_t p_resis = fld3d_make_view(p_f, _RESIS);
  fld3d_t p_Jcc = fld3d_make_view(p_f, _CURRX);
  fld3d_t p_Unext, p_Uprev, p_Ucurr;
  fld3d_t p_W, p_cmsv, p_ymask, p_zmask;
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

  patch_rmaskn(p_rmask, p_zmask);

  patch_pushfluid_c(p_Unext, dt, p_Uprev, p_Ucurr, p_W, p_cmsv,
		    p_ymask, p_zmask, stage);

  switch (mhd->par.magdiffu) {
  case MAGDIFFU_NL1:
    patch_calc_resis_nl1_c(p_resis);
    break;
  case MAGDIFFU_CONST:
    patch_calc_resis_const_c(p_resis, p_Jcc, p_Ucurr, p_zmask);
    //calc_resis_const_c(mhd, m_curr);
    break;
  default:
    assert(0);
  }

  patch_push_ej(p_Unext, dt, p_Ucurr, p_W, p_zmask);

  calce_c(mhd, dt, m_curr);
  bpush_c(mhd, dt, m_prev, m_next);
}

// ======================================================================
// ggcm_mhd_step subclass "c2"
//
// this class will do full predictor / corrector steps,
// ie., including primvar() etc.

// ----------------------------------------------------------------------
// ggcm_mhd_step_c2_get_dt

static double
ggcm_mhd_step_c2_get_dt(struct ggcm_mhd_step *step, struct mrc_fld *x)
{
  struct ggcm_mhd *mhd = step->mhd;

  ggcm_mhd_fill_ghosts(mhd, x, _RR1, mhd->time);
  zmaskn(mhd, mhd->fld, _ZMASK, x, _YMASK, mhd->fld);
  // assert(strcmp(mrc_fld_type(mhd->fld), "float") == 0);
  return pde_mhd_get_dt_scons(mhd, x, x, _ZMASK);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_c2_pred

static void
ggcm_mhd_step_c2_pred(struct ggcm_mhd_step *step)
{
  struct ggcm_mhd *mhd = step->mhd;

  primvar_c(mhd, _RR1);
  zmaskn(mhd, mhd->fld, _ZMASK, mhd->fld, _YMASK, mhd->fld);

  mrc_fld_data_t dth = .5f * step->mhd->dt;
  static int PR;
  if (!PR) {
    PR = prof_register("pred_c", 1., 0, 0);
  }
  prof_start(PR);
  pushstage_c(step->mhd, dth, _RR1, _RR1, _RR2, false);
  prof_stop(PR);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_c2_corr

static void
ggcm_mhd_step_c2_corr(struct ggcm_mhd_step *step)
{
  primvar_c(step->mhd, _RR2);
  static int PR;
  if (!PR) {
    PR = prof_register("corr_c", 1., 0, 0);
  }
  prof_start(PR);
  pushstage_c(step->mhd, step->mhd->dt, _RR1, _RR2, _RR1, true);
  prof_stop(PR);
  
  // --- check for NaNs and small density
  // (still controlled by do_badval_checks)
  badval_checks_sc(step->mhd, step->mhd->fld, step->mhd->fld);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_c2_run

static void
ggcm_mhd_step_c2_run(struct ggcm_mhd_step *step, struct mrc_fld *x)
{
  struct ggcm_mhd *mhd = step->mhd;

  assert(x == mhd->fld);

  ggcm_mhd_fill_ghosts(mhd, x, _RR1, mhd->time);
  ggcm_mhd_step_c2_pred(step);

  ggcm_mhd_fill_ghosts(mhd, x, _RR2, mhd->time + mhd->bndt);
  ggcm_mhd_step_c2_corr(step);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_c2_setup

static void
ggcm_mhd_step_c2_setup(struct ggcm_mhd_step *step)
{
  struct ggcm_mhd *mhd = step->mhd;

  pde_setup(mhd->fld);
  pde_mhd_setup(mhd);
  pde_mhd_compat_setup(mhd);

  mhd->ymask = mrc_fld_make_view(mhd->fld, _YMASK, _YMASK + 1);
  mrc_fld_set(mhd->ymask, 1.);

  ggcm_mhd_step_setup_member_objs_sub(step);
  ggcm_mhd_step_setup_super(step);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_c2_setup_flds

static void
ggcm_mhd_step_c2_setup_flds(struct ggcm_mhd_step *step)
{
  struct ggcm_mhd_step_c2 *sub = ggcm_mhd_step_c2(step);
  struct ggcm_mhd *mhd = step->mhd;

  pde_mhd_set_options(mhd, &sub->opt);
  mrc_fld_set_type(mhd->fld, FLD_TYPE);
  mrc_fld_set_param_int(mhd->fld, "nr_ghosts", 2);
  mrc_fld_dict_add_int(mhd->fld, "mhd_type", MT_SEMI_CONSERVATIVE);
  mrc_fld_set_param_int(mhd->fld, "nr_comps", _NR_FLDS);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_c2_get_e_ec

static void
ggcm_mhd_step_c2_get_e_ec(struct ggcm_mhd_step *step, struct mrc_fld *Eout,
                          struct mrc_fld *state_vec)
{
  // the state vector should already be FLD_TYPE, but Eout is the data type
  // of the output
  struct mrc_fld *E = mrc_fld_get_as(Eout, FLD_TYPE);
  struct mrc_fld *x = mrc_fld_get_as(state_vec, FLD_TYPE);

  mrc_fld_foreach(x, i, j, k, 0, 1) {
    F3(E, 0, i,j,k) = F3(x, _FLX, i,j,k);
    F3(E, 1, i,j,k) = F3(x, _FLY, i,j,k);
    F3(E, 2, i,j,k) = F3(x, _FLZ, i,j,k);
  } mrc_fld_foreach_end;

  mrc_fld_put_as(E, Eout);
  // FIXME, should use _put_as, but don't want copy-back
  if (strcmp(mrc_fld_type(state_vec), FLD_TYPE) != 0) {
    mrc_fld_destroy(x);
  }
} 

// ----------------------------------------------------------------------
// ggcm_mhd_step_c2_diag_item_zmask_run

static void
ggcm_mhd_step_c2_diag_item_zmask_run(struct ggcm_mhd_step *step,
				    struct ggcm_mhd_diag_item *item,
				    struct mrc_io *io, struct mrc_fld *f,
				    int diag_type, float plane)
{
  ggcm_mhd_diag_c_write_one_field(io, f, _ZMASK, "zmask", 1., diag_type, plane);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_c2_diag_item_rmask_run

static void
ggcm_mhd_step_c2_diag_item_rmask_run(struct ggcm_mhd_step *step,
				    struct ggcm_mhd_diag_item *item,
				    struct mrc_io *io, struct mrc_fld *f,
				    int diag_type, float plane)
{
  ggcm_mhd_diag_c_write_one_field(io, f, _RMASK, "rmask", 1., diag_type, plane);
}

// ----------------------------------------------------------------------
// subclass description

#define VAR(x) (void *)offsetof(struct ggcm_mhd_step_c2, x)
static struct param ggcm_mhd_step_c2_descr[] = {
  { "eqn"                , VAR(opt.eqn)            , PARAM_SELECT(OPT_EQN,
								  opt_eqn_descr)                },
  { "mhd_primvar"        , VAR(opt.mhd_primvar)    , PARAM_SELECT(OPT_MHD_C,
								  opt_mhd_descr)                },
  { "mhd_primbb"         , VAR(opt.mhd_primbb)     , PARAM_SELECT(OPT_MHD_C,
								  opt_mhd_descr)                },
  { "mhd_zmaskn"         , VAR(opt.mhd_zmaskn)     , PARAM_SELECT(OPT_MHD_C,
								  opt_mhd_descr)                },
  { "mhd_rmaskn"         , VAR(opt.mhd_rmaskn)     , PARAM_SELECT(OPT_MHD_C,
								  opt_mhd_descr)                },
  { "mhd_newstep"        , VAR(opt.mhd_newstep)    , PARAM_SELECT(OPT_MHD_C,
								  opt_mhd_descr)                },
  { "mhd_pushpred"       , VAR(opt.mhd_pushpred)   , PARAM_SELECT(OPT_MHD_C,
								  opt_mhd_descr)                },
  { "mhd_pushcorr"       , VAR(opt.mhd_pushcorr)   , PARAM_SELECT(OPT_MHD_C,
								  opt_mhd_descr)                },
  { "mhd_pushfluid1"     , VAR(opt.mhd_pushfluid1) , PARAM_SELECT(OPT_MHD_C,
								  opt_mhd_descr)                },
  { "mhd_pushfluid2"     , VAR(opt.mhd_pushfluid2) , PARAM_SELECT(OPT_MHD_C,
								  opt_mhd_descr)                },
  { "mhd_pushfield1"     , VAR(opt.mhd_pushfield1) , PARAM_SELECT(OPT_MHD_C,
								  opt_mhd_descr)                },
  { "mhd_pushfield2"     , VAR(opt.mhd_pushfield2) , PARAM_SELECT(OPT_MHD_C,
								  opt_mhd_descr)                },
  { "mhd_push_ej"        , VAR(opt.mhd_push_ej)    , PARAM_SELECT(OPT_MHD_C,
								  opt_mhd_descr)                },
  { "mhd_pfie3"          , VAR(opt.mhd_pfie3)      , PARAM_SELECT(OPT_MHD_C,
								  opt_mhd_descr)                },
  { "mhd_bpush1"         , VAR(opt.mhd_bpush1)     , PARAM_SELECT(OPT_MHD_C,
								  opt_mhd_descr)                },
  { "mhd_calce"          , VAR(opt.mhd_calce)      , PARAM_SELECT(OPT_MHD_C,
								  opt_mhd_descr)                },
  { "mhd_calc_resis"     , VAR(opt.mhd_calc_resis) , PARAM_SELECT(OPT_MHD_C,
								  opt_mhd_descr)                },
  
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_step subclass "c2_*"

struct ggcm_mhd_step_ops ggcm_mhd_step_c2_ops = {
  .name                = ggcm_mhd_step_c2_name,
  .size                = sizeof(struct ggcm_mhd_step_c2),
  .param_descr         = ggcm_mhd_step_c2_descr,
  .get_dt              = ggcm_mhd_step_c2_get_dt,
  .run                 = ggcm_mhd_step_c2_run,
  .setup               = ggcm_mhd_step_c2_setup,
  .setup_flds          = ggcm_mhd_step_c2_setup_flds,
  .get_e_ec            = ggcm_mhd_step_c2_get_e_ec,
  .diag_item_zmask_run = ggcm_mhd_step_c2_diag_item_zmask_run,
  .diag_item_rmask_run = ggcm_mhd_step_c2_diag_item_rmask_run,
};
