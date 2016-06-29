
#undef HAVE_OPENGGCM_FORTRAN

#include "ggcm_mhd_step_private.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_diag_private.h"

#include <string.h>

#include "pde/pde_defs.h"

// mhd options

#define OPT_EQN OPT_EQN_MHD_SCONS

#include "pde/pde_mhd_compat.c"
#include "pde/pde_mhd_line.c"
#include "pde/pde_mhd_convert.c"
#include "pde/pde_mhd_reconstruct.c"
#include "pde/pde_mhd_riemann.c"
#include "pde/pde_mhd_push_ej.c"
#include "pde/pde_mhd_rmaskn.c"
#include "pde/pde_mhd_calc_resis.c"
#include "pde/pde_mhd_calce.c"
#include "pde/pde_mhd_bpush.c"
#include "pde/pde_mhd_stage.c"
#include "pde/pde_mhd_get_dt.c"
#include "pde/pde_mhd_badval_checks.c"

//FIXME, when using hydro_rusanov / no pushpp, things go wrong when > timelo

static int s_opt_enforce_rrmin;

// ======================================================================
// ggcm_mhd_step subclass "c3"

struct ggcm_mhd_step_c3 {
  struct mhd_options opt;

  struct mrc_fld *zmask;
  struct mrc_fld *rmask;
  
  struct mrc_fld *f_W;
  struct mrc_fld *f_Uhalf;
  struct mrc_fld *f_F[3];
  struct mrc_fld *f_E;

  bool enforce_rrmin;
};

#define ggcm_mhd_step_c3(step) mrc_to_subobj(step, struct ggcm_mhd_step_c3)

// TODO:
// - handle various resistivity models
// - handle limit2, limit3
// - handle lowmask

#define REPS (1.e-10f)

// 1-d state vars statically rather than having to pass them around

static fld1d_state_t l_U, l_Ul, l_Ur, l_W, l_Wl, l_Wr, l_F;

// ----------------------------------------------------------------------
// ggcm_mhd_step_c3_setup

static void
ggcm_mhd_step_c3_setup(struct ggcm_mhd_step *step)
{
  struct ggcm_mhd_step_c3 *sub = ggcm_mhd_step_c3(step);
  struct ggcm_mhd *mhd = step->mhd;

  assert(mhd);

  pde_setup(mhd->fld);
  // FIXME, very hacky way of making the 1d state fields 5-component
  s_n_comps = 5;
  pde_mhd_setup(mhd);
  pde_mhd_compat_setup(mhd);

  fld1d_state_setup(&l_U);
  fld1d_state_setup(&l_Ul);
  fld1d_state_setup(&l_Ur);
  fld1d_state_setup(&l_W);
  fld1d_state_setup(&l_Wl);
  fld1d_state_setup(&l_Wr);
  fld1d_state_setup(&l_F);

  if (s_opt_background) {
    mhd->b0 = ggcm_mhd_get_3d_fld(mhd, 3);
  }
  mhd->ymask = ggcm_mhd_get_3d_fld(mhd, 1);
  mrc_fld_set(mhd->ymask, 1.);
  sub->zmask = ggcm_mhd_get_3d_fld(mhd, 1);
  sub->rmask = ggcm_mhd_get_3d_fld(mhd, 1);

  sub->f_W = ggcm_mhd_get_3d_fld(mhd, s_n_comps);
  sub->f_Uhalf = ggcm_mhd_get_3d_fld(mhd, 8);
  mrc_fld_dict_add_int(sub->f_Uhalf, "mhd_type", MT_SEMI_CONSERVATIVE);
  for (int d = 0; d < 3; d++) {
    sub->f_F[d] = ggcm_mhd_get_3d_fld(mhd, s_n_comps);
  }
  sub->f_E = ggcm_mhd_get_3d_fld(mhd, 3);

  ggcm_mhd_step_setup_member_objs_sub(step);
  ggcm_mhd_step_setup_super(step);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_c3_destroy

static void
ggcm_mhd_step_c3_destroy(struct ggcm_mhd_step *step)
{
  struct ggcm_mhd_step_c3 *sub = ggcm_mhd_step_c3(step);
  struct ggcm_mhd *mhd = step->mhd;

  ggcm_mhd_put_3d_fld(mhd, mhd->ymask);
  ggcm_mhd_put_3d_fld(mhd, sub->zmask);
  ggcm_mhd_put_3d_fld(mhd, sub->rmask);

  ggcm_mhd_put_3d_fld(mhd, sub->f_Uhalf);
  for (int d = 0; d < 3; d++) {
    ggcm_mhd_put_3d_fld(mhd, sub->f_F[d]);
  }
  ggcm_mhd_put_3d_fld(mhd, sub->f_E);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_c3_setup_flds

static void
ggcm_mhd_step_c3_setup_flds(struct ggcm_mhd_step *step)
{
  struct ggcm_mhd_step_c3 *sub = ggcm_mhd_step_c3(step);
  struct ggcm_mhd *mhd = step->mhd;

  pde_mhd_set_options(mhd, &sub->opt);
  s_opt_enforce_rrmin = sub->enforce_rrmin;

  mrc_fld_set_type(mhd->fld, FLD_TYPE);
  mrc_fld_set_param_int(mhd->fld, "nr_ghosts", 2);
  mrc_fld_dict_add_int(mhd->fld, "mhd_type", MT_SEMI_CONSERVATIVE);
  mrc_fld_set_param_int(mhd->fld, "nr_comps", 8);
}

// ----------------------------------------------------------------------
// patch_zmaskn_x

#define _BT(p_U, d, i,j,k)  (F3S(p_U, BX+d, i,j,k) + (s_opt_background ? F3S(p_b0, d, i,j,k) : 0))

static void
patch_zmaskn_x(struct ggcm_mhd *mhd, fld3d_t p_zmask, fld3d_t p_ymask,
	       fld3d_t p_U, fld3d_t p_b0)
{
  mrc_fld_data_t va02i = 1.f / sqr(mhd->par.speedlimit / mhd->vvnorm);
  mrc_fld_data_t eps   = 1e-15f;

  fld3d_foreach(ix,iy,iz, 1, 1) {
    mrc_fld_data_t bb = (sqr(.5f * (_BT(p_U, 0, ix,iy,iz) + _BT(p_U, 0, ix+di,iy,iz))) +
			 sqr(.5f * (_BT(p_U, 1, ix,iy,iz) + _BT(p_U, 1, ix,iy+dj,iz))) +
			 sqr(.5f * (_BT(p_U, 2, ix,iy,iz) + _BT(p_U, 2, ix,iy,iz+dk))));
    mrc_fld_data_t rrm = mrc_fld_max(eps, bb * va02i);
    F3S(p_zmask, 0, ix,iy,iz) = F3S(p_ymask, 0, ix,iy,iz) *
      mrc_fld_min(1.f, F3S(p_U, RR, ix,iy,iz) / rrm);
  } fld3d_foreach_end;
}

#undef _BT

// ======================================================================
// (hydro) predictor

// ----------------------------------------------------------------------
// patch_flux_pred

static void
patch_flux_pred(struct ggcm_mhd_step *step, fld3d_t p_F[3], fld3d_t p_U)
{
  pde_for_each_dir(dir) {
    pde_for_each_line(dir, j, k, 0) {
      int ib = 0, ie = s_ldims[dir];
      mhd_line_get_state(l_U, p_U, j, k, dir, ib - 1, ie + 1);
      mhd_prim_from_cons(l_W, l_U, ib - 1, ie + 1);
      mhd_reconstruct(l_Ul, l_Ur, l_Wl, l_Wr, l_W, (fld1d_t) {}, ib, ie + 1);
      mhd_riemann(l_F, l_Ul, l_Ur, l_Wl, l_Wr, ib, ie + 1);
      mhd_line_put_state(l_F, p_F[dir], j, k, dir, ib, ie + 1);
    }
  }
}

// ======================================================================
// (hydro) corrector

// ----------------------------------------------------------------------
// limit_hz

static inline mrc_fld_data_t
limit_hz(mrc_fld_data_t a2, mrc_fld_data_t aa, mrc_fld_data_t a1)
{
  const mrc_fld_data_t reps = 0.003;
  const mrc_fld_data_t seps = -0.001;
  const mrc_fld_data_t teps = 1.e-25;

  // Harten/Zwas type switch
  mrc_fld_data_t d1 = aa - a2;
  mrc_fld_data_t d2 = a1 - aa;
  mrc_fld_data_t s1 = fabsf(d1);
  mrc_fld_data_t s2 = fabsf(d2);
  mrc_fld_data_t f1 = fabsf(a1) + fabsf(a2) + fabsf(aa);
  mrc_fld_data_t s5 = s1 + s2 + reps*f1 + teps;
  mrc_fld_data_t r3 = fabsf(s1 - s2) / s5; // edge condition
  mrc_fld_data_t f2 = seps * f1 * f1;
  if (d1 * d2 < f2) {
    r3 = 1.f;
  }
  r3 = r3 * r3;
  r3 = r3 * r3;
  r3 = fminf(2.f * r3, 1.);

  return r3;
}

// ----------------------------------------------------------------------
// line_flux_corr

static void
line_flux_corr(fld3d_t p_F, fld3d_t p_U,
	       int j, int k, int dir, int ib, int ie)
{
  static fld1d_state_t l_Fcc, l_Flo, l_lim1;
  if (!fld1d_state_is_setup(l_Fcc)) {
    fld1d_state_setup(&l_Fcc);
    fld1d_state_setup(&l_Flo);
    fld1d_state_setup(&l_lim1);
  }

  // calculate low order fluxes
  mhd_line_get_state(l_U, p_U, j, k, dir, ib - 2, ie + 2);
  mhd_prim_from_cons(l_W, l_U, ib - 2, ie + 2);
  mhd_reconstruct(l_Ul, l_Ur, l_Wl, l_Wr, l_W, (fld1d_t) {}, ib, ie + 1);
  mhd_riemann(l_Flo, l_Ul, l_Ur, l_Wl, l_Wr, ib, ie + 1);

  // find cell centered fluxes
  for (int i = ib - 2; i < ie + 2; i++) {
    fluxes_mhd_scons(&F1S(l_Flo, 0, i), &F1S(l_U, 0, i), &F1S(l_W, 0, i), i);
  }

  // limit1 (Harten-Zwas)
  for (int i = ib; i < ie; i++) {
    mrc_fld_data_t lim1_pp = limit_hz(F1S(l_W, PP, i-1), F1S(l_W, PP, i), F1S(l_W, PP, i+1));
    for (int m = 0; m < 5; m++) {
      F1S(l_lim1, m, i) = mrc_fld_max(limit_hz(F1S(l_U, m, i-1), F1S(l_U, m, i), F1S(l_U, m, i+1)), 
				      lim1_pp);
    }
  }

  // find high-order flux and blend
  for (int i = ib; i < ie + 1; i++) {
    for (int m = 0; m < 5; m++) {
      mrc_fld_data_t fhx = ((7.f / 12.f) * (F1S(l_Fcc, m, i-1) + F1S(l_Fcc, m, i  )) -
			    (1.f / 12.f) * (F1S(l_Fcc, m, i-2) + F1S(l_Fcc, m, i+1)));
      mrc_fld_data_t cx = mrc_fld_max(F1S(l_lim1, m, i-1), F1S(l_lim1, m, i));
      F1S(l_F, m, i) = cx * F1S(l_Flo, m, i) + (1.f - cx) * fhx;
    }
  }

  mhd_line_put_state(l_F, p_F, j, k, dir, ib, ie + 1);
}

// ----------------------------------------------------------------------
// patch_flux_corr

static void
patch_flux_corr(struct ggcm_mhd_step *step, fld3d_t p_F[3], fld3d_t p_U)
{
  pde_for_each_dir(dir) {
    pde_for_each_line(dir, j, k, 0) {
      int ib = 0, ie = s_ldims[dir];
      line_flux_corr(p_F[dir], p_U, j, k, dir, ib, ie);
    }
  }
}

// ----------------------------------------------------------------------
// patch_push_pp

static void
patch_push_pp(fld3d_t p_U, mrc_fld_data_t dt, fld3d_t p_W, fld3d_t p_zmask)
{
  mrc_fld_data_t dth = -.5f * dt;

  fld3d_foreach(i,j,k, 0, 0) {
    mrc_fld_data_t z = dth * F3S(p_zmask, 0, i,j,k);
    F3S(p_U, RVX, i,j,k) += z * PDE_INV_DX(i) * (F3S(p_W, PP, i+di,j,k) - F3S(p_W, PP, i-di,j,k));
    F3S(p_U, RVY, i,j,k) += z * PDE_INV_DY(j) * (F3S(p_W, PP, i,j+dj,k) - F3S(p_W, PP, i,j-dj,k));
    F3S(p_U, RVZ, i,j,k) += z * PDE_INV_DZ(k) * (F3S(p_W, PP, i,j,k+dk) - F3S(p_W, PP, i,j,k-dk));
  } fld3d_foreach_end;
}

#define _BT(p_U, d, i,j,k)  (F3S(p_U, BX+d, i,j,k) + (s_opt_background ? F3S(p_b0, d, i,j,k) : 0))

// ve = v - d_i J
static inline void
calc_ve_x_B(mrc_fld_data_t ttmp[2], fld3d_t x, fld3d_t prim,
	    fld3d_t curr, fld3d_t p_dB, fld3d_t p_b0,
	    int i, int j, int k,
	    int XX, int YY, int ZZ, int I, int J, int K,
	    int JX1, int JY1, int JZ1, int JX2, int JY2, int JZ2,
	    mrc_fld_data_t dt)
{
  mrc_fld_data_t vcurrYY = CC_TO_EC(curr, YY, i, j, k, XX);
  mrc_fld_data_t vcurrZZ = CC_TO_EC(curr, ZZ, i, j, k, XX);
  
  // FIXME, need to check index/shift
  mrc_fld_data_t bd2m[3] = { PDE_DX(i-1), PDE_DY(j-1), PDE_DZ(k-1) };
  mrc_fld_data_t bd2[3] = { PDE_DX(i), PDE_DY(j), PDE_DZ(k) };
  mrc_fld_data_t vbZZ;
  // edge centered velocity
  mrc_fld_data_t vvYY = CC_TO_EC(prim, VX + YY, i,j,k, XX) - s_d_i * vcurrYY;
  if (vvYY > 0.f) {
    vbZZ = _BT(x, ZZ, i-JX1,j-JY1,k-JZ1) +
      F3S(p_dB, 1, i-JX1,j-JY1,k-JZ1) * (bd2m[YY] - dt*vvYY);
  } else {
    vbZZ = _BT(x, ZZ, i,j,k) -
      F3S(p_dB, 1, i,j,k) * (bd2[YY] + dt*vvYY);
  }
  ttmp[0] = vbZZ * vvYY;
  
  mrc_fld_data_t vbYY;
  // edge centered velocity
  mrc_fld_data_t vvZZ = CC_TO_EC(prim, VX + ZZ, i,j,k, XX) - s_d_i * vcurrZZ;
  if (vvZZ > 0.f) {
    vbYY = _BT(x, YY, i-JX2,j-JY2,k-JZ2) +
      F3S(p_dB, 0, i-JX2,j-JY2,k-JZ2) * (bd2m[ZZ] - dt*vvZZ);
  } else {
    vbYY = _BT(x, YY, i,j,k) -
      F3S(p_dB, 0, i,j,k) * (bd2[ZZ] + dt*vvZZ);
  }
  ttmp[1] = vbYY * vvZZ;
}

static inline void
patch_bcthy3z_NL1(struct ggcm_mhd_step *step, int XX, int YY, int ZZ, int I, int J, int K,
		  int JX1, int JY1, int JZ1, int JX2, int JY2, int JZ2,
		  fld3d_t E, mrc_fld_data_t dt, fld3d_t x, fld3d_t prim,
		  fld3d_t curr, fld3d_t rmask, fld3d_t p_b0)
{
  struct ggcm_mhd *mhd = step->mhd;
  static fld3d_t p_dB;
  if (!fld3d_is_setup(p_dB)) {
    fld3d_setup_tmp(&p_dB, 2);
  }
  fld3d_t p_B = fld3d_make_view(x, BX);

  mrc_fld_data_t diffmul = 1.f;
  if (mhd->time < mhd->par.diff_timelo) { // no anomalous res at startup
    diffmul = 0.f;
  }

  // average dz_By
  patch_calc_avg_dz_By(p_dB, p_B, p_b0, XX, YY, ZZ, JX1, JY1, JZ1, JX2, JY2, JZ2);

  // edge centered E = - ve x B (+ dissipation)
  
  fld3d_foreach(i,j,k, 0, 1) {
    mrc_fld_data_t ttmp[2];
    calc_ve_x_B(ttmp, x, prim, curr, p_dB, p_b0, i, j, k, XX, YY, ZZ, I, J, K,
		JX1, JY1, JZ1, JX2, JY2, JZ2, dt);
    
    mrc_fld_data_t t1m = _BT(x, ZZ, i+JX1,j+JY1,k+JZ1) - _BT(x, ZZ, i,j,k);
    mrc_fld_data_t t1p = fabsf(_BT(x, ZZ, i+JX1,j+JY1,k+JZ1)) + fabsf(_BT(x, ZZ, i,j,k));
    mrc_fld_data_t t2m = _BT(x, YY, i+JX2,j+JY2,k+JZ2) - _BT(x, YY, i,j,k);
    mrc_fld_data_t t2p = fabsf(_BT(x, YY, i+JX2,j+JY2,k+JZ2)) + fabsf(_BT(x, YY, i,j,k));
    mrc_fld_data_t tp = t1p + t2p + REPS;
    mrc_fld_data_t tpi = diffmul / tp;
    mrc_fld_data_t d1 = sqr(t1m * tpi);
    mrc_fld_data_t d2 = sqr(t2m * tpi);
    if (d1 < mhd->par.diffth) d1 = 0.;
    if (d2 < mhd->par.diffth) d2 = 0.;
    ttmp[0] -= d1 * t1m * F3S(rmask, 0, i,j,k);
    ttmp[1] -= d2 * t2m * F3S(rmask, 0, i,j,k);
    //    M3(f, _RESIS, i,j,k, p) += fabsf(d1+d2) * ZMASK(zmask, i,j,k, p);
    F3S(E, XX, i,j,k) = - (ttmp[0] - ttmp[1]);
  } fld3d_foreach_end;
}

#undef _BT

static inline void
patch_bcthy3z_const(int XX, int YY, int ZZ, int I, int J, int K,
		    int JX1, int JY1, int JZ1, int JX2, int JY2, int JZ2,
		    fld3d_t E, mrc_fld_data_t dt, fld3d_t x, fld3d_t prim,
		    fld3d_t curr, fld3d_t resis, fld3d_t b0)
{
  static fld3d_t p_dB;
  if (!fld3d_is_setup(p_dB)) {
    fld3d_setup_tmp(&p_dB, 2);
  }
  fld3d_t p_B = fld3d_make_view(x, BX);
  
  // average dz_By
  patch_calc_avg_dz_By(p_dB, p_B, b0, XX, YY, ZZ, JX1, JY1, JZ1, JX2, JY2, JZ2);
  
  // edge centered E = - ve x B (+ dissipation)
  fld3d_foreach(i,j,k, 0, 1) {
    mrc_fld_data_t ttmp[2];
    calc_ve_x_B(ttmp, x, prim, curr, p_dB, b0, i, j, k, XX, YY, ZZ, I, J, K,
		JX1, JY1, JZ1, JX2, JY2, JZ2, dt);
    
    mrc_fld_data_t vcurrXX = CC_TO_EC(curr, XX, i,j,k, XX);
    mrc_fld_data_t vresis = CC_TO_EC(resis, 0, i,j,k, XX);
    F3S(E, XX, i,j,k) = - (ttmp[0] - ttmp[1]) + vresis * vcurrXX;
  } fld3d_foreach_end;
}

static void
xpatch_calce(struct ggcm_mhd_step *step, fld3d_t E, mrc_fld_data_t dt,
	    fld3d_t x, fld3d_t prim, fld3d_t zmask, fld3d_t rmask, fld3d_t b0,
	    int p)
{
  struct ggcm_mhd *mhd = step->mhd;
  static fld3d_t curr;
  if (!fld3d_is_setup(curr)) {
    fld3d_setup_tmp(&curr, 3);
  }

  switch (mhd->par.magdiffu) {
  case MAGDIFFU_NL1:
    patch_calc_current_cc(curr, x, zmask);

    patch_bcthy3z_NL1(step, 0,1,2, 0,dj,dk, 0,dj,0, 0,0,dk, E, dt, x, prim, curr, rmask, b0);
    patch_bcthy3z_NL1(step, 1,2,0, di,0,dk, 0,0,dk, di,0,0, E, dt, x, prim, curr, rmask, b0);
    patch_bcthy3z_NL1(step, 2,0,1, di,di,0, di,0,0, 0,dj,0, E, dt, x, prim, curr, rmask, b0);
    break;
    
  case MAGDIFFU_CONST: {
    static fld3d_t resis;
    if (!fld3d_is_setup(resis)) {
      fld3d_setup_tmp(&resis, 1);
    }

    patch_calc_current_cc(curr, x, zmask);
    patch_res1_const(resis);

    patch_bcthy3z_const(0,1,2, 0,dj,dk, 0,dj,0, 0,0,dk, E, dt, x, prim, curr, resis, b0);
    patch_bcthy3z_const(1,2,0, di,0,dk, 0,0,dk, di,0,0, E, dt, x, prim, curr, resis, b0);
    patch_bcthy3z_const(2,0,1, di,dj,0, di,0,0, 0,dj,0, E, dt, x, prim, curr, resis, b0);
    break;
  }    
  default:
    assert(0);
  }
}

// ----------------------------------------------------------------------
// patch_enforce_rrmin_sc
//
// nudge rr and uu such that rr >= rrmin if needed

static void
patch_enforce_rrmin_sc(struct ggcm_mhd *mhd, fld3d_t p_U, int p)
{
  if (!s_opt_enforce_rrmin) {
    return;
  }

  mrc_fld_data_t rrmin = mhd->par.rrmin / mhd->rrnorm;
  mrc_fld_data_t gamma_m1 = s_gamma - 1.f;
  
  fld3d_foreach(i,j,k, 0, 0) {
    mrc_fld_data_t rr = F3S(p_U, RR, i,j,k);
    if (rr < rrmin) {
      // get pressure
      mrc_fld_data_t rrvv = (sqr(F3S(p_U, RVX, i,j,k)) +
			     sqr(F3S(p_U, RVY, i,j,k)) +
			     sqr(F3S(p_U, RVZ, i,j,k)));
      mrc_fld_data_t uu = F3S(p_U, UU, i,j,k);
      mrc_fld_data_t pp = gamma_m1 * (uu - .5f * rrvv / rr);
      
      // set new values using rrmin
      mrc_fld_data_t new_rr = rrmin;
      mrc_fld_data_t new_uu = pp / gamma_m1 + .5f * rrvv / new_rr;
      F3S(p_U, RR, i,j,k) = new_rr;
      F3S(p_U, UU, i,j,k) = new_uu;
      
      mprintf("!! Note: enforcing min density at (x=%g y=%g z=%g): "
	      "rr %lg -> %lg, uu %lg -> %lg\n",
	      PDE_CRDX_CC(i), PDE_CRDY_CC(j), PDE_CRDZ_CC(k), rr, new_rr, uu, new_uu);
    }
  } fld3d_foreach_end;
}

// ----------------------------------------------------------------------
// patch_pushstage_pt1

static void
patch_pushstage_pt1(struct ggcm_mhd_step *step, fld3d_t p_Ucurr, fld3d_t p_Wcurr,
		    fld3d_t p_F[3], bool limit, int p)
{
  // primvar, badval
  patch_prim_from_cons(p_Wcurr, p_Ucurr, 2);
  patch_badval_checks_sc(p_Ucurr, p_Wcurr);
  
  // find hydro fluxes
  // FIXME: we could use the fact that we calculate primitive variables already
  if (limit) {
    patch_flux_corr(step, p_F, p_Ucurr);
  } else {
    patch_flux_pred(step, p_F, p_Ucurr);
  }
}

// ----------------------------------------------------------------------
// patch_pushstage_pt2

static void
patch_pushstage_pt2(struct ggcm_mhd_step *step, fld3d_t p_Unext, mrc_fld_data_t dt,
		    fld3d_t p_Ucurr, fld3d_t p_Wcurr,
		    fld3d_t p_E, fld3d_t p_F[3], fld3d_t p_ymask, fld3d_t p_zmask,
		    fld3d_t p_rmask, fld3d_t p_b0, int stage, int p)
{
  struct ggcm_mhd *mhd = step->mhd;

  // update hydro quantities
  mhd_update_finite_volume(mhd, p_Unext, p_F, p_ymask, dt, 0, 0);
  // update momentum (grad p)
  patch_push_pp(p_Unext, dt, p_Wcurr, p_zmask);
  if (stage == 0) {
    patch_zmaskn_x(mhd, p_zmask, p_ymask, p_Ucurr, p_b0);
  }
  // update momentum (J x B) and energy
  patch_push_ej_b0(p_Unext, dt, p_Ucurr, p_Wcurr, p_zmask, p_b0);
  // enforce rrmin
  patch_enforce_rrmin_sc(mhd, p_Unext, p);

  // find E
  patch_rmaskn_c(p_rmask, p_zmask);
  xpatch_calce(step, p_E, dt, p_Ucurr, p_Wcurr, p_zmask, p_rmask, p_b0, p);
}

// ----------------------------------------------------------------------
// pushstage

static void
pushstage(struct ggcm_mhd_step *step, struct mrc_fld *f_Unext,
	  mrc_fld_data_t dt, struct mrc_fld *f_Ucurr, struct mrc_fld *f_Wcurr,
	  int stage, bool limit)
{
  struct ggcm_mhd_step_c3 *sub = ggcm_mhd_step_c3(step);
  struct ggcm_mhd *mhd = step->mhd;

  fld3d_t p_Unext, p_Ucurr, p_Wcurr, p_ymask, p_zmask, p_rmask, p_b0;
  fld3d_t p_F[3], p_E;
  fld3d_setup(&p_Unext, f_Unext);
  fld3d_setup(&p_Ucurr, f_Ucurr);
  fld3d_setup(&p_Wcurr, f_Wcurr);
  fld3d_setup(&p_ymask, mhd->ymask);
  fld3d_setup(&p_zmask, sub->zmask);
  fld3d_setup(&p_rmask, sub->rmask);
  fld3d_setup(&p_b0, mhd->b0);
  for (int d = 0; d < 3; d++) {
    fld3d_setup(&p_F[d], sub->f_F[d]);
  }
  fld3d_setup(&p_E, sub->f_E);

  // primvar, badval, reconstruct
  pde_for_each_patch(p) {
    fld3d_t *patches[] = { &p_Ucurr, &p_Wcurr, &p_F[0], &p_F[1], &p_F[2], NULL };
    fld3d_get_list(p, patches);
    patch_pushstage_pt1(step, p_Ucurr, p_Wcurr, p_F, limit, p);
    fld3d_put_list(p, patches);
  }

  // correct hydro fluxes
  ggcm_mhd_correct_fluxes(mhd, sub->f_F);

  // add MHD terms, find E
  pde_for_each_patch(p) {
    fld3d_t *mhd_patches[] = { &p_Unext, &p_Ucurr, &p_Wcurr, 
			       &p_E, &p_F[0], &p_F[1], &p_F[2],
			       &p_ymask, &p_zmask, &p_rmask, NULL };

    fld3d_get_list(p, mhd_patches);
    if (s_opt_background) {
      fld3d_get(&p_b0, p);
    }

    patch_pushstage_pt2(step, p_Unext, dt, p_Ucurr, p_Wcurr,
			p_E, p_F, p_ymask, p_zmask, p_rmask, p_b0, stage, p);

    fld3d_put_list(p, mhd_patches);
    if (s_opt_background) {
      fld3d_put(&p_b0, p);
    }
  }

  // correct E field
  ggcm_mhd_correct_E(mhd, sub->f_E);

  pde_for_each_patch(p) {
    fld3d_t *update_ct_patches[] = { &p_Unext, &p_E, NULL };

    fld3d_get_list(p, update_ct_patches);
    // update B using E
    patch_update_ct(p_Unext, dt, p_E);
    fld3d_put_list(p, update_ct_patches);
  }
}

static double
ggcm_mhd_step_c3_get_dt(struct ggcm_mhd_step *step, struct mrc_fld *x)
{
  struct ggcm_mhd_step_c3 *sub = ggcm_mhd_step_c3(step);
  struct ggcm_mhd *mhd = step->mhd;

  if (step->do_nwst) {
    ggcm_mhd_fill_ghosts(mhd, x, 0, mhd->time);
    fld3d_t p_zmask, p_ymask, p_U, p_b0;
    fld3d_setup(&p_zmask, sub->zmask);
    fld3d_setup(&p_ymask, mhd->ymask);
    fld3d_setup(&p_U, x);
    fld3d_setup(&p_b0, mhd->b0);

    pde_for_each_patch(p) {
      fld3d_t *zmaskn_patches[] = { &p_zmask, &p_ymask, &p_U, NULL };

      fld3d_get_list(p, zmaskn_patches);
      if (s_opt_background) {
	fld3d_get(&p_b0, p);
      }

      patch_zmaskn_x(mhd, p_zmask, p_ymask, p_U, p_b0);

      fld3d_put_list(p, zmaskn_patches);
      if (s_opt_background) {
	fld3d_put(&p_b0, p);
      }
    }
    double dtn = pde_mhd_get_dt_scons_v2(mhd, x, sub->zmask, 0);

    // --- update timestep
    dtn = mrc_fld_min(1., dtn); // FIXME, only kept for compatibility

    if (dtn > 1.02 * mhd->dt || dtn < mhd->dt / 1.01) {
      mpi_printf(ggcm_mhd_comm(mhd), "switched dt %g <- %g\n", dtn, mhd->dt);

      if (mhd->istep > 0 &&
          (dtn < 0.5 * mhd->dt || dtn > 2.0 * mhd->dt)) {            
        mpi_printf(ggcm_mhd_comm(mhd), "!!! dt changed by > a factor of 2. "
                   "Dying now!\n");
        ggcm_mhd_wrongful_death(mhd, mhd->fld, 2);
      }
      
      return dtn;
    }
  }

  return mhd->dt;
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_c3_run

static void
ggcm_mhd_step_c3_run(struct ggcm_mhd_step *step, struct mrc_fld *f_U)
{
  struct ggcm_mhd_step_c3 *sub = ggcm_mhd_step_c3(step);
  struct ggcm_mhd *mhd = step->mhd;
  struct mrc_fld *f_Uhalf = sub->f_Uhalf, *f_W = sub->f_W;

  s_mhd_time = mhd->time; 

  // --- PREDICTOR
  // set U_half = U^n, then advance to n+1/2.
  // WARNING: If we're fixing up neg pressure/density in primvar, we should probably do
  // the copy later
  ggcm_mhd_fill_ghosts(mhd, f_U, 0, mhd->time);
  mrc_fld_copy(f_Uhalf, f_U);
  pushstage(step, f_Uhalf, .5f * mhd->dt, f_U, f_W, 0, false);

  // --- CORRECTOR
  ggcm_mhd_fill_ghosts(mhd, f_Uhalf, 0, mhd->time + mhd->bndt);
  int limit = mhd->time >= mhd->par.timelo;
  pushstage(step, f_U, mhd->dt, f_Uhalf, f_W, 1, limit);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_c3_get_e_ec

static void
ggcm_mhd_step_c3_get_e_ec(struct ggcm_mhd_step *step, struct mrc_fld *Eout,
                          struct mrc_fld *f_U)
{
  struct ggcm_mhd_step_c3 *sub = ggcm_mhd_step_c3(step);
  struct ggcm_mhd *mhd = step->mhd;
  struct mrc_fld *f_W = sub->f_W;
  // the state vector should already be FLD_TYPE, but Eout is the data type
  // of the output
  struct mrc_fld *f_E = mrc_fld_get_as(Eout, FLD_TYPE);

  fld3d_t p_U, p_W, p_E, p_ymask, p_zmask, p_rmask, p_b0;
  fld3d_setup(&p_U, f_U);
  fld3d_setup(&p_W, f_W);
  fld3d_setup(&p_E, f_E);
  fld3d_setup(&p_ymask, mhd->ymask);
  fld3d_setup(&p_zmask, sub->zmask);
  fld3d_setup(&p_rmask, sub->rmask);
  fld3d_setup(&p_b0, mhd->b0);

  ggcm_mhd_fill_ghosts(mhd, f_U, 0, mhd->time);
  pde_for_each_patch(p) {
    fld3d_t *get_e_ec_patches[] = { &p_E, &p_U, &p_W, &p_ymask, &p_zmask, &p_rmask, NULL };
    fld3d_get_list(p, get_e_ec_patches);
    if (s_opt_background) {
      fld3d_get(&p_b0, p);
    }

    patch_prim_from_cons(p_W, p_U, 2);
    patch_zmaskn_x(mhd, p_zmask, p_ymask, p_U, p_b0); // FIXME, name conflict
    patch_calce(step, p_E, mhd->dt, p_U, p_W, p_zmask, p_rmask, p_b0, p);

    fld3d_put_list(p, get_e_ec_patches);
    if (s_opt_background) {
      fld3d_put(&p_b0, p);
    }
  }
  //  ggcm_mhd_fill_ghosts_E(mhd, E);
  
  mrc_fld_put_as(f_E, Eout);
} 

// ----------------------------------------------------------------------
// ggcm_mhd_step_c3_diag_item_zmask_run

static void
ggcm_mhd_step_c3_diag_item_zmask_run(struct ggcm_mhd_step *step,
				    struct ggcm_mhd_diag_item *item,
				    struct mrc_io *io, struct mrc_fld *f,
				    int diag_type, float plane)
{
  struct ggcm_mhd_step_c3 *sub = ggcm_mhd_step_c3(step);
  ggcm_mhd_diag_c_write_one_field(io, sub->zmask, 0, "zmask", 1., diag_type, plane);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_c3_diag_item_rmask_run

static void
ggcm_mhd_step_c3_diag_item_rmask_run(struct ggcm_mhd_step *step,
				    struct ggcm_mhd_diag_item *item,
				    struct mrc_io *io, struct mrc_fld *f,
				    int diag_type, float plane)
{
  struct ggcm_mhd_step_c3 *sub = ggcm_mhd_step_c3(step);
  ggcm_mhd_diag_c_write_one_field(io, sub->rmask, 0, "rmask", 1., diag_type, plane);
}

// ----------------------------------------------------------------------
// subclass description

#define VAR(x) (void *)offsetof(struct ggcm_mhd_step_c3, x)
static struct param ggcm_mhd_step_c3_descr[] = {
  { "eqn"                , VAR(opt.eqn)            , PARAM_SELECT(OPT_EQN,
								  opt_eqn_descr)                },
  { "limiter"            , VAR(opt.limiter)        , PARAM_SELECT(OPT_LIMITER_FLAT,
								  opt_limiter_descr)            },
  { "riemann"            , VAR(opt.riemann)        , PARAM_SELECT(OPT_RIEMANN_RUSANOV,
								  opt_riemann_descr)            },
  { "background"         , VAR(opt.background)     , PARAM_BOOL(false)                          },
  { "limiter_mc_beta"    , VAR(opt.limiter_mc_beta), PARAM_DOUBLE(2.)                           },

  { "enforce_rrmin"      , VAR(enforce_rrmin)      , PARAM_BOOL(false)                          },
  
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_step subclass "c3_*"

struct ggcm_mhd_step_ops ggcm_mhd_step_c3_ops = {
  .name                = ggcm_mhd_step_c3_name,
  .size                = sizeof(struct ggcm_mhd_step_c3),
  .param_descr         = ggcm_mhd_step_c3_descr,
  .setup               = ggcm_mhd_step_c3_setup,
  .get_dt              = ggcm_mhd_step_c3_get_dt,
  .run                 = ggcm_mhd_step_c3_run,
  .destroy             = ggcm_mhd_step_c3_destroy,
  .setup_flds          = ggcm_mhd_step_c3_setup_flds,
  .get_e_ec            = ggcm_mhd_step_c3_get_e_ec,
  .diag_item_zmask_run = ggcm_mhd_step_c3_diag_item_zmask_run,
  .diag_item_rmask_run = ggcm_mhd_step_c3_diag_item_rmask_run,
};
