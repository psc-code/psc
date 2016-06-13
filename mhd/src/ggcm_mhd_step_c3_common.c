
#include "ggcm_mhd_step_private.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_crds.h"
#include "ggcm_mhd_diag_private.h"
#include "mhd_util.h"

#include <mrc_domain.h>
#include <mrc_profile.h>
#include <mrc_io.h>

#include <math.h>
#include <string.h>

#include "pde/pde_defs.h"

static bool s_opt_bc_reconstruct = false;

// mhd options

#define OPT_EQN OPT_EQN_MHD_SCONS

#include "pde/pde_setup.c"
#include "pde/pde_mhd_setup.c"
#include "pde/pde_mhd_line.c"
#include "pde/pde_mhd_convert.c"
#include "pde/pde_mhd_reconstruct.c"
#include "pde/pde_mhd_divb_glm.c"
#include "pde/pde_mhd_riemann.c"
#include "pde/pde_mhd_stage.c"
#include "pde/pde_mhd_get_dt.c"

#include "mhd_3d.c"
#include "mhd_sc.c"

//FIXME, when using hydro_rusanov / no pushpp, things go wrong when > timelo

#define _BT(p_U, d, i,j,k)  (F3S(p_U, BX+d, i,j,k) + (s_opt_background ? F3S(p_b0, d, i,j,k) : 0))

// ======================================================================
// ggcm_mhd_step subclass "c3"

struct ggcm_mhd_step_c3 {
  struct mhd_options opt;

  struct mrc_fld *zmask;
  struct mrc_fld *rmask;
  
  struct mrc_fld *f_W;
  struct mrc_fld *f_Uhalf;
  struct mrc_fld *f_fluxes[3];
  struct mrc_fld *f_E;

  bool enforce_rrmin;
};

#define ggcm_mhd_step_c3(step) mrc_to_subobj(step, struct ggcm_mhd_step_c3)

// TODO:
// - handle various resistivity models
// - handle limit2, limit3
// - handle lowmask

#define REPS (1.e-10f)

enum {
  LIMIT_NONE,
  LIMIT_1,
};

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
    sub->f_fluxes[d] = ggcm_mhd_get_3d_fld(mhd, s_n_comps);
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
    ggcm_mhd_put_3d_fld(mhd, sub->f_fluxes[d]);
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

  mrc_fld_set_type(mhd->fld, FLD_TYPE);
  mrc_fld_set_param_int(mhd->fld, "nr_ghosts", 2);
  mrc_fld_dict_add_int(mhd->fld, "mhd_type", MT_SEMI_CONSERVATIVE);
  mrc_fld_set_param_int(mhd->fld, "nr_comps", 8);
}

// ----------------------------------------------------------------------
// patch_primvar

static void
patch_primvar(fld3d_t p_W, fld3d_t p_U, int p)
{
  int dir = 0;
  pde_for_each_line(dir, j, k, 2) {
    int ib = -2, ie = s_ldims[0] + 2;
    mhd_line_get_state(l_U, p_U, j, k, dir, ib, ie);
    mhd_prim_from_cons(l_W, l_U, ib, ie);
    mhd_line_put_state(l_W, p_W, j, k, dir, ib, ie);
  }
}

// ----------------------------------------------------------------------
// patch_zmaskn

static void
patch_zmaskn(struct ggcm_mhd *mhd, fld3d_t p_zmask, fld3d_t p_ymask,
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

// ======================================================================
// (hydro) predictor

// ----------------------------------------------------------------------
// line_flux_pred_pt1

static void
line_flux_pred_pt1(fld3d_t p_U, int j, int k, int dir, int ib, int ie)
{
  mhd_line_get_state(l_U, p_U, j, k, dir, ib - 1, ie + 1);
  mhd_prim_from_cons(l_W, l_U, ib - 1, ie + 1);
  mhd_reconstruct(l_Ul, l_Ur, l_Wl, l_Wr, l_W, (fld1d_t) {}, ib, ie + 1);
}

// ----------------------------------------------------------------------
// line_flux_pred_pt2

static void
line_flux_pred_pt2(fld3d_t p_F, int j, int k, int dir, int ib, int ie)
{
  mhd_riemann(l_F, l_Ul, l_Ur, l_Wl, l_Wr, ib, ie + 1);
  mhd_line_put_state(l_F, p_F, j, k, dir, ib, ie + 1);
}

// ----------------------------------------------------------------------
// patch_flux_pred

static void
patch_flux_pred(struct ggcm_mhd_step *step, fld3d_t p_F[3], fld3d_t p_U)
{
  pde_for_each_dir(dir) {
    pde_for_each_line(dir, j, k, 0) {
      int ib = 0, ie = s_ldims[dir];
      line_flux_pred_pt1(p_U, j, k, dir, ib, ie);
      line_flux_pred_pt2(p_F[dir], j, k, dir, ib, ie);
    }
  }
}

// ----------------------------------------------------------------------
// patch_flux_pred_bc_reconstruct
//
// does the same as patch_flux_pred, but allows for setting b.c. on
// reconstructed fields along the way

static void
patch_flux_pred_bc_reconstruct(struct ggcm_mhd_step *step, fld3d_t p_F[3], fld3d_t p_U, int p)
{
  // even though we do everything for one patch only, we still need
  // the full mrc_fld's for setting the boundary reconstructed values,
  // since there's no mechanism to pass a single-patch field
  static struct mrc_fld *f_Ul[3], *f_Ur[3];
  static fld3d_t p_Ul[3], p_Ur[3];
  if (!f_Ul[0]) {
    for (int d = 0; d < 3; d++) {
      f_Ul[d] = ggcm_mhd_get_3d_fld(step->mhd, s_n_comps);
      f_Ur[d] = ggcm_mhd_get_3d_fld(step->mhd, s_n_comps);
      fld3d_setup(&p_Ul[d]);
      fld3d_setup(&p_Ur[d]);
    }
  }

  // reconstruct
  for (int d = 0; d < 3; d++) {
    fld3d_get(&p_Ul[d], f_Ul[d], p);
    fld3d_get(&p_Ur[d], f_Ur[d], p);
  }
  
  pde_for_each_dir(dir) {
    pde_for_each_line(dir, j, k, 0) {
      int ib = 0, ie = s_ldims[dir];
      line_flux_pred_pt1(p_U, j, k, dir, ib, ie);
      mhd_line_put_state(l_Ul, p_Ul[dir], j, k, dir, ib, ie + 1);
      mhd_line_put_state(l_Ur, p_Ul[dir], j, k, dir, ib, ie + 1);
    }
  }
  
  for (int d = 0; d < 3; d++) {
    fld3d_put(&p_Ul[d], f_Ul[d], p);
    fld3d_put(&p_Ur[d], f_Ur[d], p);
  }

  // set boundary on reconstructed values
  ggcm_mhd_fill_ghosts_reconstr(step->mhd, f_Ul, f_Ur, p);
  
  // riemann solve
  for (int d = 0; d < 3; d++) {
    fld3d_get(&p_Ul[d], f_Ul[d], p);
    fld3d_get(&p_Ur[d], f_Ur[d], p);
  }

  pde_for_each_dir(dir) {
    pde_for_each_line(dir, j, k, 0) {
      int ib = 0, ie = s_ldims[dir];
      mhd_line_get_state(l_Ul, p_Ul[dir], j, k, dir, ib, ie + 1);
      mhd_line_get_state(l_Ur, p_Ur[dir], j, k, dir, ib, ie + 1);
      mhd_prim_from_cons(l_Wl, l_Ul, ib, ie + 1);
      mhd_prim_from_cons(l_Wr, l_Ur, ib, ie + 1);
      line_flux_pred_pt2(p_F[dir], j, k, dir, ib, ie);
    }
  }

  for (int d = 0; d < 3; d++) {
    fld3d_get(&p_Ul[d], f_Ul[d], p);
    fld3d_get(&p_Ur[d], f_Ur[d], p);
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
patch_push_pp(fld3d_t p_U, mrc_fld_data_t dt, fld3d_t pW, fld3d_t pzmask)
{
  mrc_fld_data_t dth = -.5f * dt;

  fld3d_foreach(i,j,k, 0, 0) {
    mrc_fld_data_t z = dth * F3S(pzmask, 0, i,j,k);
    F3S(p_U, RVX, i,j,k) += z * PDE_INV_DX(i) * (F3S(pW, PP, i+di,j,k) - F3S(pW, PP, i-di,j,k));
    F3S(p_U, RVY, i,j,k) += z * PDE_INV_DY(j) * (F3S(pW, PP, i,j+dj,k) - F3S(pW, PP, i,j-dj,k));
    F3S(p_U, RVZ, i,j,k) += z * PDE_INV_DZ(k) * (F3S(pW, PP, i,j,k+dk) - F3S(pW, PP, i,j,k-dk));
  } fld3d_foreach_end;
}

// ----------------------------------------------------------------------
// patch_calc_current_ec
//
// edge centered current density

static void
patch_calc_current_ec(fld3d_t j_ec, fld3d_t x)
{
  fld3d_foreach(i,j,k, 1, 2) {
    F3S(j_ec, 0, i,j,k) =
      (F3S(x, BZ, i,j,k) - F3S(x, BZ, i,j-dj,k)) * PDE_INV_DYF(j) -
      (F3S(x, BY, i,j,k) - F3S(x, BY, i,j,k-dk)) * PDE_INV_DZF(k);
    F3S(j_ec, 1, i,j,k) =
      (F3S(x, BX, i,j,k) - F3S(x, BX, i,j,k-dk)) * PDE_INV_DZF(k) -
      (F3S(x, BZ, i,j,k) - F3S(x, BZ, i-di,j,k)) * PDE_INV_DXF(i);
    F3S(j_ec, 2, i,j,k) =
      (F3S(x, BY, i,j,k) - F3S(x, BY, i-di,j,k)) * PDE_INV_DXF(i) -
      (F3S(x, BX, i,j,k) - F3S(x, BX, i,j-dj,k)) * PDE_INV_DYF(j);
  } fld3d_foreach_end;
}

// ----------------------------------------------------------------------
// patch_calc_current_cc
//
// cell-centered current density

static void
patch_calc_current_cc(fld3d_t j_cc, fld3d_t x, fld3d_t zmask)
{ 
  static fld3d_t j_ec;
  if (!fld3d_is_setup(j_ec)) {
    fld3d_setup_tmp(&j_ec, 3);
  }
  
  // get j on edges
  patch_calc_current_ec(j_ec, x);
  
  // then average to cell centers
  fld3d_foreach(i,j,k, 2, 1) {
    mrc_fld_data_t s = .25f * F3S(zmask, 0,  i,j,k);
    F3S(j_cc, 0, i,j,k) = s * (F3S(j_ec, 0, i   ,j+dj,k+dk) + F3S(j_ec, 0, i   ,j   ,k+dk) +
			       F3S(j_ec, 0, i   ,j+dj,k   ) + F3S(j_ec, 0, i   ,j   ,k   ));
    F3S(j_cc, 1, i,j,k) = s * (F3S(j_ec, 1, i+di,j   ,k+dk) + F3S(j_ec, 1, i+di,j   ,k   ) +
			       F3S(j_ec, 1, i   ,j   ,k+dk) + F3S(j_ec, 1, i   ,j   ,k   ));
    F3S(j_cc, 2, i,j,k) = s * (F3S(j_ec, 2, i+di,j+dj,k   ) + F3S(j_ec, 2, i   ,j+dj,k   ) + 
			       F3S(j_ec, 2, i+di,j   ,k   ) + F3S(j_ec, 2, i   ,j   ,k   ));
  } fld3d_foreach_end;
}

// ----------------------------------------------------------------------
// patch_calc_Bt_cc
//
// cell-averaged Btotal (ie., add B0 back in, if applicable)

static void _mrc_unused
patch_calc_Bt_cc(fld3d_t p_b, fld3d_t p_U, fld3d_t p_b0, int l, int r)
{
  fld3d_foreach(i,j,k, l, r) {
    F3S(p_b, 0, i,j,k) = .5f * (_BT(p_U, 0, i,j,k) + _BT(p_U, 0, i+di,j,k));
    F3S(p_b, 1, i,j,k) = .5f * (_BT(p_U, 1, i,j,k) + _BT(p_U, 1, i,j+dj,k));
    F3S(p_b, 2, i,j,k) = .5f * (_BT(p_U, 2, i,j,k) + _BT(p_U, 2, i,j,k+dk));
  } fld3d_foreach_end;
}

// ----------------------------------------------------------------------
// patch_push_ej

static void
patch_push_ej(fld3d_t x_next, mrc_fld_data_t dt, fld3d_t x_curr, fld3d_t prim,
	      fld3d_t zmask, fld3d_t b0)
{
  static fld3d_t j_ec, b_cc;
  if (!fld3d_is_setup(j_ec)) {
    fld3d_setup_tmp(&j_ec, 3);
    fld3d_setup_tmp(&b_cc, 3);
  }

  // FIXME/OPT, cell centered current is calculated here, and later again in calce()
  patch_calc_current_ec(j_ec, x_curr);
  patch_calc_Bt_cc(b_cc, x_curr, b0, 1, 1);

  mrc_fld_data_t s1 = .25f * dt;
  fld3d_foreach(i,j,k, 0, 0) {
    mrc_fld_data_t z = F3S(zmask, 0, i,j,k);
    mrc_fld_data_t s2 = s1 * z;
    mrc_fld_data_t cx = (F3S(j_ec, 0, i   ,j+dj,k+dk) + F3S(j_ec, 0, i  ,j   ,k+dk) +
			 F3S(j_ec, 0, i   ,j+dj,k   ) + F3S(j_ec, 0, i  ,j   ,k   ));
    mrc_fld_data_t cy = (F3S(j_ec, 1, i+di,j   ,k+dk) + F3S(j_ec, 1, i  ,j   ,k+dk) +
			 F3S(j_ec, 1, i+di,j   ,k   ) + F3S(j_ec, 1, i  ,j   ,k   ));
    mrc_fld_data_t cz = (F3S(j_ec, 2, i+di,j+dj,k   ) + F3S(j_ec, 2, i  ,j+dj,k   ) +
			 F3S(j_ec, 2, i+di,j   ,k   ) + F3S(j_ec, 2, i  ,j   ,k   ));
    mrc_fld_data_t ffx = s2 * (cy * F3S(b_cc, 2, i,j,k) - cz * F3S(b_cc, 1, i,j,k));
    mrc_fld_data_t ffy = s2 * (cz * F3S(b_cc, 0, i,j,k) - cx * F3S(b_cc, 2, i,j,k));
    mrc_fld_data_t ffz = s2 * (cx * F3S(b_cc, 1, i,j,k) - cy * F3S(b_cc, 0, i,j,k));
    mrc_fld_data_t duu = (ffx * F3S(prim, VX, i,j,k) +
			  ffy * F3S(prim, VY, i,j,k) +
			  ffz * F3S(prim, VZ, i,j,k));
    
    F3S(x_next, RVX, i,j,k) += ffx;
    F3S(x_next, RVY, i,j,k) += ffy;
    F3S(x_next, RVZ, i,j,k) += ffz;
    F3S(x_next, UU , i,j,k) += duu;
  } fld3d_foreach_end;
}

// ----------------------------------------------------------------------
// patch_rmaskn

static void
patch_rmaskn(struct ggcm_mhd *mhd, fld3d_t rmask, fld3d_t zmask, int p)
{
  mrc_fld_data_t diffco = mhd->par.diffco;
  mrc_fld_data_t diff_swbnd = mhd->par.diff_swbnd;
  int diff_obnd = mhd->par.diff_obnd;
  int gdims[3];
  mrc_domain_get_global_dims(mhd->domain, gdims);

  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);

  struct mrc_patch_info info;
  mrc_domain_get_local_patch_info(mhd->domain, p, &info);

  // _ZMASK not set at -2 ghost
  fld3d_foreach(i,j,k, 1, 2) {
    F3S(rmask, 0, i,j,k) = 0.f;
    mrc_fld_data_t xxx = MRC_MCRDX(crds, i, p);
    if (xxx < diff_swbnd)
      continue;
    if (j + info.off[1] < diff_obnd)
      continue;
    if (k + info.off[2] < diff_obnd)
      continue;
    if (i + info.off[0] >= gdims[0] - diff_obnd)
      continue;
    if (j + info.off[1] >= gdims[1] - diff_obnd)
      continue;
    if (k + info.off[2] >= gdims[2] - diff_obnd)
      continue;

    F3S(rmask, 0, i,j,k) = diffco * F3S(zmask, 0, i,j,k);
  } fld3d_foreach_end;
}

// ----------------------------------------------------------------------
// patch_res1_const

static void
patch_res1_const(struct ggcm_mhd *mhd, fld3d_t resis, int p)
{
  // resistivity comes in ohm*m
  int diff_obnd = mhd->par.diff_obnd;
  mrc_fld_data_t eta0i = 1. / mhd->resnorm;
  mrc_fld_data_t diffsphere2 = sqr(mhd->par.diffsphere);
  mrc_fld_data_t diff = mhd->par.diffco * eta0i;

  int gdims[3];
  mrc_domain_get_global_dims(mhd->domain, gdims);
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);

  struct mrc_patch_info info;
  mrc_domain_get_local_patch_info(mhd->domain, p, &info);
  
  fld3d_foreach(i,j,k, 1, 1) {
    F3S(resis, 0, i,j,k) = 0.f;
    mrc_fld_data_t r2 = (sqr(MRC_MCRDX(crds, i, p)) +
			 sqr(MRC_MCRDY(crds, j, p)) +
			 sqr(MRC_MCRDZ(crds, k, p)));
    if (r2 < diffsphere2)
      continue;
    if (j + info.off[1] < diff_obnd)
      continue;
    if (k + info.off[2] < diff_obnd)
      continue;
    if (i + info.off[0] >= gdims[0] - diff_obnd)
      continue;
    if (j + info.off[1] >= gdims[1] - diff_obnd)
      continue;
    if (k + info.off[2] >= gdims[2] - diff_obnd)
      continue;
    
    F3S(resis, 0, i,j,k) = diff;
  } fld3d_foreach_end;
}

static inline mrc_fld_data_t
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
patch_calc_avg_dz_By(fld3d_t tmp, fld3d_t x, fld3d_t p_b0,
		     int XX, int YY, int ZZ,
		     int JX1, int JY1, int JZ1, int JX2, int JY2, int JZ2)
{
  fld3d_foreach(i,j,k, 1, 2) {
    // FIXME, check offset -1
    mrc_fld_data_t bd1[3] = { PDE_INV_DX(i-1), PDE_INV_DY(j-1), PDE_INV_DZ(k-1) };
    
    F3S(tmp, 0, i,j,k) = bd1[ZZ] * 
      (_BT(x, YY, i,j,k) - _BT(x, YY, i-JX2,j-JY2,k-JZ2));
    F3S(tmp, 1, i,j,k) = bd1[YY] * 
      (_BT(x, ZZ, i,j,k) - _BT(x, ZZ, i-JX1,j-JY1,k-JZ1));
  } fld3d_foreach_end;
  
  // .5 * harmonic average if same sign
  fld3d_foreach(i,j,k, 1, 1) {
    mrc_fld_data_t s1, s2;
    // dz_By on y face
    s1 = F3S(tmp, 0, i+JX2,j+JY2,k+JZ2) * F3S(tmp, 0, i,j,k);
    s2 = F3S(tmp, 0, i+JX2,j+JY2,k+JZ2) + F3S(tmp, 0, i,j,k);
    F3S(tmp, 2, i,j,k) = bcthy3f(s1, s2);
    // dy_Bz on z face
    s1 = F3S(tmp, 1, i+JX1,j+JY1,k+JZ1) * F3S(tmp, 1, i,j,k);
    s2 = F3S(tmp, 1, i+JX1,j+JY1,k+JZ1) + F3S(tmp, 1, i,j,k);
    F3S(tmp, 3, i,j,k) = bcthy3f(s1, s2);
  } fld3d_foreach_end;
}

#define _CC_TO_EC(f, m, i,j,k, I,J,K)			\
  ({							\
    (.25f * (F3S(f, m, i-di*I,j-dj*J,k-dk*K) +		\
	     F3S(f, m, i-di*I,j     ,k     ) +		\
	     F3S(f, m, i     ,j-dj*J,k     ) +		\
	     F3S(f, m, i     ,j     ,k-dk*K)));})

// ve = v - d_i J
static inline void
calc_ve_x_B(mrc_fld_data_t ttmp[2], fld3d_t x, fld3d_t prim,
	    fld3d_t curr, fld3d_t tmp, fld3d_t p_b0,
	    int i, int j, int k,
	    int XX, int YY, int ZZ, int I, int J, int K,
	    int JX1, int JY1, int JZ1, int JX2, int JY2, int JZ2,
	    mrc_fld_data_t dt)
{
  mrc_fld_data_t vcurrYY = _CC_TO_EC(curr, YY, i, j, k, I, J, K);
  mrc_fld_data_t vcurrZZ = _CC_TO_EC(curr, ZZ, i, j, k, I, J, K);
  
  // FIXME, need to check index/shift
  mrc_fld_data_t bd2m[3] = { PDE_DX(i-1), PDE_DY(j-1), PDE_DZ(k-1) };
  mrc_fld_data_t bd2[3] = { PDE_DX(i), PDE_DY(j), PDE_DZ(k) };
  mrc_fld_data_t vbZZ;
  // edge centered velocity
  mrc_fld_data_t vvYY = _CC_TO_EC(prim, VX + YY, i,j,k, I,J,K) - s_d_i * vcurrYY;
  if (vvYY > 0.f) {
    vbZZ = _BT(x, ZZ, i-JX1,j-JY1,k-JZ1) +
      F3S(tmp, 3, i-JX1,j-JY1,k-JZ1) * (bd2m[YY] - dt*vvYY);
  } else {
    vbZZ = _BT(x, ZZ, i,j,k) -
      F3S(tmp, 3, i,j,k) * (bd2[YY] + dt*vvYY);
  }
  ttmp[0] = vbZZ * vvYY;
  
  mrc_fld_data_t vbYY;
  // edge centered velocity
  mrc_fld_data_t vvZZ = _CC_TO_EC(prim, VX + ZZ, i,j,k, I,J,K) - s_d_i * vcurrZZ;
  if (vvZZ > 0.f) {
    vbYY = _BT(x, YY, i-JX2,j-JY2,k-JZ2) +
      F3S(tmp, 2, i-JX2,j-JY2,k-JZ2) * (bd2m[ZZ] - dt*vvZZ);
  } else {
    vbYY = _BT(x, YY, i,j,k) -
      F3S(tmp, 2, i,j,k) * (bd2[ZZ] + dt*vvZZ);
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
  static fld3d_t tmp;
  if (!fld3d_is_setup(tmp)) {
    fld3d_setup_tmp(&tmp, 4);
  }

  mrc_fld_data_t diffmul = 1.f;
  if (mhd->time < mhd->par.diff_timelo) { // no anomalous res at startup
    diffmul = 0.f;
  }

  // average dz_By
  patch_calc_avg_dz_By(tmp, x, p_b0, XX, YY, ZZ, JX1, JY1, JZ1, JX2, JY2, JZ2);

  // edge centered E = - ve x B (+ dissipation)
  
  fld3d_foreach(i,j,k, 0, 1) {
    mrc_fld_data_t ttmp[2];
    calc_ve_x_B(ttmp, x, prim, curr, tmp, p_b0, i, j, k, XX, YY, ZZ, I, J, K,
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

static inline void
patch_bcthy3z_const(int XX, int YY, int ZZ, int I, int J, int K,
		    int JX1, int JY1, int JZ1, int JX2, int JY2, int JZ2,
		    fld3d_t E, mrc_fld_data_t dt, fld3d_t x, fld3d_t prim,
		    fld3d_t curr, fld3d_t resis, fld3d_t b0)
{
  static fld3d_t tmp;
  if (!fld3d_is_setup(tmp)) {
    fld3d_setup_tmp(&tmp, 4);
  }
  
  // average dz_By
  patch_calc_avg_dz_By(tmp, x, b0, XX, YY, ZZ, JX1, JY1, JZ1, JX2, JY2, JZ2);
  
  // edge centered E = - ve x B (+ dissipation)
  fld3d_foreach(i,j,k, 0, 1) {
    mrc_fld_data_t ttmp[2];
    calc_ve_x_B(ttmp, x, prim, curr, tmp, b0, i, j, k, XX, YY, ZZ, I, J, K,
		JX1, JY1, JZ1, JX2, JY2, JZ2, dt);
    
    mrc_fld_data_t vcurrXX = _CC_TO_EC(curr, XX, i,j,k, I,J,K);
    mrc_fld_data_t vresis = _CC_TO_EC(resis, 0, i,j,k, I,J,K);
    F3S(E, XX, i,j,k) = - (ttmp[0] - ttmp[1]) + vresis * vcurrXX;
  } fld3d_foreach_end;
}

static void
patch_calce(struct ggcm_mhd_step *step, fld3d_t E, mrc_fld_data_t dt,
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
    patch_res1_const(mhd, resis, p);

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
// patch_update_ct

static void
patch_update_ct(struct ggcm_mhd *mhd, fld3d_t x, fld3d_t E,
		mrc_fld_data_t dt, int p)
{
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);

  mrc_fld_data_t r_db_dt_sq = sqr(mhd->par.r_db_dt);

  fld3d_foreach(i,j,k, 0, 1) {
    float crd_fc[3];
    mrc_crds_at_fc(crds, i,j,k, p, 0, crd_fc);
    if (sqr(crd_fc[0]) + sqr(crd_fc[1]) + sqr(crd_fc[2]) >= r_db_dt_sq) {
      F3S(x, BX, i,j,k) -= dt * (PDE_INV_DY(j) * (F3S(E, 2, i,j+dj,k) - F3S(E, 2, i,j,k)) -
				 PDE_INV_DZ(k) * (F3S(E, 1, i,j,k+dk) - F3S(E, 1, i,j,k)));
    }
    mrc_crds_at_fc(crds, i,j,k, p, 1, crd_fc);
    if (sqr(crd_fc[0]) + sqr(crd_fc[1]) + sqr(crd_fc[2]) >= r_db_dt_sq) {
      F3S(x, BY, i,j,k) -= dt * (PDE_INV_DZ(k) * (F3S(E, 0, i,j,k+dk) - F3S(E, 0, i,j,k)) -
				 PDE_INV_DX(i) * (F3S(E, 2, i+di,j,k) - F3S(E, 2, i,j,k)));
    }
    mrc_crds_at_fc(crds, i,j,k, p, 2, crd_fc);
    if (sqr(crd_fc[0]) + sqr(crd_fc[1]) + sqr(crd_fc[2]) >= r_db_dt_sq) {
      F3S(x, BZ, i,j,k) -= dt * (PDE_INV_DX(i) * (F3S(E, 1, i+di,j,k) - F3S(E, 1, i,j,k)) -
				 PDE_INV_DY(j) * (F3S(E, 0, i,j+dj,k) - F3S(E, 0, i,j,k)));
    }
  } fld3d_foreach_end;
}

// ----------------------------------------------------------------------
// pushstage_c

static void
pushstage_c(struct ggcm_mhd_step *step, struct mrc_fld *x_next,
	    mrc_fld_data_t dt, struct mrc_fld *f_Ucurr, struct mrc_fld *f_Wcurr,
	    int limit)
{
  struct ggcm_mhd_step_c3 *sub = ggcm_mhd_step_c3(step);
  struct ggcm_mhd *mhd = step->mhd;

  fld3d_t _x_next, p_Ucurr, p_Wcurr, ymask, zmask, rmask, p_b0, fluxes[3];
  fld3d_t E;
  fld3d_setup(&_x_next);
  fld3d_setup(&p_Ucurr);
  fld3d_setup(&p_Wcurr);
  fld3d_setup(&ymask);
  fld3d_setup(&zmask);
  fld3d_setup(&p_b0);
  for (int d = 0; d < 3; d++) {
    fld3d_setup(&fluxes[d]);
  }
  fld3d_setup(&E);

  for (int p = 0; p < mrc_fld_nr_patches(f_Ucurr); p++) {
    pde_patch_set(p);
    fld3d_get(&p_Ucurr, f_Ucurr, p);
    for (int d = 0; d < 3; d++) {
      fld3d_get(&fluxes[d], sub->f_fluxes[d], p);
    }

    if (limit == LIMIT_NONE || mhd->time < mhd->par.timelo) {
      if (s_opt_bc_reconstruct) {
	patch_flux_pred_bc_reconstruct(step, fluxes, p_Ucurr, p);
      } else {
	patch_flux_pred(step, fluxes, p_Ucurr);
      }
    } else { // !LIMIT_NONE
      patch_flux_corr(step, fluxes, p_Ucurr);
    }

    fld3d_put(&p_Ucurr, f_Ucurr, p);
    for (int d = 0; d < 3; d++) {
      fld3d_put(&fluxes[d], sub->f_fluxes[d], p);
    }
  }

  ggcm_mhd_correct_fluxes(mhd, sub->f_fluxes);

  for (int p = 0; p < mrc_fld_nr_patches(f_Ucurr); p++) {
    pde_patch_set(p);

    fld3d_get(&p_Ucurr, f_Ucurr, p);
    fld3d_get(&_x_next, x_next, p);
    fld3d_get(&p_Wcurr, f_Wcurr, p);
    fld3d_get(&ymask, mhd->ymask, p);
    fld3d_get(&zmask, sub->zmask, p);
    fld3d_get(&rmask, sub->rmask, p);
    for (int d = 0; d < 3; d++) {
      fld3d_get(&fluxes[d], sub->f_fluxes[d], p);
    }
    if (s_opt_background) {
      fld3d_get(&p_b0, mhd->b0, p);
    }
    fld3d_get(&E, sub->f_E, p);

    mhd_update_finite_volume(mhd, _x_next, fluxes, ymask, dt, 0, 0);
    patch_push_pp(_x_next, dt, p_Wcurr, zmask);
    patch_push_ej(_x_next, dt, p_Ucurr, p_Wcurr, zmask, p_b0);

    patch_rmaskn(mhd, rmask, zmask, p);
    patch_calce(step, E, dt, p_Ucurr, p_Wcurr, zmask, rmask, p_b0, p);

    fld3d_put(&p_Ucurr, f_Ucurr, p);
    fld3d_put(&_x_next, x_next, p);
    fld3d_put(&ymask, mhd->ymask, p);
    fld3d_put(&zmask, sub->zmask, p);
    fld3d_put(&rmask, sub->rmask, p);
    fld3d_put(&p_Wcurr, f_Wcurr, p);
    for (int d = 0; d < 3; d++) {
      fld3d_put(&fluxes[d], sub->f_fluxes[d], p);
    }
    if (s_opt_background) {
      fld3d_put(&p_b0, mhd->b0, p);
    }
    fld3d_put(&E, sub->f_E, p);
  }

  //  ggcm_mhd_fill_ghosts_E(mhd, E);
  if (mhd->amr > 0) {
    mrc_ddc_amr_apply(mhd->ddc_amr_E, sub->f_E);
  }

  for (int p = 0; p < mrc_fld_nr_patches(x_next); p++) {
    pde_patch_set(p);
    fld3d_get(&_x_next, x_next, p);
    fld3d_get(&E, sub->f_E, p);

    patch_update_ct(mhd, _x_next, E, dt, p);

    fld3d_put(&_x_next, x_next, p);
    fld3d_put(&E, sub->f_E, p);
  }
}

static double
ggcm_mhd_step_c3_get_dt(struct ggcm_mhd_step *step, struct mrc_fld *x)
{
  struct ggcm_mhd_step_c3 *sub = ggcm_mhd_step_c3(step);
  struct ggcm_mhd *mhd = step->mhd;
  struct mrc_fld *ymask = mhd->ymask, *zmask = sub->zmask;

  if (step->do_nwst) {
    ggcm_mhd_fill_ghosts(mhd, x, 0, mhd->time);
    fld3d_t p_zmask, p_ymask, p_U, p_b0;
    for (int p = 0; p < mrc_fld_nr_patches(zmask); p++) {
      fld3d_get(&p_zmask, zmask, p);
      fld3d_get(&p_ymask, ymask, p);
      fld3d_get(&p_U, x, p);
      fld3d_get(&p_b0, mhd->b0, p);
      patch_zmaskn(mhd, p_zmask, p_ymask, p_U, p_b0);
      fld3d_put(&p_zmask, zmask, p);
      fld3d_put(&p_ymask, ymask, p);
      fld3d_put(&p_U, x, p);
      fld3d_put(&p_b0, mhd->b0, p);
    }
    double dtn = pde_mhd_get_dt_scons(mhd, x, zmask, 0);

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
  struct mrc_fld *f_ymask = mhd->ymask, *f_zmask = sub->zmask;

  static int pr_A, pr_B;
  if (!pr_A) {
    pr_A = prof_register("c3_pred", 0, 0, 0);
    pr_B = prof_register("c3_corr", 0, 0, 0);
  }

  fld3d_t p_U, p_W, p_zmask, p_ymask, p_b0;

  // --- PREDICTOR
  prof_start(pr_A);
  ggcm_mhd_fill_ghosts(mhd, f_U, 0, mhd->time);

  // primvar
  for (int p = 0; p < mrc_fld_nr_patches(f_U); p++) {
    fld3d_get(&p_U, f_U, p);
    fld3d_get(&p_W, f_W, p);
    patch_primvar(p_W, p_U, p);
    fld3d_put(&p_U, f_U, p);
    fld3d_put(&p_W, f_W, p);
  }

  // --- check for NaNs and negative pressures
  // (still controlled by do_badval_checks)
  badval_checks_sc(mhd, f_U, f_W);

  // zmaskn
  for (int p = 0; p < mrc_fld_nr_patches(f_zmask); p++) {
    fld3d_get(&p_zmask, f_zmask, p);
    fld3d_get(&p_ymask, f_ymask, p);
    fld3d_get(&p_U, f_U, p);
    fld3d_get(&p_b0, mhd->b0, p);
    patch_zmaskn(mhd, p_zmask, p_ymask, p_U, p_b0);
    fld3d_put(&p_zmask, f_zmask, p);
    fld3d_put(&p_ymask, f_ymask, p);
    fld3d_put(&p_U, f_U, p);
    fld3d_put(&p_b0, mhd->b0, p);
  }

  // set x_half = x^n, then advance to n+1/2
  mrc_fld_copy(f_Uhalf, f_U);
  pushstage_c(step, f_Uhalf, .5f * mhd->dt, f_U, f_W, LIMIT_NONE);
  if (sub->enforce_rrmin) {
    enforce_rrmin_sc(mhd, f_Uhalf);
  }
  prof_stop(pr_A);

  // --- CORRECTOR
  prof_start(pr_B);
  ggcm_mhd_fill_ghosts(mhd, f_Uhalf, 0, mhd->time + mhd->bndt);
  for (int p = 0; p < mrc_fld_nr_patches(f_U); p++) {
    fld3d_get(&p_U, f_Uhalf, p);
    fld3d_get(&p_W, f_W, p);
    patch_primvar(p_W, p_U, p);
    fld3d_put(&p_U, f_Uhalf, p);
    fld3d_put(&p_W, f_W, p);
  }
  // --- check for NaNs and negative pressures
  // (still controlled by do_badval_checks)
  badval_checks_sc(mhd, f_Uhalf, f_W);
  pushstage_c(step, f_U, mhd->dt, f_Uhalf, f_W, LIMIT_1);
  if (sub->enforce_rrmin) {
    enforce_rrmin_sc(mhd, f_Uhalf);
  }
  prof_stop(pr_B);
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

  ggcm_mhd_fill_ghosts(mhd, f_U, 0, mhd->time);
  for (int p = 0; p < mrc_fld_nr_patches(f_E); p++) {
    pde_patch_set(p);
    fld3d_get(&p_E, f_E, p);
    fld3d_get(&p_U, f_U, p);
    fld3d_get(&p_W, f_W, p);
    fld3d_get(&p_ymask, mhd->ymask, p);
    fld3d_get(&p_zmask, sub->zmask, p);
    fld3d_get(&p_rmask, sub->rmask, p);
    fld3d_get(&p_b0, mhd->b0, p);

    patch_primvar(p_W, p_U, p);
    patch_zmaskn(mhd, p_zmask, p_ymask, p_U, p_b0);
    patch_calce(step, p_E, mhd->dt, p_U, p_W, p_zmask, p_rmask, p_b0, p);

    fld3d_put(&p_E, f_E, p);
    fld3d_put(&p_U, f_U, p);
    fld3d_put(&p_W, f_W, p);
    fld3d_put(&p_ymask, mhd->ymask, p);
    fld3d_put(&p_zmask, sub->zmask, p);
    fld3d_put(&p_rmask, sub->rmask, p);
    fld3d_put(&p_b0, mhd->b0, p);
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
