
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
#include "pde/pde_mhd_pushfluid.c"
#include "pde/pde_mhd_push_ej.c"
#include "pde/pde_mhd_rmaskn.c"
#include "pde/pde_mhd_calc_resis.c"
#include "pde/pde_mhd_calce.c"
#include "pde/pde_mhd_bpush.c"
#include "pde/pde_mhd_stage.c"
#include "pde/pde_mhd_get_dt.c"
#include "pde/pde_mhd_badval_checks.c"

//FIXME, when using hydro_rusanov / no pushpp, things go wrong when > timelo

// ======================================================================
// ggcm_mhd_step subclass "c3"

struct ggcm_mhd_step_c3 {
  struct mhd_options opt;

  struct mrc_fld *f_zmask;
  struct mrc_fld *f_rmask;
  
  struct mrc_fld *f_W;
  struct mrc_fld *f_Uhalf;
  struct mrc_fld *f_F[3];
  struct mrc_fld *f_E;
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
// ggcm_mhd_step_c3_setup_flds

static void
ggcm_mhd_step_c3_setup_flds(struct ggcm_mhd_step *step)
{
  struct ggcm_mhd_step_c3 *sub = ggcm_mhd_step_c3(step);
  struct ggcm_mhd *mhd = step->mhd;

  pde_mhd_set_options(mhd, &sub->opt);
  mrc_fld_set_type(mhd->fld, FLD_TYPE);
  mrc_fld_set_param_int(mhd->fld, "nr_ghosts", 2);
  mrc_fld_dict_add_int(mhd->fld, "mhd_type", MT_SCONS_FC);
  mrc_fld_set_param_int(mhd->fld, "nr_comps", 8);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_c3_setup

static void
ggcm_mhd_step_c3_setup(struct ggcm_mhd_step *step)
{
  struct ggcm_mhd_step_c3 *sub = ggcm_mhd_step_c3(step);
  struct ggcm_mhd *mhd = step->mhd;

  pde_mhd_setup(mhd, 5);
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

  sub->f_zmask = ggcm_mhd_get_3d_fld(mhd, 1);
  sub->f_rmask = ggcm_mhd_get_3d_fld(mhd, 1);

  sub->f_W = ggcm_mhd_get_3d_fld(mhd, s_n_comps);
  sub->f_Uhalf = ggcm_mhd_get_3d_fld(mhd, 8);
  mrc_fld_dict_add_int(sub->f_Uhalf, "mhd_type", MT_SCONS_FC);
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
  ggcm_mhd_put_3d_fld(mhd, sub->f_zmask);
  ggcm_mhd_put_3d_fld(mhd, sub->f_rmask);
  ggcm_mhd_put_3d_fld(mhd, sub->f_Uhalf);
  for (int d = 0; d < 3; d++) {
    ggcm_mhd_put_3d_fld(mhd, sub->f_F[d]);
  }
  ggcm_mhd_put_3d_fld(mhd, sub->f_E);

  pde_free();
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_c3_get_dt

static double
ggcm_mhd_step_c3_get_dt(struct ggcm_mhd_step *step, struct mrc_fld *x)
{
  struct ggcm_mhd *mhd = step->mhd;

  // FIXME, the fill_ghosts is necessary (should be in other steps, too)
  // if we do it here, we should avoid doing it again in run() -- but note
  // we're only doing it here if do_nwst is true.
  ggcm_mhd_fill_ghosts(mhd, x, mhd->time_code);
  return pde_mhd_get_dt_scons(mhd, x, mhd->ymask);
}

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
    fluxes_mhd_scons(&F1S(l_Fcc, 0, i), &F1S(l_U, 0, i), &F1S(l_W, 0, i), i);
  }

  // limit1 (Harten-Zwas)
  for (int i = ib - 1; i < ie + 1; i++) {
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
// patch_pushstage_pt1

static void
patch_pushstage_pt1(struct ggcm_mhd_step *step, fld3d_t p_Ucurr, fld3d_t p_W,
		    fld3d_t p_F[3], bool limit, int p)
{
  // primvar, badval
  patch_prim_from_cons(p_W, p_Ucurr, 2);
  patch_badval_checks_sc(p_Ucurr, p_W);
  
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
		    fld3d_t p_Ucurr, fld3d_t p_W,
		    fld3d_t p_E, fld3d_t p_F[3], fld3d_t p_ymask, fld3d_t p_zmask,
		    fld3d_t p_rmask, int stage, int p)
{
  struct ggcm_mhd *mhd = step->mhd;

  // update hydro quantities
  mhd_update_finite_volume(mhd, p_Unext, p_F, p_ymask, dt, 0, 0);

  if (stage == 0) {
    patch_calc_zmask(p_zmask, p_Ucurr, p_ymask);
  }
  // update momentum (grad p)
  pushpp_c(p_Unext, p_W, p_zmask, dt);
  // update momentum (J x B) and energy
  patch_push_ej(p_Unext, dt, p_Ucurr, p_W, p_zmask);

  // find E
  patch_rmaskn_c(p_rmask, p_zmask);
  patch_calc_e(p_E, dt, p_Ucurr, p_W, p_zmask, p_rmask);
}

// ----------------------------------------------------------------------
// pushstage

static void
pushstage(struct ggcm_mhd_step *step, struct mrc_fld *f_Unext,
	  mrc_fld_data_t dt, struct mrc_fld *f_Ucurr, struct mrc_fld *f_W,
	  int stage)
{
  struct ggcm_mhd_step_c3 *sub = ggcm_mhd_step_c3(step);
  struct ggcm_mhd *mhd = step->mhd;

  fld3d_t p_Unext, p_Ucurr, p_W, p_ymask, p_zmask, p_rmask;
  fld3d_t p_F[3], p_E;
  fld3d_setup(&p_Unext, f_Unext);
  fld3d_setup(&p_Ucurr, f_Ucurr);
  fld3d_setup(&p_W    , f_W);
  fld3d_setup(&p_ymask, mhd->ymask);
  fld3d_setup(&p_zmask, sub->f_zmask);
  fld3d_setup(&p_rmask, sub->f_rmask);
  for (int d = 0; d < 3; d++) {
    fld3d_setup(&p_F[d], sub->f_F[d]);
  }
  fld3d_setup(&p_E, sub->f_E);
  pde_mhd_p_aux_setup_b0(mhd->b0);

  bool limit = stage != 0 && s_mhd_time > s_timelo;

  // primvar, badval, reconstruct
  pde_for_each_patch(p) {
    fld3d_t *patches[] = { &p_Ucurr, &p_W, &p_F[0], &p_F[1], &p_F[2], NULL };
    fld3d_get_list(p, patches);
    patch_pushstage_pt1(step, p_Ucurr, p_W, p_F, limit, p);
    fld3d_put_list(p, patches);
  }

  // correct hydro fluxes
  ggcm_mhd_correct_fluxes(mhd, sub->f_F);

  // add MHD terms, find E
  pde_for_each_patch(p) {
    fld3d_t *mhd_patches[] = { &p_Unext, &p_Ucurr, &p_W, 
			       &p_E, &p_F[0], &p_F[1], &p_F[2],
			       &p_ymask, &p_zmask, &p_rmask, NULL };

    fld3d_get_list(p, mhd_patches);
    pde_mhd_p_aux_get(p);

    patch_pushstage_pt2(step, p_Unext, dt, p_Ucurr, p_W,
			p_E, p_F, p_ymask, p_zmask, p_rmask, stage, p);

    fld3d_put_list(p, mhd_patches);
    pde_mhd_p_aux_put(p);
  }

  // correct E field
  ggcm_mhd_correct_E(mhd, sub->f_E);

  pde_for_each_patch(p) {
    fld3d_t *update_ct_patches[] = { &p_Unext, &p_E, NULL };

    fld3d_get_list(p, update_ct_patches);
    // update B using E
    patch_bpush1(p_Unext, dt, p_Unext, p_E);
    fld3d_put_list(p, update_ct_patches);
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_c3_run

static void
ggcm_mhd_step_c3_run(struct ggcm_mhd_step *step, struct mrc_fld *f_U)
{
  struct ggcm_mhd_step_c3 *sub = ggcm_mhd_step_c3(step);
  struct ggcm_mhd *mhd = step->mhd;

  struct mrc_fld *f_Uhalf = sub->f_Uhalf, *f_W = sub->f_W;

  s_mhd_time = mhd->time_code * mhd->tnorm; 

  // set f_Uhalf = f_U
  mrc_fld_copy(f_Uhalf, f_U);
  ggcm_mhd_fill_ghosts(mhd, f_U, mhd->time_code);
  pushstage(step, f_Uhalf, .5f * mhd->dt_code, f_U, f_W, 0);

  // f_U += dt * rhs(f_Uhalf)
  ggcm_mhd_fill_ghosts(mhd, f_Uhalf, mhd->time_code + .5 * mhd->dt_code);
  pushstage(step, f_U, mhd->dt_code, f_Uhalf, f_W, 1);
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

  fld3d_t p_U, p_W, p_E, p_ymask, p_zmask, p_rmask;
  fld3d_setup(&p_U, f_U);
  fld3d_setup(&p_W, f_W);
  fld3d_setup(&p_E, f_E);
  fld3d_setup(&p_ymask, mhd->ymask);
  fld3d_setup(&p_zmask, sub->f_zmask);
  fld3d_setup(&p_rmask, sub->f_rmask);
  pde_mhd_p_aux_setup_b0(mhd->b0);

  ggcm_mhd_fill_ghosts(mhd, f_U, mhd->time_code);
  pde_for_each_patch(p) {
    fld3d_t *get_e_ec_patches[] = { &p_E, &p_U, &p_W, &p_ymask, &p_zmask, &p_rmask, NULL };
    fld3d_get_list(p, get_e_ec_patches);
    pde_mhd_p_aux_get(p);

    patch_prim_from_cons(p_W, p_U, 2);
    patch_calc_zmask(p_zmask, p_U, p_ymask);
    patch_calc_e(p_E, mhd->dt_code, p_U, p_W, p_zmask, p_rmask);

    fld3d_put_list(p, get_e_ec_patches);
    pde_mhd_p_aux_put(p);
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
  ggcm_mhd_diag_c_write_one_field(io, sub->f_zmask, 0, "zmask", 1., diag_type, plane);
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
  ggcm_mhd_diag_c_write_one_field(io, sub->f_rmask, 0, "rmask", 1., diag_type, plane);
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
  { "mhd_primvar"        , VAR(opt.mhd_primvar)    , PARAM_SELECT(OPT_MHD_C,
								  opt_mhd_descr)                },
  { "mhd_primbb"         , VAR(opt.mhd_primbb)     , PARAM_SELECT(OPT_MHD_C,
								  opt_mhd_descr)                },
  { "mhd_zmaskn"         , VAR(opt.mhd_zmaskn)     , PARAM_SELECT(OPT_MHD_C,
								  opt_mhd_descr)                },
  { "mhd_newstep"        , VAR(opt.mhd_newstep)    , PARAM_SELECT(OPT_MHD_C_V2,
								  opt_mhd_descr)                },
  { "mhd_push_ej"        , VAR(opt.mhd_push_ej)    , PARAM_SELECT(OPT_MHD_C,
								  opt_mhd_descr)                },
  { "mhd_bpush1"         , VAR(opt.mhd_bpush1)     , PARAM_SELECT(OPT_MHD_C,
								  opt_mhd_descr)                },
  
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
