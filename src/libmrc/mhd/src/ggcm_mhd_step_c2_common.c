
#undef HAVE_OPENGGCM_FORTRAN

#include "ggcm_mhd_step_private.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_diag_private.h"

#include <string.h>

#include "pde/pde_defs.h"

#define OPT_GGCM_CRDS OPT_GGCM_CRDS_LEGACY

// mhd options

#define OPT_EQN OPT_EQN_MHD_SCONS

#include "pde/pde_mhd_compat.c"
#include "pde/pde_mhd_get_dt.c"
#include "pde/pde_mhd_rmaskn.c"
#include "pde/pde_mhd_pushfluid.c"
#include "pde/pde_mhd_push_ej.c"
#include "pde/pde_mhd_calc_resis.c"
#include "pde/pde_mhd_calce.c"
#include "pde/pde_mhd_bpush.c"
#include "pde/pde_mhd_badval_checks.c"

// TODO:
// - handle remaining resistivity models
// - handle limit2, limit3
// - handle lowmask

// ======================================================================
// ggcm_mhd_step subclass "c2"
//
// this class will do full predictor / corrector steps,
// ie., including primvar() etc.

struct ggcm_mhd_step_c2 {
  struct mhd_options opt;

  struct mrc_fld *f_zmask;
  struct mrc_fld *f_Uhalf;
  struct mrc_fld *f_E;
};

#define ggcm_mhd_step_c2(step) mrc_to_subobj(step, struct ggcm_mhd_step_c2)

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
  mrc_fld_dict_add_int(mhd->fld, "mhd_type", MT_SCONS_FC);
  mrc_fld_set_param_int(mhd->fld, "nr_comps", 8);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_c2_setup

static void
ggcm_mhd_step_c2_setup(struct ggcm_mhd_step *step)
{
  struct ggcm_mhd_step_c2 *sub = ggcm_mhd_step_c2(step);
  struct ggcm_mhd *mhd = step->mhd;

  pde_mhd_setup(mhd, mrc_fld_nr_comps(mhd->fld));
  pde_mhd_compat_setup(mhd);

  mhd->ymask = ggcm_mhd_get_3d_fld(mhd, 1);
  mrc_fld_set(mhd->ymask, 1.);

  sub->f_zmask = ggcm_mhd_get_3d_fld(mhd, 1);

  sub->f_Uhalf = ggcm_mhd_get_3d_fld(mhd, 8);
  mrc_fld_dict_add_int(sub->f_Uhalf, "mhd_type", MT_SCONS_FC);

  sub->f_E = ggcm_mhd_get_3d_fld(mhd, 3);

  ggcm_mhd_step_setup_member_objs_sub(step);
  ggcm_mhd_step_setup_super(step);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_c2_destroy

static void
ggcm_mhd_step_c2_destroy(struct ggcm_mhd_step *step)
{
  struct ggcm_mhd_step_c2 *sub = ggcm_mhd_step_c2(step);
  struct ggcm_mhd *mhd = step->mhd;

  ggcm_mhd_put_3d_fld(mhd, mhd->ymask);
  ggcm_mhd_put_3d_fld(mhd, sub->f_zmask);
  ggcm_mhd_put_3d_fld(mhd, sub->f_Uhalf);
  ggcm_mhd_put_3d_fld(mhd, sub->f_E);

  pde_free();
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_c2_get_dt

static double
ggcm_mhd_step_c2_get_dt(struct ggcm_mhd_step *step, struct mrc_fld *x)
{
  return pde_mhd_get_dt_scons(step->mhd, x, step->mhd->ymask);
}

// ----------------------------------------------------------------------
// patch_pushstage

static void
patch_pushstage(fld3d_t p_Unext, mrc_fld_data_t dt, fld3d_t p_Ucurr, fld3d_t p_ymask,
		fld3d_t p_zmask, fld3d_t p_E, int stage, int limit)
{
  static fld3d_t p_W, p_cmsv, p_rmask, p_Jcc, p_B;
  fld3d_setup_tmp_compat(&p_W, 5, _RR);
  fld3d_setup_tmp_compat(&p_cmsv, 1, _CMSV);
  fld3d_setup_tmp_compat(&p_rmask, 1, _RMASK);
  fld3d_setup_tmp_compat(&p_Jcc, 3, _CURRX);
  fld3d_setup_tmp_compat(&p_B, 3, _BX);

  patch_primvar(p_W, p_Ucurr, p_cmsv);
  patch_badval_checks_sc(p_Ucurr, p_W); // FIXME, incorporate

  if (stage == 0) {
    patch_calc_zmask(p_zmask, p_Ucurr, p_ymask);
  }

  patch_rmaskn(p_rmask, p_zmask);

  if (limit) {
    vgrs(p_B, 0, 0.f); vgrs(p_B, 1, 0.f); vgrs(p_B, 2, 0.f);
    assert(!s_do_limit2);
    assert(!s_do_limit3);
    limit1_c(p_W, PP, p_B);
  }

  pushfv_c(p_Unext, p_Unext, p_Ucurr, RR , p_W, p_cmsv, p_ymask, dt, limit, p_B);
  pushfv_c(p_Unext, p_Unext, p_Ucurr, RVX, p_W, p_cmsv, p_ymask, dt, limit, p_B);
  pushfv_c(p_Unext, p_Unext, p_Ucurr, RVY, p_W, p_cmsv, p_ymask, dt, limit, p_B);
  pushfv_c(p_Unext, p_Unext, p_Ucurr, RVZ, p_W, p_cmsv, p_ymask, dt, limit, p_B);
  pushfv_c(p_Unext, p_Unext, p_Ucurr, UU , p_W, p_cmsv, p_ymask, dt, limit, p_B);

  pushpp_c(p_Unext, p_W, p_zmask, dt);

  patch_push_ej(p_Unext, dt, p_Ucurr, p_W, p_zmask);

  patch_calc_e(p_E, dt, p_Ucurr, p_W, p_zmask, p_rmask);
  patch_bpush1(p_Unext, dt, p_Unext, p_E);
}

// ----------------------------------------------------------------------
// pushstage

static void
pushstage(struct mrc_fld *f_Unext, mrc_fld_data_t dt, struct mrc_fld *f_Ucurr,
	  struct mrc_fld *f_ymask, struct mrc_fld *f_zmask, struct mrc_fld *f_E,
	  int stage)
{
  fld3d_t p_Unext, p_Ucurr, p_ymask, p_zmask, p_E;
  fld3d_setup(&p_Unext, f_Unext);
  fld3d_setup(&p_Ucurr, f_Ucurr);
  fld3d_setup(&p_ymask, f_ymask);
  fld3d_setup(&p_zmask, f_zmask);
  fld3d_setup(&p_E    , f_E);

  bool limit = stage != 0 && s_mhd_time > s_timelo;

  pde_for_each_patch(p) {
    fld3d_t *patches[] = { &p_Unext, &p_Ucurr, &p_ymask, &p_zmask, &p_E, NULL };
    fld3d_get_list(p, patches);
    patch_pushstage(p_Unext, dt, p_Ucurr, p_ymask, p_zmask, p_E, stage, limit);
    fld3d_put_list(p, patches);
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_c2_run

static void
ggcm_mhd_step_c2_run(struct ggcm_mhd_step *step, struct mrc_fld *f_U)
{
  struct ggcm_mhd_step_c2 *sub = ggcm_mhd_step_c2(step);
  struct ggcm_mhd *mhd = step->mhd;

  struct mrc_fld *f_Uhalf = sub->f_Uhalf;
  struct mrc_fld *f_ymask = mhd->ymask, *f_zmask = sub->f_zmask, *f_E = sub->f_E;

  // FIXME? It's not going to make a difference, but this is the
  // time at the beginning of the whole step, rather than the time of the current state
  s_mhd_time = mhd->time_code * mhd->tnorm; 

  // set f_Uhalf = f_U
  mrc_fld_copy(f_Uhalf, f_U);
  // then advance f_Uhalf += .5f * dt * rhs(f_U)
  ggcm_mhd_fill_ghosts(mhd, f_U, mhd->time_code);
  pushstage(f_Uhalf, .5f * mhd->dt_code, f_U, f_ymask, f_zmask, f_E, 0);

  // f_U += dt * rhs(f_Uhalf)
  ggcm_mhd_fill_ghosts(mhd, f_Uhalf, mhd->time_code +.5f *  mhd->dt_code);
  pushstage(f_U, mhd->dt_code, f_Uhalf, f_ymask, f_zmask, f_E, 1);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_c2_get_e_ec

static void
ggcm_mhd_step_c2_get_e_ec(struct ggcm_mhd_step *step, struct mrc_fld *Eout,
                          struct mrc_fld *state_vec)
{
  struct ggcm_mhd_step_c2 *sub = ggcm_mhd_step_c2(step);

  // the state vector should already be FLD_TYPE, but Eout is the data type
  // of the output
  struct mrc_fld *E = mrc_fld_get_as(Eout, FLD_TYPE);
  struct mrc_fld *f_E = sub->f_E;

  fld3d_foreach(i, j, k, 0, 1) {
    F3(E, 0, i,j,k) = F3(f_E, 0, i,j,k);
    F3(E, 1, i,j,k) = F3(f_E, 1, i,j,k);
    F3(E, 2, i,j,k) = F3(f_E, 2, i,j,k);
  } fld3d_foreach_end;

  mrc_fld_put_as(E, Eout);
} 

// ----------------------------------------------------------------------
// ggcm_mhd_step_c2_diag_item_zmask_run

static void
ggcm_mhd_step_c2_diag_item_zmask_run(struct ggcm_mhd_step *step,
				    struct ggcm_mhd_diag_item *item,
				    struct mrc_io *io, struct mrc_fld *f,
				    int diag_type, float plane)
{
  struct ggcm_mhd_step_c2 *sub = ggcm_mhd_step_c2(step);
  ggcm_mhd_diag_c_write_one_field(io, sub->f_zmask, 0, "zmask", 1., diag_type, plane);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_c2_diag_item_rmask_run

static void
ggcm_mhd_step_c2_diag_item_rmask_run(struct ggcm_mhd_step *step,
				    struct ggcm_mhd_diag_item *item,
				    struct mrc_io *io, struct mrc_fld *f,
				    int diag_type, float plane)
{
  // rmask is only allocated temporarily and as one patch only, this
  // would have to change if we want to output it
  assert(0);
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
  .setup               = ggcm_mhd_step_c2_setup,
  .destroy             = ggcm_mhd_step_c2_destroy,
  .setup_flds          = ggcm_mhd_step_c2_setup_flds,
  .get_dt              = ggcm_mhd_step_c2_get_dt,
  .run                 = ggcm_mhd_step_c2_run,
  .get_e_ec            = ggcm_mhd_step_c2_get_e_ec,
  .diag_item_zmask_run = ggcm_mhd_step_c2_diag_item_zmask_run,
  .diag_item_rmask_run = ggcm_mhd_step_c2_diag_item_rmask_run,
};
