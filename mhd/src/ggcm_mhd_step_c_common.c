
#include "ggcm_mhd_step_private.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_defs_extra.h"
#include "ggcm_mhd_diag_private.h"

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

#define OPT_TMP OPT_TMP_COMPAT
#define OPT_GGCM_CRDS OPT_GGCM_CRDS_LEGACY

// mhd options

#define OPT_EQN OPT_EQN_MHD_SCONS

#include "pde/pde_mhd_compat.c"
#include "pde/pde_mhd_get_dt.c"
#include "pde/pde_mhd_push.c"

// TODO:
// - handle remaining resistivity models
// - handle limit2, limit3
// - handle lowmask

// ======================================================================
// ggcm_mhd_step subclass "c"
//
// this class will do full predictor / corrector steps,
// ie., including primvar() etc.

struct ggcm_mhd_step_c {
  struct mhd_options opt;
};

#define ggcm_mhd_step_c(step) mrc_to_subobj(step, struct ggcm_mhd_step_c)

// ----------------------------------------------------------------------
// ggcm_mhd_step_c_setup_flds

static void
ggcm_mhd_step_c_setup_flds(struct ggcm_mhd_step *step)
{
  struct ggcm_mhd_step_c *sub = ggcm_mhd_step_c(step);
  struct ggcm_mhd *mhd = step->mhd;

  pde_mhd_set_options(mhd, &sub->opt);
  mrc_fld_set_type(mhd->fld, FLD_TYPE);
  mrc_fld_set_param_int(mhd->fld, "nr_ghosts", 2);
#if OPT_STAGGER == OPT_STAGGER_GGCM
  mrc_fld_dict_add_int(mhd->fld, "mhd_type", MT_SCONS_FC_GGCM);
#else
  mrc_fld_dict_add_int(mhd->fld, "mhd_type", MT_SCONS_FC);
#endif
  mrc_fld_set_param_int(mhd->fld, "nr_comps", _NR_FLDS);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_c_setup

static void
ggcm_mhd_step_c_setup(struct ggcm_mhd_step *step)
{
  struct ggcm_mhd *mhd = step->mhd;

  pde_mhd_setup(mhd, 5);
  pde_mhd_compat_setup(mhd);

  mhd->ymask = mrc_fld_make_view(mhd->fld, _YMASK, _YMASK + 1);
  mrc_fld_set(mhd->ymask, 1.);

  ggcm_mhd_step_setup_member_objs_sub(step);
  ggcm_mhd_step_setup_super(step);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_c_destroy

static void
ggcm_mhd_step_c_destroy(struct ggcm_mhd_step *step)
{
  mrc_fld_destroy(step->mhd->ymask);

  pde_free();
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_c_get_dt

static double
ggcm_mhd_step_c_get_dt(struct ggcm_mhd_step *step, struct mrc_fld *x)
{
  return pde_mhd_get_dt_scons(step->mhd, x, step->mhd->ymask);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_c_run

static void
ggcm_mhd_step_c_run(struct ggcm_mhd_step *step, struct mrc_fld *x)
{
  struct ggcm_mhd *mhd = step->mhd;

  assert(x == mhd->fld);
  int mhd_type;
  mrc_fld_get_param_int(x, "mhd_type", &mhd_type);

  struct mrc_fld *x_n = mrc_fld_make_view(x, _RR1, _RR1 + 8);
  mrc_fld_dict_add_int(x_n, "mhd_type", mhd_type);

  struct mrc_fld *x_star = mrc_fld_make_view(x, _RR2, _RR2 + 8);
  mrc_fld_dict_add_int(x_star, "mhd_type", mhd_type);

  // FIXME? It's not going to make a difference, but this is the
  // time at the beginning of the whole step, rather than the time of the current state
  s_mhd_time = mhd->time_code * mhd->tnorm; 

  ggcm_mhd_fill_ghosts(mhd, x_n, mhd->time_code);
  pde_mhd_pushstage(x, .5f * mhd->dt_code, 0);

  ggcm_mhd_fill_ghosts(mhd, x_star, mhd->time_code + .5 * mhd->dt_code);
  pde_mhd_pushstage(x, mhd->dt_code, 1);

  mrc_fld_destroy(x_n);
  mrc_fld_destroy(x_star);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_c_get_e_ec

static void
ggcm_mhd_step_c_get_e_ec(struct ggcm_mhd_step *step, struct mrc_fld *Eout,
                         struct mrc_fld *state_vec)
{
  // the state vector should already be FLD_TYPE, but Eout is the data type
  // of the output
  struct mrc_fld *E = mrc_fld_get_as(Eout, FLD_TYPE);
  struct mrc_fld *x = mrc_fld_get_as(state_vec, FLD_TYPE);

#if OPT_STAGGER == OPT_STAGGER_GGCM
  fld3d_foreach(i, j, k, 1, 0) {
#else
  fld3d_foreach(i, j, k, 0, 1) {
#endif
    F3(E, 0, i,j,k) = F3(x, _FLX, i,j,k);
    F3(E, 1, i,j,k) = F3(x, _FLY, i,j,k);
    F3(E, 2, i,j,k) = F3(x, _FLZ, i,j,k);
  } fld3d_foreach_end;

  mrc_fld_put_as(E, Eout);
  // FIXME, should use _put_as, but don't want copy-back
  if (strcmp(mrc_fld_type(state_vec), FLD_TYPE) != 0) {
    mrc_fld_destroy(x);
  }
} 

// ----------------------------------------------------------------------
// ggcm_mhd_step_c_diag_item_zmask_run

static void
ggcm_mhd_step_c_diag_item_zmask_run(struct ggcm_mhd_step *step,
				    struct ggcm_mhd_diag_item *item,
				    struct mrc_io *io, struct mrc_fld *f,
				    int diag_type, float plane)
{
  ggcm_mhd_diag_c_write_one_field(io, f, _ZMASK, "zmask", 1., diag_type, plane);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_c_diag_item_rmask_run

static void
ggcm_mhd_step_c_diag_item_rmask_run(struct ggcm_mhd_step *step,
				    struct ggcm_mhd_diag_item *item,
				    struct mrc_io *io, struct mrc_fld *f,
				    int diag_type, float plane)
{
  ggcm_mhd_diag_c_write_one_field(io, f, _RMASK, "rmask", 1., diag_type, plane);
}

// ----------------------------------------------------------------------
// subclass description

#define VAR(x) (void *)offsetof(struct ggcm_mhd_step_c, x)
static struct param ggcm_mhd_step_c_descr[] = {
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
// ggcm_mhd_step subclass "c_*"

struct ggcm_mhd_step_ops ggcm_mhd_step_c_ops = {
  .name                = ggcm_mhd_step_c_name,
  .size                = sizeof(struct ggcm_mhd_step_c),
  .param_descr         = ggcm_mhd_step_c_descr,
  .setup               = ggcm_mhd_step_c_setup,
  .destroy             = ggcm_mhd_step_c_destroy,
  .setup_flds          = ggcm_mhd_step_c_setup_flds,
  .get_dt              = ggcm_mhd_step_c_get_dt,
  .run                 = ggcm_mhd_step_c_run,
  .get_e_ec            = ggcm_mhd_step_c_get_e_ec,
  .diag_item_zmask_run = ggcm_mhd_step_c_diag_item_zmask_run,
  .diag_item_rmask_run = ggcm_mhd_step_c_diag_item_rmask_run,
};
