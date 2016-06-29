
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

// mhd options

#define OPT_EQN OPT_EQN_MHD_SCONS

#include "pde/pde_mhd_compat.c"
#include "pde/pde_mhd_get_dt.c"
#include "pde/pde_mhd_rmaskn.c"
#include "pde/pde_mhd_pushfluid.c"
#include "pde/pde_mhd_pushfield.c"
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

  struct mrc_fld *f_Uhalf;
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
  mrc_fld_dict_add_int(mhd->fld, "mhd_type", MT_SEMI_CONSERVATIVE);
  mrc_fld_set_param_int(mhd->fld, "nr_comps", _NR_FLDS);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_c2_setup

static void
ggcm_mhd_step_c2_setup(struct ggcm_mhd_step *step)
{
  struct ggcm_mhd_step_c2 *sub = ggcm_mhd_step_c2(step);
  struct ggcm_mhd *mhd = step->mhd;

  pde_setup(mhd->fld);
  pde_mhd_setup(mhd);
  pde_mhd_compat_setup(mhd);

  mhd->ymask = mrc_fld_make_view(mhd->fld, _YMASK, _YMASK + 1);
  mrc_fld_set(mhd->ymask, 1.);

  sub->f_Uhalf = ggcm_mhd_get_3d_fld(mhd, 8);
  mrc_fld_dict_add_int(sub->f_Uhalf, "mhd_type", MT_SEMI_CONSERVATIVE);

  ggcm_mhd_step_setup_member_objs_sub(step);
  ggcm_mhd_step_setup_super(step);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_c2_destroy

static void
ggcm_mhd_step_c2_destroy(struct ggcm_mhd_step *step)
{
  mrc_fld_destroy(step->mhd->ymask);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_c2_get_dt

static double
ggcm_mhd_step_c2_get_dt(struct ggcm_mhd_step *step, struct mrc_fld *x)
{
  return pde_mhd_get_dt_scons(step->mhd, x);
}

// ----------------------------------------------------------------------
// patch_push

static void
patch_push(fld3d_t p_Unext, mrc_fld_data_t dt, fld3d_t p_Ucurr,
	   fld3d_t p_W, fld3d_t p_cmsv,
	   fld3d_t p_ymask, fld3d_t p_zmask,
	   int stage)
{
  fld3d_t p_rmask = fld3d_make_tmp(1, _RMASK), p_resis = fld3d_make_tmp(1, _RESIS);
  fld3d_t p_Jcc = fld3d_make_tmp(3, _CURRX);

  patch_rmaskn(p_rmask, p_zmask);
  patch_pushfluid(p_Unext, dt, p_Unext, p_Ucurr, p_W,
		  p_cmsv, p_ymask, p_zmask, stage);
  patch_pushfield(p_Unext, dt, p_Unext, p_Ucurr, p_W,
		  p_zmask, p_rmask, p_resis, p_Jcc, stage);
}

// ----------------------------------------------------------------------
// patch_pushstage

static void
patch_pushstage(fld3d_t p_Unext, mrc_fld_data_t dt, fld3d_t p_Ucurr, fld3d_t p_ymask,
		fld3d_t p_zmask, int stage)
{
  fld3d_t p_W     = fld3d_make_tmp(5, _RR);
  fld3d_t p_cmsv  = fld3d_make_tmp(1, _CMSV);

  patch_primvar(p_W, p_Ucurr, p_cmsv);
  patch_badval_checks_sc(p_Ucurr, p_W); // FIXME, incorporate

  if (stage == 0) {
    fld3d_t p_bcc = fld3d_make_tmp(3, _BX);
    patch_primbb(p_bcc, p_Ucurr);
    patch_zmaskn(p_zmask, p_W, p_bcc, p_ymask);
  }

  patch_push(p_Unext, dt, p_Ucurr, p_W, p_cmsv,
	     p_ymask, p_zmask, stage);
}

// ----------------------------------------------------------------------
// pushstage

static void
pushstage(struct mrc_fld *f_Unext, mrc_fld_data_t dt, struct mrc_fld *f_Ucurr,
	  struct mrc_fld *f_ymask, struct mrc_fld *f_zmask, int stage)
{
  fld3d_t p_Unext, p_Ucurr, p_ymask, p_zmask;
  fld3d_setup(&p_Unext, f_Unext);
  fld3d_setup(&p_Ucurr, f_Ucurr);
  fld3d_setup(&p_ymask, f_ymask);
  fld3d_setup(&p_zmask, f_zmask);

  pde_for_each_patch(p) {
    fld3d_t *patches[] = { &p_Unext, &p_Ucurr, &p_ymask, &p_zmask, NULL };
    fld3d_get_list(p, patches);
    patch_pushstage(p_Unext, dt, p_Ucurr, p_ymask, p_zmask, stage);
    fld3d_put_list(p, patches);
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_c2_run

static void
ggcm_mhd_step_c2_run(struct ggcm_mhd_step *step, struct mrc_fld *x)
{
  struct ggcm_mhd_step_c2 *sub = ggcm_mhd_step_c2(step);
  struct ggcm_mhd *mhd = step->mhd;

  struct mrc_fld *f_ymask = mhd->ymask;
  struct mrc_fld *f_zmask = mrc_fld_make_view(x, _ZMASK, _ZMASK + 1);
  struct mrc_fld *f_U = mrc_fld_make_view(x, _RR1, _RR1 + 8);
  mrc_fld_dict_add_int(f_U, "mhd_type", MT_SEMI_CONSERVATIVE);
  struct mrc_fld *f_Uhalf = sub->f_Uhalf;

  // FIXME? It's not going to make a difference, but this is the
  // time at the beginning of the whole step, rather than the time of the current state
  s_mhd_time = mhd->time; 

  // set f_Uhalf = f_U
  mrc_fld_copy(f_Uhalf, f_U);
  // then advance f_Uhalf += .5f * dt * rhs(f_U)
  ggcm_mhd_fill_ghosts(mhd, f_U, 0, mhd->time);
  pushstage(f_Uhalf, .5f * mhd->dt, f_U, f_ymask, f_zmask, 0);

  // f_U += dt * rhs(f_Uhalf)
  ggcm_mhd_fill_ghosts(mhd, f_Uhalf, 0, mhd->time + mhd->bndt);
  pushstage(f_U, mhd->dt, f_Uhalf, f_ymask, f_zmask, 1);

  mrc_fld_destroy(f_U);
  mrc_fld_destroy(f_zmask);
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

  fld3d_foreach(i, j, k, 0, 1) {
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
  .setup               = ggcm_mhd_step_c2_setup,
  .destroy             = ggcm_mhd_step_c2_destroy,
  .setup_flds          = ggcm_mhd_step_c2_setup_flds,
  .get_dt              = ggcm_mhd_step_c2_get_dt,
  .run                 = ggcm_mhd_step_c2_run,
  .get_e_ec            = ggcm_mhd_step_c2_get_e_ec,
  .diag_item_zmask_run = ggcm_mhd_step_c2_diag_item_zmask_run,
  .diag_item_rmask_run = ggcm_mhd_step_c2_diag_item_rmask_run,
};
