
#include <ggcm_mhd_step_private.h>

#include <ggcm_mhd_private.h>
#include <ggcm_mhd_defs.h>
#include <ggcm_mhd_diag.h>

#include <mrc_domain.h>
#include <mrc_ddc.h>
#include <mrc_ts.h>
#include <mrc_io.h>

#include <mrc_fld_as_double.h>

#include "pde/pde_defs.h"

#define OPT_EQN OPT_EQN_HD

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

static int ldims[3];

// ======================================================================
// ggcm_mhd_step subclass "vl"

struct ggcm_mhd_step_vl {
  struct mhd_options opt;

  bool debug_dump;

  fld1d_state_t U;
  fld1d_state_t U_l;
  fld1d_state_t U_r;
  fld1d_state_t W;
  fld1d_state_t W_l;
  fld1d_state_t W_r;
  fld1d_state_t F;
};

#define ggcm_mhd_step_vl(step) mrc_to_subobj(step, struct ggcm_mhd_step_vl)

// ======================================================================

static void
fluxes_pred(struct ggcm_mhd_step *step, struct mrc_fld *flux[3], struct mrc_fld *x)
{
  struct ggcm_mhd_step_vl *sub = ggcm_mhd_step_vl(step);
  fld1d_state_t U = sub->U, U_l = sub->U_l, U_r = sub->U_r;
  fld1d_state_t W = sub->W, W_l = sub->W_l, W_r = sub->W_r;
  fld1d_state_t F = sub->F;

  for (int p = 0; p < mrc_fld_nr_patches(x); p++) {
    pde_for_each_dir(dir) {
      int ldim = s_ldims[dir];
      pde_for_each_line(dir, j, k, 0) {
	mhd_get_line_state(U, x, j, k, dir, p, -2, ldim + 2);
	mhd_prim_from_cons(W, U, -2, ldim + 2); // for up to plm reconstruction
	mhd_reconstruct_pcm(U_l, U_r, W_l, W_r, W, (fld1d_t) {}, 0, ldim + 1);
	mhd_riemann(F, U_l, U_r, W_l, W_r, 0, ldim + 1);
	mhd_put_line_state(flux[dir], F, j, k, dir, p, 0, ldim + 1);
      }
    }
  }
}

static void
fluxes_corr(struct ggcm_mhd_step *step, struct mrc_fld *flux[3], struct mrc_fld *x)
{
  struct ggcm_mhd_step_vl *sub = ggcm_mhd_step_vl(step);
  fld1d_state_t U = sub->U, U_l = sub->U_l, U_r = sub->U_r;
  fld1d_state_t W = sub->W, W_l = sub->W_l, W_r = sub->W_r;
  fld1d_state_t F = sub->F;

  for (int p = 0; p < mrc_fld_nr_patches(x); p++) {
    pde_for_each_dir(dir) {
      int ldim = s_ldims[dir];
      pde_for_each_line(dir, j, k, 0) {
	mhd_get_line_state(U, x, j, k, dir, p, -2, ldim + 2);
	mhd_prim_from_cons(W, U, -2, ldim + 2); // for up to plm reconstruction
	mhd_reconstruct(U_l, U_r, W_l, W_r, W, (fld1d_t) {}, 0, ldim + 1);
	mhd_riemann(F, U_l, U_r, W_l, W_r, 0, ldim + 1);
	mhd_put_line_state(flux[dir], F, j, k, dir, p, 0, ldim + 1);
      }
    }
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_vl_setup

static void
ggcm_mhd_step_vl_setup(struct ggcm_mhd_step *step)
{
  struct ggcm_mhd_step_vl *sub = ggcm_mhd_step_vl(step);
  struct ggcm_mhd *mhd = step->mhd;

  pde_mhd_setup(mhd, mrc_fld_nr_comps(mhd->fld));

  fld1d_state_setup(&sub->U);
  fld1d_state_setup(&sub->U_l);
  fld1d_state_setup(&sub->U_r);
  fld1d_state_setup(&sub->W);
  fld1d_state_setup(&sub->W_l);
  fld1d_state_setup(&sub->W_r);
  fld1d_state_setup(&sub->F);

  ggcm_mhd_step_setup_member_objs_sub(step);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_vl_destroy

static void
ggcm_mhd_step_vl_destroy(struct ggcm_mhd_step *step)
{
  pde_free();
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_vl_get_dt

static double
ggcm_mhd_step_vl_get_dt(struct ggcm_mhd_step *step, struct mrc_fld *x)
{
  return pde_mhd_get_dt(step->mhd, x);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_vl_run

static void
ggcm_mhd_step_vl_run(struct ggcm_mhd_step *step, struct mrc_fld *x)
{
  struct ggcm_mhd_step_vl *sub = ggcm_mhd_step_vl(step);
  struct ggcm_mhd *mhd = step->mhd;

  ldims[0] = mrc_fld_spatial_dims(x)[0];
  ldims[1] = mrc_fld_spatial_dims(x)[1];
  ldims[2] = mrc_fld_spatial_dims(x)[2];

  struct mrc_fld *x_half = ggcm_mhd_get_3d_fld(mhd, 8);
  mrc_fld_dict_add_int(x_half, "mhd_type", MT_FCONS_CC); // FIXME
  struct mrc_fld *flux[3] = { ggcm_mhd_get_3d_fld(mhd, 5),
			      ggcm_mhd_get_3d_fld(mhd, 5),
			      ggcm_mhd_get_3d_fld(mhd, 5), };

  mrc_fld_data_t dt = mhd->dt_code;

  // PREDICTOR

  ggcm_mhd_fill_ghosts(mhd, x, mhd->time_code);

  if (sub->debug_dump) {
    static struct ggcm_mhd_diag *diag;
    static int cnt;
    if (!diag) {
      diag = ggcm_mhd_diag_create(ggcm_mhd_comm(mhd));
      ggcm_mhd_diag_set_type(diag, "c");
      ggcm_mhd_diag_set_param_obj(diag, "mhd", mhd);
      ggcm_mhd_diag_set_param_string(diag, "run", "dbg");
      ggcm_mhd_diag_set_param_string(diag, "fields", "rr1:rv1:uu1:b1:rr:v:pp:b:divb");
      ggcm_mhd_diag_setup(diag);
      ggcm_mhd_diag_view(diag);
    }
    ggcm_mhd_fill_ghosts(mhd, x, mhd->time_code);
    ggcm_mhd_diag_run_now(diag, x, DIAG_TYPE_3D, cnt++);
  }

  fld3d_t _x_half, _x, _flux[3];
  fld3d_setup(&_x_half, x_half);
  fld3d_setup(&_x, x);
  for (int d = 0; d < 3; d++) {
    fld3d_setup(&_flux[d], flux[d]);
  }

  // ghosts have already been set
  mrc_fld_copy_range(x_half, x, 0, 5);
  fluxes_pred(step, flux, x);
  for (int p = 0; p < mrc_fld_nr_patches(x); p++) {
    pde_patch_set(p);
    fld3d_get(&_x_half, p);
    for (int d = 0; d < 3; d++) {
      fld3d_get(&_flux[d], p);
    }

    mhd_update_finite_volume(mhd, _x_half, _flux, (fld3d_t) {}, .5f * dt, 0, 0);

    fld3d_put(&_x_half, p);
    for (int d = 0; d < 3; d++) {
      fld3d_put(&_flux[d], p);
    }
  }

  // CORRECTOR

  ggcm_mhd_fill_ghosts(mhd, x_half, mhd->time_code);
  fluxes_corr(step, flux, x);
  ggcm_mhd_correct_fluxes(mhd, flux);
  for (int p = 0; p < mrc_fld_nr_patches(x); p++) {
    pde_patch_set(p);
    fld3d_get(&_x, p);
    for (int d = 0; d < 3; d++) {
      fld3d_get(&_flux[d], p);
    }

    mhd_update_finite_volume(mhd, _x, _flux, (fld3d_t) {}, dt, 0, 0);

    fld3d_put(&_x, p);
    for (int d = 0; d < 3; d++) {
      fld3d_put(&_flux[d], p);
    }
  }

  // clean up

  ggcm_mhd_put_3d_fld(mhd, x_half);
  ggcm_mhd_put_3d_fld(mhd, flux[0]);
  ggcm_mhd_put_3d_fld(mhd, flux[1]);
  ggcm_mhd_put_3d_fld(mhd, flux[2]);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_vl_setup_flds

static void
ggcm_mhd_step_vl_setup_flds(struct ggcm_mhd_step *step)
{
  struct ggcm_mhd_step_vl *sub = ggcm_mhd_step_vl(step);
  struct ggcm_mhd *mhd = step->mhd;

  pde_mhd_set_options(mhd, &sub->opt);

  mrc_fld_set_type(mhd->fld, FLD_TYPE);
  mrc_fld_set_param_int(mhd->fld, "nr_ghosts", 2);
  mrc_fld_dict_add_int(mhd->fld, "mhd_type", MT_FCONS_CC);  // FIXME
  mrc_fld_set_param_int(mhd->fld, "nr_comps", 8); // FIXME, should be 5, but needs testing
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_vl subclass description

#define VAR(x) (void *)offsetof(struct ggcm_mhd_step_vl, x)
static struct param ggcm_mhd_step_vl_descr[] = {
  { "eqn"                , VAR(opt.eqn)            , PARAM_SELECT(OPT_EQN,
								  opt_eqn_descr)                },
  { "limiter"            , VAR(opt.limiter)        , PARAM_SELECT(OPT_LIMITER_FLAT,
								  opt_limiter_descr)            },
  { "riemann"            , VAR(opt.riemann)        , PARAM_SELECT(OPT_RIEMANN_RUSANOV,
								  opt_riemann_descr)            },

  { "debug_dump"         , VAR(debug_dump)         , PARAM_BOOL(false)                          },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_step_vl_ops
//
// This scheme is a hydro-solver only -- still working on 8-component fields,
// but ignoring the magnetic field completely.

struct ggcm_mhd_step_ops ggcm_mhd_step_vl_ops = {
  .name             = "vl",
  .size             = sizeof(struct ggcm_mhd_step_vl),
  .param_descr      = ggcm_mhd_step_vl_descr,
  .setup            = ggcm_mhd_step_vl_setup,
  .destroy          = ggcm_mhd_step_vl_destroy,
  .run              = ggcm_mhd_step_vl_run,
  .setup_flds       = ggcm_mhd_step_vl_setup_flds,
  .get_dt           = ggcm_mhd_step_vl_get_dt,
};

