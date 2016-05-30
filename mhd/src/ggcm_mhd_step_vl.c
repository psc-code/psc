
#include <ggcm_mhd_step_private.h>

#include <ggcm_mhd_private.h>
#include <ggcm_mhd_defs.h>
#include <ggcm_mhd_diag.h>

#include <mrc_domain.h>
#include <mrc_ddc.h>
#include <mrc_ts.h>
#include <mrc_io.h>

#include <mrc_fld_as_double.h>
#define F1(f, m, i) MRC_D2(f, m, i)

#include "pde/pde_setup.c"
#include "pde/pde_mhd_setup.c"
#include "pde/pde_mhd_line.c"
#include "pde/pde_mhd_convert.c"
#include "pde/pde_mhd_reconstruct.c"
#include "pde/pde_mhd_riemann.c"

#include "mhd_3d.c"

static int ldims[3];

// ======================================================================
// ggcm_mhd_step subclass "vl"

struct ggcm_mhd_step_vl {
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
flux_pred(struct ggcm_mhd_step *step, struct mrc_fld *flux[3], struct mrc_fld *x, struct mrc_fld *B_cc,
	  int ldim, int nghost, int j, int k, int dir, int p)
{
  struct ggcm_mhd_step_vl *sub = ggcm_mhd_step_vl(step);
  fld1d_state_t U = sub->U, U_l = sub->U_l, U_r = sub->U_r;
  fld1d_state_t W = sub->W, W_l = sub->W_l, W_r = sub->W_r;
  fld1d_state_t F = sub->F;

  pick_line_sc(U, x, ldim, 2, 2, j, k, dir, p);
  mhd_prim_from_sc(W, U, ldim, 2, 2); // for up to plm reconstruction
  mhd_reconstruct_pcm_run_sc(U_l, U_r, W_l, W_r, W, NULL,
			     ldim, 1, 1, dir);
  mhd_riemann_rusanov_run_hydro(F, U_l, U_r, W_l, W_r, ldim, 0, 1, dir);
  put_line_sc(flux[dir], F, ldim, 0, 1, j, k, dir, p);
}

static void
flux_corr(struct ggcm_mhd_step *step, struct mrc_fld *flux[3], struct mrc_fld *x, struct mrc_fld *B_cc,
	  int ldim, int nghost, int j, int k, int dir, int p)
{
  struct ggcm_mhd_step_vl *sub = ggcm_mhd_step_vl(step);
  fld1d_state_t U = sub->U, U_l = sub->U_l, U_r = sub->U_r;
  fld1d_state_t W = sub->W, W_l = sub->W_l, W_r = sub->W_r;
  fld1d_state_t F = sub->F;

  pick_line_sc(U, x, ldim, 2, 2, j, k, dir, p);
  mhd_prim_from_sc(W, U, ldim, 2, 2); // for up to plm reconstruction
  mhd_reconstruct_plm_run_sc(U_l, U_r, W_l, W_r, W, NULL,
			     ldim, 1, 1, dir);
  mhd_riemann_rusanov_run_hydro(F, U_l, U_r, W_l, W_r, ldim, 0, 1, dir);
  put_line_sc(flux[dir], F, ldim, 0, 1, j, k, dir, p);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_vl_setup

static void
ggcm_mhd_step_vl_setup(struct ggcm_mhd_step *step)
{
  struct ggcm_mhd_step_vl *sub = ggcm_mhd_step_vl(step);
  struct ggcm_mhd *mhd = step->mhd;

  pde_setup(mhd->fld);
  pde_mhd_setup(mhd);

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
// newstep_hydro

static mrc_fld_data_t
newstep_hydro(struct ggcm_mhd *mhd, struct mrc_fld *x)
{
  int gdims[3];
  mrc_domain_get_global_dims(x->_domain, gdims);

  mrc_fld_data_t gamma = mhd->par.gamm;
  mrc_fld_data_t gamma_minus_1 = gamma - 1.;
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);

  mrc_fld_data_t min_dt = 1e10;
  for (int p = 0; p < mrc_fld_nr_patches(x); p++) {
    double dx[3]; mrc_crds_get_dx(crds, p, dx);

    mrc_fld_foreach(x, i,j,k, 0, 0) {
      mrc_fld_data_t rri = 1.f / RR(x, i,j,k);
      mrc_fld_data_t vx = RVX(x, i,j,k) * rri;
      mrc_fld_data_t vy = RVY(x, i,j,k) * rri;
      mrc_fld_data_t vz = RVZ(x, i,j,k) * rri;
      
      /* compute sound speed */
      mrc_fld_data_t pp = mrc_fld_max(gamma_minus_1*(EE(x, i,j,k) - .5f*RR(x, i,j,k) * (sqr(vx) + sqr(vy) + sqr(vz))), 1e-15);
      mrc_fld_data_t cs = mrc_fld_sqrt(gamma * pp * rri);
      
      /* compute min dt based on maximum wave velocity */
      if (gdims[0] > 1) { min_dt = mrc_fld_min(min_dt, dx[0] / (mrc_fld_abs(vx) + cs)); }
      if (gdims[1] > 1) { min_dt = mrc_fld_min(min_dt, dx[1] / (mrc_fld_abs(vy) + cs)); }
      if (gdims[2] > 1) { min_dt = mrc_fld_min(min_dt, dx[2] / (mrc_fld_abs(vz) + cs)); }
    } mrc_fld_foreach_end;
  }

  mrc_fld_data_t cfl = mhd->par.thx;
  mrc_fld_data_t local_dt = cfl * min_dt;
  mrc_fld_data_t global_dt;
  MPI_Allreduce(&local_dt, &global_dt, 1, MPI_MRC_FLD_DATA_T, MPI_MIN, ggcm_mhd_comm(mhd));
  return global_dt;
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
  struct mrc_fld *flux[3] = { ggcm_mhd_get_3d_fld(mhd, 5),
			      ggcm_mhd_get_3d_fld(mhd, 5),
			      ggcm_mhd_get_3d_fld(mhd, 5), };

  // CFL CONDITION

  ggcm_mhd_fill_ghosts(mhd, x, 0, mhd->time);

  mhd->dt = newstep_hydro(mhd, x);
  mrc_fld_data_t dt = mhd->dt;

  // PREDICTOR

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
    ggcm_mhd_fill_ghosts(mhd, x, 0, mhd->time);
    ggcm_mhd_diag_run_now(diag, x, DIAG_TYPE_3D, cnt++);
  }

  // ghosts have already been set
  mrc_fld_copy_range(x_half, x, 0, 5);
  mhd_fluxes(step, flux, x, NULL, 0, 99, flux_pred);
  update_finite_volume_uniform(mhd, x_half, flux, NULL, .5 * dt, 0, 0, false);

  // CORRECTOR

  ggcm_mhd_fill_ghosts(mhd, x_half, 0, mhd->time);
  mhd_fluxes(step, flux, x_half, NULL, 0, 99, flux_corr);
  update_finite_volume_uniform(mhd, x, flux, NULL, dt, 0, 0, true);

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
  struct ggcm_mhd *mhd = step->mhd;

  mrc_fld_set_type(mhd->fld, FLD_TYPE);
  mrc_fld_set_param_int(mhd->fld, "nr_ghosts", 2);
  mrc_fld_dict_add_int(mhd->fld, "mhd_type", MT_FULLY_CONSERVATIVE);
  mrc_fld_set_param_int(mhd->fld, "nr_comps", 8); // FIXME, should be 5, but needs testing
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_vl subclass description

#define VAR(x) (void *)offsetof(struct ggcm_mhd_step_vl, x)
static struct param ggcm_mhd_step_vl_descr[] = {
  { "debug_dump"      , VAR(debug_dump)      , PARAM_BOOL(false)            },
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
  .run              = ggcm_mhd_step_vl_run,
  .setup_flds       = ggcm_mhd_step_vl_setup_flds,
};

