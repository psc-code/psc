
#include "ggcm_mhd_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_defs_extra.h"
#include "ggcm_mhd_crds.h"
#include "ggcm_mhd_crds_private.h"
#include "ggcm_mhd_step.h"
#include "ggcm_mhd_diag.h"
#include "ggcm_mhd_bnd.h"
#include "ggcm_mhd_ic.h"

#include <mrc_domain.h>
#include <mrc_ddc.h>
#include <mrc_ts.h>
#include <mrc_ts_monitor.h>
#include <mrc_io.h>
#include <mrc_profile.h>
#include <mrc_physics.h>

#include <assert.h>
#include <string.h>
#include <math.h>

#define ggcm_mhd_ops(mhd) ((struct ggcm_mhd_ops *) mhd->obj.ops)

static const char *fldname[_NR_FLDS] = {
  [_RR1 ] = "rr1",
  [_RV1X] = "rv1x",
  [_RV1Y] = "rv1y",
  [_RV1Z] = "rv1z",
  [_UU1 ] = "uu1",
  [_B1X ] = "b1x",
  [_B1Y ] = "b1y",
  [_B1Z ] = "b1z",

  [_RR2 ] = "rr2",
  [_RV2X] = "rv2x",
  [_RV2Y] = "rv2y",
  [_RV2Z] = "rv2z",
  [_UU2 ] = "uu2",
  [_B2X ] = "b2x",
  [_B2Y ] = "b2y",
  [_B2Z ] = "b2z",

  [_YMASK] = "ymask",
  [_ZMASK] = "zmask",
  [_CMSV ] = "cmsv",

  [_RR  ] = "rr",
  [_PP  ] = "pp",
  [_VX  ] = "vx",
  [_VY  ] = "vy",
  [_VZ  ] = "vz",
  [_BX  ] = "bx",
  [_BY  ] = "by",
  [_BZ  ] = "bz",

  [_TMP1] = "tmp1",
  [_TMP2] = "tmp2",
  [_TMP3] = "tmp3",
  [_TMP4] = "tmp4",

  [_FLX ] = "ex",
  [_FLY ] = "ey",
  [_FLZ ] = "ez",

  [_CX  ] = "cx",
  [_CY  ] = "cy",
  [_CZ  ] = "cz",

  [_XTRA1] = "xtra1",
  [_XTRA2] = "xtra2",

  [_RESIS] = "resis",

  [_CURRX] = "currx",
  [_CURRY] = "curry",
  [_CURRZ] = "currz",

  [_RMASK] = "rmask",
};

// ----------------------------------------------------------------------
// ggcm_mhd methods

static void
_ggcm_mhd_create(struct ggcm_mhd *mhd)
{
  mrc_domain_set_type(mhd->domain, "simple");
  ggcm_mhd_bnd_set_name(mhd->bnd1, "bnd1");
  
  ggcm_mhd_crds_set_param_obj(mhd->crds, "domain", mhd->domain);
  ggcm_mhd_step_set_param_obj(mhd->step, "mhd", mhd);
  ggcm_mhd_diag_set_param_obj(mhd->diag, "mhd", mhd);
  ggcm_mhd_bnd_set_param_obj(mhd->bnd, "mhd", mhd);
  ggcm_mhd_bnd_set_param_obj(mhd->bnd1, "mhd", mhd);
  ggcm_mhd_ic_set_param_obj(mhd->ic, "mhd", mhd);

  mrc_fld_set_name(mhd->fld, "ggcm_mhd_fld");
  mrc_fld_set_param_obj(mhd->fld, "domain", mhd->domain);
  mrc_fld_set_param_int(mhd->fld, "nr_spatial_dims", 3);
  mrc_fld_dict_add_obj(mhd->fld, "mhd", mhd);
}

// ----------------------------------------------------------------------
// ggcm_mhd_set_state
//
// update Fortran common blocks from C ggcm_mhd state

void
ggcm_mhd_set_state(struct ggcm_mhd *mhd)
{
  struct ggcm_mhd_ops *ops = ggcm_mhd_ops(mhd);
  if (ops->set_state) {
    ops->set_state(mhd);
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_pre_step

void
ggcm_mhd_pre_step(struct ggcm_mhd *mhd, struct mrc_ts *ts, struct mrc_fld *fld)
{
  struct ggcm_mhd_ops *ops = ggcm_mhd_ops(mhd);
  if (ops->pre_step) {
    ops->pre_step(mhd, ts, fld);
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_post_step

void
ggcm_mhd_post_step(struct ggcm_mhd *mhd, struct mrc_ts *ts, struct mrc_fld *fld)
{
  struct ggcm_mhd_ops *ops = ggcm_mhd_ops(mhd);
  if (ops->post_step) {
    ops->post_step(mhd, ts, fld);
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_setup_internal

static void
ggcm_mhd_setup_internal(struct ggcm_mhd *mhd)
{
  const int *ghost_dims = mrc_fld_ghost_dims(mhd->fld);
  const int *dims = mrc_fld_dims(mhd->fld);
  int shift = 0;
  for (int d = 0; d < 3; d++) {
    // local domain size
    mhd->im[d] = dims[d+shift];
    // local domain size incl ghost points
    mhd->img[d] = ghost_dims[d+shift];
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_read

static void
_ggcm_mhd_read(struct ggcm_mhd *mhd, struct mrc_io *io)
{
  ggcm_mhd_read_member_objs(mhd, io);

  ggcm_mhd_setup_internal(mhd);
}

// ----------------------------------------------------------------------
// ggcm_mhd_setup_normalization

static void
ggcm_mhd_setup_normalization(struct ggcm_mhd *mhd)
{
  double x0 = mhd->par.norm_length;
  double b0 = mhd->par.norm_B;
  double rr0 = mhd->par.norm_density;
  double mu00 = mhd->par.norm_mu0;
  
  mhd->xxnorm = x0;
  mhd->bbnorm = b0;
  mhd->rrnorm = rr0;
  mhd->vvnorm = b0 / sqrt(mu00 * rr0);
  mhd->ppnorm = sqr(mhd->vvnorm) * rr0;
  mhd->ccnorm = b0 / (x0 * mu00);
  mhd->eenorm = mhd->vvnorm * b0;
  mhd->resnorm = mhd->eenorm / mhd->ccnorm;
  mhd->tnorm = x0 / mhd->vvnorm;
  mhd->qqnorm = mhd->ccnorm / (mhd->rrnorm * mhd->vvnorm);

  MPI_Comm comm = ggcm_mhd_comm(mhd);
  mpi_printf(comm, "NORMALIZATION: based on x0 = %g m\n", x0);
  mpi_printf(comm, "NORMALIZATION:          B0 = %g T\n", b0);
  mpi_printf(comm, "NORMALIZATION:         rr0 = %g kg/m^3\n", rr0);
  mpi_printf(comm, "NORMALIZATION:        mu00 = %g N/A^2\n", mu00);
  mpi_printf(comm, "NORMALIZATION:    mu0_code = %g N/A^2\n", mhd->par.mu0_code);
  mpi_printf(comm, "NORMALIZATION: xxnorm  = %g m\n", mhd->xxnorm);
  mpi_printf(comm, "NORMALIZATION: bbnorm  = %g T\n", mhd->bbnorm);
  mpi_printf(comm, "NORMALIZATION: rrnorm  = %g 1/m^3\n", mhd->rrnorm);
  mpi_printf(comm, "NORMALIZATION: vvnorm  = %g m/s\n", mhd->vvnorm);
  mpi_printf(comm, "NORMALIZATION: ppnorm  = %g Pa\n", mhd->ppnorm);
  mpi_printf(comm, "NORMALIZATION: ccnorm  = %g A/m^2\n", mhd->ccnorm);
  mpi_printf(comm, "NORMALIZATION: eenorm  = %g V/m\n", mhd->eenorm);
  mpi_printf(comm, "NORMALIZATION: resnorm = %g \n", mhd->resnorm);
  mpi_printf(comm, "NORMALIZATION: tnorm   = %g s\n", mhd->tnorm);
  mpi_printf(comm, "NORMALIZATION: qqnorm  = %g C\n", mhd->qqnorm);

  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  mrc_crds_set_param_double(crds, "norm_length", mhd->xxnorm);
  mrc_crds_set_param_double(crds, "norm_length_scale", mhd->par.xxnorm0);

  mhd->xxnorm /= mhd->par.xxnorm0;
  mhd->bbnorm /= mhd->par.bbnorm0;
  mhd->rrnorm /= mhd->par.rrnorm0;
  mhd->vvnorm /= mhd->par.vvnorm0;
  mhd->ppnorm /= mhd->par.ppnorm0;
  mhd->ccnorm /= mhd->par.ccnorm0;
  mhd->eenorm /= mhd->par.eenorm0;
  mhd->resnorm /= mhd->par.resnorm0;
  mhd->tnorm /= mhd->par.tnorm0;
  mhd->qqnorm /= mhd->par.qqnorm0;

  mpi_printf(comm, "NORMALIZATION: the following are the internally used scale factors\n"
	     "to convert code units to external units, including prefix. So, e.g., \n"
	     "internal (normalized) B to external B in nT\n");
  mpi_printf(comm, "NORMALIZATION: int. xxnorm  = %g\n", mhd->xxnorm);
  mpi_printf(comm, "NORMALIZATION: int. bbnorm  = %g\n", mhd->bbnorm);
  mpi_printf(comm, "NORMALIZATION: int. rrnorm  = %g\n", mhd->rrnorm);
  mpi_printf(comm, "NORMALIZATION: int. vvnorm  = %g\n", mhd->vvnorm);
  mpi_printf(comm, "NORMALIZATION: int. ppnorm  = %g\n", mhd->ppnorm);
  mpi_printf(comm, "NORMALIZATION: int. ccnorm  = %g\n", mhd->ccnorm);
  mpi_printf(comm, "NORMALIZATION: int. eenorm  = %g\n", mhd->eenorm);
  mpi_printf(comm, "NORMALIZATION: int. resnorm = %g\n", mhd->resnorm);
  mpi_printf(comm, "NORMALIZATION: int. tnorm   = %g\n", mhd->tnorm);
  mpi_printf(comm, "NORMALIZATION: int. qqnorm  = %g\n", mhd->qqnorm);
}

// ----------------------------------------------------------------------
// ggcm_mhd_setup_gk_norm

static void
ggcm_mhd_setup_gk_norm(struct ggcm_mhd *mhd)
{
  if (mhd->par.gk_norm) {
    if (mhd->par.d_i == 0.) {
      // FIXME, seems it'd be better to not call this if we aren't using a gkeyll stepper
      mpi_printf(ggcm_mhd_comm(mhd), "WARNING: ggcm_mhd_setup_gk_norm(): d_i == 0\n");
      return;
    }
    double mu0 = 1.;
    double rr_e = mhd->par.gk_norm_rr / mhd->rrnorm / (1. + mhd->par.gk_norm_mi_over_me);
    double rr_i = mhd->par.gk_norm_rr / mhd->rrnorm - rr_e;
    assert(mhd->par.d_i > 0.);
    double ion_mass = 1.;
    double electron_mass = ion_mass / mhd->par.gk_norm_mi_over_me;
    double e = ion_mass / mhd->par.d_i / sqrt(mu0 * rr_i);
    double pressure_ratio = mhd->par.gk_norm_ppi_over_ppe;
    
    ggcm_mhd_set_param_double(mhd, "gk_speed_of_light", mhd->par.gk_norm_speed_of_light / mhd->vvnorm);
    ggcm_mhd_set_param_float_array(mhd, "gk_charge", 2, (float[2]) { -e, e });
    ggcm_mhd_set_param_float_array(mhd, "gk_mass", 2, (float[2]) { electron_mass, ion_mass });
    ggcm_mhd_set_param_float_array(mhd, "gk_pressure_ratios", 2, 
				   (float[2]) { 1./(1.+pressure_ratio), 1./(1.+1./pressure_ratio) });
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_setup

static void
_ggcm_mhd_setup(struct ggcm_mhd *mhd)
{
  ggcm_mhd_setup_normalization(mhd);
  ggcm_mhd_setup_gk_norm(mhd);

  ggcm_mhd_step_setup_flds(mhd->step);
  for (int m = 0; m < mrc_fld_nr_comps(mhd->fld); m++) {
    mrc_fld_set_comp_name(mhd->fld, m, fldname[m]);
  }

  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  // only set the sw on the domain's crds if they're not already set
  if (1||crds->sw == 0) {
    mrc_crds_set_param_int(crds, "sw", mhd->fld->_nr_ghosts);
  }

  if (mhd->amr > 0) {
    ggcm_mhd_setup_amr_domain(mhd);
  }

  ggcm_mhd_setup_member_objs(mhd);
  ggcm_mhd_setup_internal(mhd);

  if (mhd->amr > 0) {
    // FIXME, all leaked
    mhd->ddc_amr_cc = ggcm_mhd_create_amr_ddc(mhd);
    mhd->ddc_amr_E  = ggcm_mhd_create_amr_ddc_E(mhd);
    for (int d = 0; d < 3; d++) {
      mhd->ddc_amr_flux[d] = ggcm_mhd_create_amr_ddc_flux(mhd, d);
    }
  }
}

void
ggcm_mhd_fill_ghosts(struct ggcm_mhd *mhd, struct mrc_fld *fld, float bntim_code)
{
  // FIXME, should go to I/O units further down in ggcm_mhd_bnd_fill_ghosts()
  float bntim = bntim_code * mhd->tnorm;
  if (mhd->amr == 0) {
    // FIXME, this really should be done in a cleaner way (pass mb, me, probably)
    int nr_comps = mrc_fld_nr_comps(fld);
    mrc_ddc_fill_ghosts_fld(mrc_domain_get_ddc(mhd->domain), 0, nr_comps, fld);
  } else {
    mrc_ddc_amr_apply(mhd->ddc_amr_cc, fld);
    // ggcm_mhd_amr_fill_ghosts_b(mhd, fld); // has been taken over by ddc_amr_cc
  }
  ggcm_mhd_bnd_fill_ghosts(mhd->bnd, fld, bntim);
  ggcm_mhd_bnd_fill_ghosts(mhd->bnd1, fld, bntim);
}

void
ggcm_mhd_fill_ghosts_E(struct ggcm_mhd *mhd, struct mrc_fld *E)
{
  // FIXME, also could do patch boundary ghost points / the AMR correction
  ggcm_mhd_bnd_fill_ghosts_E(mhd->bnd, E);
  ggcm_mhd_bnd_fill_ghosts_E(mhd->bnd1, E);
}

void
ggcm_mhd_fill_ghosts_reconstr(struct ggcm_mhd *mhd, struct mrc_fld *U_l[],
			      struct mrc_fld *U_r[], int p)
{
  ggcm_mhd_bnd_fill_ghosts_reconstr(mhd->bnd, U_l, U_r, p);
  ggcm_mhd_bnd_fill_ghosts_reconstr(mhd->bnd1, U_l, U_r, p);
}

void
ggcm_mhd_correct_fluxes(struct ggcm_mhd *mhd, struct mrc_fld *fluxes[3])
{
  if (mhd->amr > 0) {
    int gdims[3];
    mrc_domain_get_global_dims(mhd->domain, gdims);

    for (int d = 0; d < 3; d++) {
      if (gdims[d] > 1) {
	mrc_ddc_amr_apply(mhd->ddc_amr_flux[d], fluxes[d]);
      }
    }
  }
}

void
ggcm_mhd_correct_E(struct ggcm_mhd *mhd, struct mrc_fld *E)
{
  if (mhd->amr > 0) {
    mrc_ddc_amr_apply(mhd->ddc_amr_E, E);
  }
}

int
ggcm_mhd_ntot(struct ggcm_mhd *mhd)
{
  const int *ghost_dims = mrc_fld_ghost_dims(mhd->fld);
  return ghost_dims[0] * ghost_dims[1] * ghost_dims[2];
}

// ----------------------------------------------------------------------
// ggcm_mhd_get_3d_fld
//
// FIXME, this should cache the fields, rather than creating/destroying
// all the time

struct mrc_fld *
ggcm_mhd_get_3d_fld(struct ggcm_mhd *mhd, int nr_comps)
{
  struct mrc_fld *f = mrc_fld_create(ggcm_mhd_comm(mhd));
  mrc_fld_set_type(f , mrc_fld_type(mhd->fld));
  mrc_fld_set_param_obj(f, "domain", mhd->fld->_domain);
  mrc_fld_set_param_int(f, "nr_spatial_dims", 3);
  mrc_fld_set_param_int(f, "nr_comps", nr_comps);
  mrc_fld_set_param_int(f, "nr_ghosts", mhd->fld->_nr_ghosts);
  mrc_fld_set_param_bool(f, "aos", mhd->fld->_aos);
  mrc_fld_set_param_bool(f, "c_order", mhd->fld->_c_order);
  mrc_fld_setup(f);

  return f;
}

// ----------------------------------------------------------------------
// ggcm_mhd_put_3d_fld

void
ggcm_mhd_put_3d_fld(struct ggcm_mhd *mhd, struct mrc_fld *f)
{
  mrc_fld_destroy(f);
}

// ----------------------------------------------------------------------
// ggcm_mhd_default_box
//
// This function can be called in a subclass's ::create() function to
// set defaults for non-GGCM, normalized MHD-in-a-box simulations
//
// TODO: This should probably be the default in the first place

void
ggcm_mhd_default_box(struct ggcm_mhd *mhd)
{
  // use normalized units
  mhd->par.norm_mu0 = 1.f;

  mhd->par.diffco = 0.f;
  mhd->par.r_db_dt = 0.f;

  ggcm_mhd_set_param_float(mhd, "isphere", 0.);
  ggcm_mhd_set_param_float(mhd, "diffsphere", 0.);
  ggcm_mhd_set_param_float(mhd, "speedlimit", 1e9);

  // default to periodic boundary conditions
  ggcm_mhd_bnd_set_type(mhd->bnd, "none");
  ggcm_mhd_bnd_set_type(mhd->bnd1, "none");
  mrc_domain_set_param_int(mhd->domain, "bcx", BC_PERIODIC);
  mrc_domain_set_param_int(mhd->domain, "bcy", BC_PERIODIC);
  mrc_domain_set_param_int(mhd->domain, "bcz", BC_PERIODIC);

  mhd->par.gk_norm = true;
}

// ======================================================================
// ggcm_mhd class

static void
ggcm_mhd_init()
{
  mrc_class_register_subclass(&mrc_class_ggcm_mhd, &ggcm_mhd_ops_box);
}

static struct mrc_param_select magdiffu_descr[] = {
  { .val = MAGDIFFU_NL1  , .str = "nl1"     },
  { .val = MAGDIFFU_RES1 , .str = "res1"    },
  { .val = MAGDIFFU_CONST, .str = "const"   },
  {},
};

#define VAR(x) (void *)offsetof(struct ggcm_mhd, x)
static struct param ggcm_mhd_descr[] = {
  { "gamma"           , VAR(par.gamm)        , PARAM_FLOAT(1.66667f) },
  { "rrmin"           , VAR(par.rrmin)       , PARAM_FLOAT(.1f)      },

  { "xxnorm0"         , VAR(par.xxnorm0)     , PARAM_DOUBLE(1.)       },
  { "bbnorm0"         , VAR(par.bbnorm0)     , PARAM_DOUBLE(1.)       },
  { "vvnorm0"         , VAR(par.vvnorm0)     , PARAM_DOUBLE(1.)       },
  { "rrnorm0"         , VAR(par.rrnorm0)     , PARAM_DOUBLE(1.)       },
  { "ppnorm0"         , VAR(par.ppnorm0)     , PARAM_DOUBLE(1.)       },
  { "ccnorm0"         , VAR(par.ccnorm0)     , PARAM_DOUBLE(1.)       },
  { "eenorm0"         , VAR(par.eenorm0)     , PARAM_DOUBLE(1.)       },
  { "resnorm0"        , VAR(par.resnorm0)    , PARAM_DOUBLE(1.)       },
  { "tnorm0"          , VAR(par.tnorm0)      , PARAM_DOUBLE(1.)       },
  { "qqnorm0"         , VAR(par.qqnorm0)     , PARAM_DOUBLE(1.)       },

  { "norm_length"     , VAR(par.norm_length) , PARAM_DOUBLE(1.)       },
  { "norm_B"          , VAR(par.norm_B)      , PARAM_DOUBLE(1.)       },
  { "norm_density"    , VAR(par.norm_density), PARAM_DOUBLE(1.)       },
  { "norm_mu0"        , VAR(par.norm_mu0)    , PARAM_DOUBLE(C_MU0)    },
  { "mu0_code"        , VAR(par.mu0_code)    , PARAM_DOUBLE(1.)       },

  { "diffconstant"    , VAR(par.diffco)      , PARAM_FLOAT(.03f)     },
  { "diffthreshold"   , VAR(par.diffth)      , PARAM_FLOAT(.75f)     },
  { "diffsphere"      , VAR(par.diffsphere)  , PARAM_FLOAT(6.f)      },
  { "speedlimit"      , VAR(par.speedlimit)  , PARAM_FLOAT(1500.f)   },
  { "thx"             , VAR(par.thx)         , PARAM_FLOAT(.40f)     },
  { "isphere"         , VAR(par.isphere)     , PARAM_FLOAT(3.0f)     },
  { "r_db_dt"         , VAR(par.r_db_dt)     , PARAM_FLOAT(1.4f)     },
  { "timelo"          , VAR(par.timelo)      , PARAM_FLOAT(0.f)      },
  { "d_i"             , VAR(par.d_i)         , PARAM_FLOAT(0.f)      },
  { "dtmin"           , VAR(par.dtmin)       , PARAM_FLOAT(.0002f)   },
  { "modnewstep"      , VAR(par.modnewstep)  , PARAM_INT(1)          },
  { "magdiffu"        , VAR(par.magdiffu)    , PARAM_SELECT(MAGDIFFU_NL1,
							    magdiffu_descr) },
  { "diff_timelo"     , VAR(par.diff_timelo) , PARAM_FLOAT(0.)       },
  { "diff_swbnd"      , VAR(par.diff_swbnd)  , PARAM_FLOAT(-1e30)    },
  { "diff_obnd"       , VAR(par.diff_obnd)   , PARAM_INT(0)          },

  { "do_badval_checks"    , VAR(par.do_badval_checks)    , PARAM_BOOL(true)   },

  { "do_limit2"           , VAR(par.do_limit2)           , PARAM_BOOL(false)  },
  { "do_limit3"           , VAR(par.do_limit3)           , PARAM_BOOL(false)  },
  { "limit_aspect_low"    , VAR(par.limit_aspect_low)    , PARAM_BOOL(false)  },
  { "calce_aspect_low"    , VAR(par.calce_aspect_low)    , PARAM_BOOL(false)  },

  { "monitor_conservation", VAR(par.monitor_conservation), PARAM_BOOL(false)  },
  { "amr_grid_file"   , VAR(amr_grid_file)   , PARAM_STRING("amr_grid.txt")   },
  { "amr"             , VAR(amr)             , PARAM_INT(0)                   },

  { "xxnorm"          , VAR(xxnorm)          , MRC_VAR_DOUBLE         },
  { "bbnorm"          , VAR(bbnorm)          , MRC_VAR_DOUBLE         },
  { "vvnorm"          , VAR(vvnorm)          , MRC_VAR_DOUBLE         },
  { "rrnorm"          , VAR(rrnorm)          , MRC_VAR_DOUBLE         },
  { "ppnorm"          , VAR(ppnorm)          , MRC_VAR_DOUBLE         },
  { "ccnorm"          , VAR(ccnorm)          , MRC_VAR_DOUBLE         },
  { "eenorm"          , VAR(eenorm)          , MRC_VAR_DOUBLE         },
  { "resnorm"         , VAR(resnorm)         , MRC_VAR_DOUBLE         },
  { "tnorm"           , VAR(tnorm)           , MRC_VAR_DOUBLE         },

  { "time_code"       , VAR(time_code)       , MRC_VAR_FLOAT         },
  { "dt_code"         , VAR(dt_code)         , MRC_VAR_FLOAT         },
  { "istep"           , VAR(istep)           , MRC_VAR_INT           },
  { "timla"           , VAR(timla)           , MRC_VAR_FLOAT         },
  { "dacttime"        , VAR(dacttime)        , MRC_VAR_DOUBLE        },

  { "domain"          , VAR(domain)          , MRC_VAR_OBJ(mrc_domain)        },
  { "fld"             , VAR(fld)             , MRC_VAR_OBJ(mrc_fld)           },
  // FIXME, need to checkpoint b0, but only if !NULL
  //  { "b0"              , VAR(b0)              , MRC_VAR_OBJ(mrc_fld)           },
  { "crds"            , VAR(crds)            , MRC_VAR_OBJ(ggcm_mhd_crds)     },
  { "step"            , VAR(step)            , MRC_VAR_OBJ(ggcm_mhd_step)     },
  { "diag"            , VAR(diag)            , MRC_VAR_OBJ(ggcm_mhd_diag)     },
  { "bnd"             , VAR(bnd)             , MRC_VAR_OBJ(ggcm_mhd_bnd)      },
  { "bnd1"            , VAR(bnd1)            , MRC_VAR_OBJ(ggcm_mhd_bnd)      },
  { "ic"              , VAR(ic)              , MRC_VAR_OBJ(ggcm_mhd_ic)       },

  // gkeyll parameters // FIXME, use SI defaults
  { "gk_speed_of_light" , VAR(par.gk_speed_of_light) , PARAM_DOUBLE(1.f)              },
  { "gk_nr_fluids"      , VAR(par.gk_nr_fluids)      , PARAM_INT(2)                   },
  { "gk_nr_moments"     , VAR(par.gk_nr_moments)     , PARAM_INT(5)                   },
  { "gk_charge"         , VAR(par.gk_charge)         , PARAM_FLOAT_ARRAY(2, 1.f)      },
  { "gk_mass"           , VAR(par.gk_mass)           , PARAM_FLOAT_ARRAY(2, 1.f)      },
  { "gk_pressure_ratios", VAR(par.gk_pressure_ratios), PARAM_FLOAT_ARRAY(2, 1.f)      },

  { "gk_norm"                , VAR(par.gk_norm)                , PARAM_BOOL(false)    },
  { "gk_norm_speed_of_light" , VAR(par.gk_norm_speed_of_light) , PARAM_DOUBLE(20.)    }, // in units of vvnorm0
  { "gk_norm_mi_over_me"     , VAR(par.gk_norm_mi_over_me)     , PARAM_DOUBLE(25.)    },
  { "gk_norm_ppi_over_ppe"   , VAR(par.gk_norm_ppi_over_ppe)   , PARAM_DOUBLE(1.)     },
  { "gk_norm_rr"             , VAR(par.gk_norm_rr)             , PARAM_DOUBLE(1.)     }, // in units of rrnorm0

  {},
};
#undef VAR

struct mrc_class_ggcm_mhd mrc_class_ggcm_mhd = {
  .name             = "ggcm_mhd",
  .size             = sizeof(struct ggcm_mhd),
  .param_descr      = ggcm_mhd_descr,
  .init             = ggcm_mhd_init,
  .create           = _ggcm_mhd_create,
  .setup            = _ggcm_mhd_setup,
  .read             = _ggcm_mhd_read,
};

// ----------------------------------------------------------------------
// ggcm_mhd_wrongful_death
//
// Dump state and MPI_Abort, used for dtmin type crashes. This must be
// called from all mhd processes.

void
ggcm_mhd_wrongful_death(struct ggcm_mhd *mhd, struct mrc_fld *x, int errcode)
{
  static int cnt = 0;
  struct ggcm_mhd_diag *diag = ggcm_mhd_diag_create(ggcm_mhd_comm(mhd));
  // this extra diag will deadlock if using the fortran diag server...
  // FIXME: this is not a reliable way to fix the problem, but it works at the moment
  assert(strcmp(ggcm_mhd_diag_type(mhd->diag), "s2") != 0 &&
         strcmp(ggcm_mhd_diag_type(mhd->diag), "f2") != 0);
  ggcm_mhd_diag_set_type(diag, "c");
  ggcm_mhd_diag_set_param_obj(diag, "mhd", mhd);
  ggcm_mhd_diag_set_param_string(diag, "run", "wrongful_death");
  ggcm_mhd_diag_set_param_string(diag, "fields", "rr1:rv1:uu1:b1:rr:pp:v:b:j:e_cc:divb");
  ggcm_mhd_diag_setup(diag);
  // ggcm_mhd_diag_view(diag);
  
  mpi_printf(ggcm_mhd_comm(mhd), "Something bad happened. Dumping state then "
             "keeling over.\n");
  // ggcm_mhd_fill_ghosts(mhd, x, mhd->time);  // is this needed?
  ggcm_mhd_diag_run_now(diag, x, DIAG_TYPE_3D, cnt++);
  ggcm_mhd_diag_shutdown(diag);
  
  MPI_Barrier(ggcm_mhd_comm(mhd));
  MPI_Abort(MPI_COMM_WORLD, errcode);
}

// ----------------------------------------------------------------------
// ts_ggcm_mhd_step_calc_rhs
//
// wrapper to be used in a mrc_ts object

void
ts_ggcm_mhd_step_calc_rhs(void *ctx, struct mrc_obj *_rhs, float time, struct mrc_obj *_fld)
{
  struct ggcm_mhd *mhd = ctx;
  struct mrc_fld *rhs = (struct mrc_fld *) _rhs;
  struct mrc_fld *fld = (struct mrc_fld *) _fld;

  mhd->time_code = time;
  ggcm_mhd_step_calc_rhs(mhd->step, rhs, fld);
}

// ----------------------------------------------------------------------
// ts_ggcm_mhd_step_run
//
// wrapper to be used in a mrc_ts object

void
ts_ggcm_mhd_step_run(void *ctx, struct mrc_ts *ts, struct mrc_obj *_x)
{
  struct ggcm_mhd *mhd = ctx;
  struct mrc_fld *x = (struct mrc_fld *) _x;

  mhd->time_code = mrc_ts_time(ts);
  mhd->dt_code = mrc_ts_dt(ts);
  mhd->istep = mrc_ts_step_number(ts);

  ggcm_mhd_step_run(mhd->step, x);

  mrc_ts_set_dt(ts, mhd->dt_code);
}

// ----------------------------------------------------------------------
// ts_ggcm_mhd_step_get_dt
//
// wrapper to be used in a mrc_ts object

double
ts_ggcm_mhd_step_get_dt(void *ctx, struct mrc_ts *ts, struct mrc_obj *_x)
{
  struct ggcm_mhd *mhd = ctx;
  struct mrc_fld *x = (struct mrc_fld *) _x;

  mhd->time_code = mrc_ts_time(ts);
  mhd->dt_code = mrc_ts_dt(ts);
  mhd->istep = mrc_ts_step_number(ts);

  double dt = ggcm_mhd_step_get_dt(mhd->step, x);
  return dt;
}

// ----------------------------------------------------------------------
// ts_ggcm_mhd_pre_step
//
// mrc_ts wrapper for pre_step

static void
ts_ggcm_mhd_pre_step(void *ctx, struct mrc_ts *ts, struct mrc_obj *_x)
{
  struct ggcm_mhd *mhd = ctx;
  struct mrc_fld *x = (struct mrc_fld *) _x;

  mhd->time_code = mrc_ts_time(ts);
  mhd->dt_code = mrc_ts_dt(ts);
  mhd->istep = mrc_ts_step_number(ts);

  ggcm_mhd_pre_step(mhd, ts, x);
}

// ----------------------------------------------------------------------
// ts_ggcm_mhd_post_step
//
// mrc_ts wrapper for post_step

static void
ts_ggcm_mhd_post_step(void *ctx, struct mrc_ts *ts, struct mrc_obj *_x)
{
  struct ggcm_mhd *mhd = ctx;
  struct mrc_fld *x = (struct mrc_fld *) _x;

  ggcm_mhd_post_step(mhd, ts, x);
}

// ----------------------------------------------------------------------
// ggcm_mhd_setup_ts

void
ggcm_mhd_setup_ts(struct ggcm_mhd *mhd, struct mrc_ts *ts)
{
  mrc_ts_set_context(ts, ggcm_mhd_to_mrc_obj(mhd));
  mrc_ts_set_solution(ts, mrc_fld_to_mrc_obj(mhd->fld));
  mrc_ts_set_rhs_function(ts, ts_ggcm_mhd_step_calc_rhs, mhd);
  mrc_ts_set_step_function(ts, ts_ggcm_mhd_step_run, mhd);
  mrc_ts_set_get_dt_function(ts, ts_ggcm_mhd_step_get_dt, mhd);
  mrc_ts_set_pre_step_function(ts, ts_ggcm_mhd_pre_step, mhd);
  mrc_ts_set_post_step_function(ts, ts_ggcm_mhd_post_step, mhd);
}

// ----------------------------------------------------------------------
// ggcm_mhd_main
//
// Helper function that does most of the work of actually running a
// ggcm_mhd based simulation.
// The subclass of ggcm_mhd, and ggcm_mhd_ic are not explicitly set,
// so the default will be used unless overridden on the command line.
// This typically "just works", since the default will be the class you
// registered before calling this function.

int
ggcm_mhd_main(int *argc, char ***argv)
{
  MPI_Init(argc, argv);
  libmrc_params_init(*argc, *argv);
  ggcm_mhd_register();

  struct ggcm_mhd *mhd = ggcm_mhd_create(MPI_COMM_WORLD);
  ggcm_mhd_set_from_options(mhd);
  ggcm_mhd_setup(mhd);
  ggcm_mhd_view(mhd);

  // set up initial condition

  mpi_printf(MPI_COMM_WORLD, "Setting initial condition...\n");
  ggcm_mhd_ic_run(mhd->ic);
  
  // run time integration

  mpi_printf(MPI_COMM_WORLD, "Starting time integration...\n");
  double time_start = MPI_Wtime();

  struct mrc_ts *ts = mrc_ts_create(mrc_domain_comm(mhd->domain));
  if (ggcm_mhd_step_has_calc_rhs(mhd->step)) {
    mrc_ts_set_type(ts, "rk2");
  } else {
    mrc_ts_set_type(ts, "step");
  }

  struct mrc_ts_monitor *mon_output =
    mrc_ts_monitor_create(mrc_ts_comm(ts));
  mrc_ts_monitor_set_type(mon_output, "ggcm");
  mrc_ts_monitor_set_name(mon_output, "mrc_ts_output");
  mrc_ts_add_monitor(ts, mon_output);

  if (mhd->par.monitor_conservation) {
    struct mrc_ts_monitor *mon_conservation =
      mrc_ts_monitor_create(mrc_ts_comm(ts));
    mrc_ts_monitor_set_type(mon_conservation, "conservation");
    mrc_ts_monitor_set_name(mon_conservation, "mrc_ts_conservation");
    mrc_ts_add_monitor(ts, mon_conservation);
  }

  mrc_ts_set_param_double(ts, "norm_time", mhd->tnorm);
  mrc_ts_set_param_double(ts, "norm_time_scale", mhd->par.tnorm0);

  mrc_ts_set_dt(ts, 1e-6);
  ggcm_mhd_setup_ts(mhd, ts);
  mrc_ts_set_from_options(ts);
  mrc_ts_view(ts);
  mrc_ts_setup(ts);
  mrc_ts_solve(ts);
  mrc_ts_view(ts);

  double time_end = MPI_Wtime();

  int gdims[3];
  mrc_domain_get_global_dims(mhd->domain, gdims);
  int gsize = gdims[0] * gdims[1] * gdims[2];

  double cpu_time = time_end - time_start;
  mpi_printf(MPI_COMM_WORLD,"elapsed time = %g sec.\n", cpu_time);
  mpi_printf(MPI_COMM_WORLD,"\ncell-steps / second = %e\n",
	     (double) gsize * mrc_ts_step_number(ts) / cpu_time);

  mrc_ts_destroy(ts);  

  ggcm_mhd_destroy(mhd);

  MPI_Finalize();
  return 0;
}

