
#include "ggcm_mhd_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_defs_extra.h"
#include "ggcm_mhd_crds.h"
#include "ggcm_mhd_crds_private.h"
#include "ggcm_mhd_crds_gen.h"
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

#include <assert.h>
#include <string.h>

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
// ggcm_mhd_get_state
//
// update C ggcm_mhd state from Fortran common blocks

void
ggcm_mhd_get_state(struct ggcm_mhd *mhd)
{
  struct ggcm_mhd_ops *ops = ggcm_mhd_ops(mhd);
  if (ops->get_state) {
    ops->get_state(mhd);
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_set_state
//
// updated Fortran common blocks from C ggcm_mhd state

void
ggcm_mhd_set_state(struct ggcm_mhd *mhd)
{
  struct ggcm_mhd_ops *ops = ggcm_mhd_ops(mhd);
  if (ops->set_state) {
    ops->set_state(mhd);
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
  if (mhd->fld->_nr_spatial_dims == 3 
   && mhd->fld->_is_aos) shift = 1;
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

static void
_ggcm_mhd_setup(struct ggcm_mhd *mhd)
{
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
ggcm_mhd_fill_ghosts(struct ggcm_mhd *mhd, struct mrc_fld *fld, int m, float bntim)
{
  if (mhd->amr == 0) {
    // FIXME, this really should be done in a cleaner way (pass mb, me, probably)
    int mhd_type;
    mrc_fld_get_param_int(fld, "mhd_type", &mhd_type);
    if (mhd_type == MT_GKEYLL) {
      int nr_comps = mrc_fld_nr_comps(fld);
      mrc_ddc_fill_ghosts_fld(mrc_domain_get_ddc(mhd->domain), m, m + nr_comps, fld);
    } else {
      mrc_ddc_fill_ghosts_fld(mrc_domain_get_ddc(mhd->domain), m, m + 8, fld);
    }
  } else {
    assert(m == 0);
    mrc_ddc_amr_apply(mhd->ddc_amr_cc, fld);
    // ggcm_mhd_amr_fill_ghosts_b(mhd, fld); // has been taken over by ddc_amr_cc
  }
  ggcm_mhd_bnd_fill_ghosts(mhd->bnd, fld, m, bntim);
  ggcm_mhd_bnd_fill_ghosts(mhd->bnd1, fld, m, bntim);
}

void
ggcm_mhd_fill_ghosts_E(struct ggcm_mhd *mhd, struct mrc_fld *E)
{
  struct ggcm_mhd_ops *ops = ggcm_mhd_ops(mhd);
  // FIXME, should this be done via ggcm_mhd_bnd instead through ggcm_mhd itself?
  // FIXME, also could do patch boundary ghost points / the AMR correction
  if (ops->fill_ghosts_E) {
    ops->fill_ghosts_E(mhd, E);
  }
  ggcm_mhd_bnd_fill_ghosts_E(mhd->bnd, E);
  ggcm_mhd_bnd_fill_ghosts_E(mhd->bnd1, E);
}

int
ggcm_mhd_ntot(struct ggcm_mhd *mhd)
{
  const int *ghost_dims = mrc_fld_ghost_dims(mhd->fld);
  return ghost_dims[0] * ghost_dims[1] * ghost_dims[2];
}

// ----------------------------------------------------------------------
// ggcm_mhd_get_crds_cc

void
ggcm_mhd_get_crds_cc(struct ggcm_mhd *mhd, int ix, int iy, int iz, int p,
		     float crd_cc[3])
{
  struct mrc_crds *crds = mrc_domain_get_crds(mhd->domain);
  crd_cc[0] = MRC_MCRDX(crds, ix, p);
  crd_cc[1] = MRC_MCRDY(crds, iy, p);
  crd_cc[2] = MRC_MCRDZ(crds, iz, p);
}

// ----------------------------------------------------------------------
// ggcm_mhd_get_crds_nc

void
ggcm_mhd_get_crds_nc(struct ggcm_mhd *mhd, int ix, int iy, int iz, int p,
		     float crd_nc[3])
{
  float crd_cc[3], crd_cc_m[3];
  ggcm_mhd_get_crds_cc(mhd, ix  ,iy  ,iz  , p, crd_cc);
  ggcm_mhd_get_crds_cc(mhd, ix-1,iy-1,iz-1, p, crd_cc_m);

  crd_nc[0] = .5f * (crd_cc_m[0] + crd_cc[0]);
  crd_nc[1] = .5f * (crd_cc_m[1] + crd_cc[1]);
  crd_nc[2] = .5f * (crd_cc_m[2] + crd_cc[2]);
}

// ----------------------------------------------------------------------
// ggcm_mhd_get_crds_fc

void
ggcm_mhd_get_crds_fc(struct ggcm_mhd *mhd, int ix, int iy, int iz, int p,
		     int d, float crd_fc[3])
{
  float crd_cc[3];
  float crd_nc[3];
  ggcm_mhd_get_crds_cc(mhd, ix, iy, iz, p, crd_cc);
  ggcm_mhd_get_crds_nc(mhd, ix, iy, iz, p, crd_nc);

  if (d == 0) {
    // Bx located at i, j+.5, k+.5
    crd_fc[0] = crd_nc[0];
    crd_fc[1] = crd_cc[1];
    crd_fc[2] = crd_cc[2];
  } else if (d == 1) {
    // By located at i+.5, j, k+.5
    crd_fc[0] = crd_cc[0];
    crd_fc[1] = crd_nc[1];
    crd_fc[2] = crd_cc[2];
  } else if (d == 2) {
    // Bz located at i+.5, j+.5, k
    crd_fc[0] = crd_cc[0];
    crd_fc[1] = crd_cc[1];
    crd_fc[2] = crd_nc[2];
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_get_crds_ec

void
ggcm_mhd_get_crds_ec(struct ggcm_mhd *mhd, int ix, int iy, int iz, int p,
		     int d, float crd_ec[3])
{
  float crd_cc[3], crd_nc[3];
  ggcm_mhd_get_crds_cc(mhd, ix,iy,iz, p, crd_cc);
  ggcm_mhd_get_crds_nc(mhd, ix,iy,iz, p, crd_nc);

  if (d == 0) {
    // Ex located at i+.5, j, k
    crd_ec[0] = crd_cc[0];
    crd_ec[1] = crd_nc[1];
    crd_ec[2] = crd_nc[2];
  } else if (d == 1) {
    // Ey located at i, j+.5, k
    crd_ec[0] = crd_nc[0];
    crd_ec[1] = crd_cc[1];
    crd_ec[2] = crd_nc[2];
  } else if (d == 2) {
    // Ez located at i, j, k+.5
    crd_ec[0] = crd_nc[0];
    crd_ec[1] = crd_nc[1];
    crd_ec[2] = crd_cc[2];
  } else {
    assert(0);
  }
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
  mhd->par.rrnorm = 1.f;
  mhd->par.ppnorm = 1.f;
  mhd->par.vvnorm = 1.f;
  mhd->par.bbnorm = 1.f;
  mhd->par.ccnorm = 1.f;
  mhd->par.eenorm = 1.f;
  mhd->par.resnorm = 1.f;
  mhd->par.tnorm = 1.f;
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

  // generate MHD solver grid from mrc_crds
  ggcm_mhd_crds_gen_set_type(mhd->crds->crds_gen, "mrc");
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
  { "bbnorm"          , VAR(par.bbnorm)      , PARAM_FLOAT(30574.f)  },
  { "vvnorm"          , VAR(par.vvnorm)      , PARAM_FLOAT(6692.98f) },
  { "rrnorm"          , VAR(par.rrnorm)      , PARAM_FLOAT(10000.f)  },
  { "ppnorm"          , VAR(par.ppnorm)      , PARAM_FLOAT(7.43866e8)},
  { "ccnorm"          , VAR(par.ccnorm)      , PARAM_FLOAT(3.81885)  },
  { "eenorm"          , VAR(par.eenorm)      , PARAM_FLOAT(204631.f) },
  { "resnorm"         , VAR(par.resnorm)     , PARAM_FLOAT(53.5848e6)},
  { "tnorm"           , VAR(par.tnorm)       , PARAM_FLOAT(.95189935)},
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
  { "dbasetime"       , VAR(par.dbasetime)   , PARAM_DOUBLE(0.)      },
  { "modnewstep"      , VAR(par.modnewstep)  , PARAM_INT(1)          },
  { "magdiffu"        , VAR(par.magdiffu)    , PARAM_SELECT(MAGDIFFU_NL1,
							    magdiffu_descr) },
  { "diff_timelo"     , VAR(par.diff_timelo) , PARAM_FLOAT(0.)       },
  { "diff_swbnd"      , VAR(par.diff_swbnd)  , PARAM_FLOAT(-1e30)    },
  { "diff_obnd"       , VAR(par.diff_obnd)   , PARAM_INT(0)          },

  { "monitor_conservation", VAR(par.monitor_conservation), PARAM_BOOL(false)  },
  { "amr_grid_file"   , VAR(amr_grid_file)   , PARAM_STRING("amr_grid.txt")   },
  { "amr"             , VAR(amr)             , PARAM_INT(0)                   },

  { "time"            , VAR(time)            , MRC_VAR_FLOAT         },
  { "dt"              , VAR(dt)              , MRC_VAR_FLOAT         },
  { "istep"           , VAR(istep)           , MRC_VAR_INT           },
  { "timla"           , VAR(timla)           , MRC_VAR_FLOAT         },
  { "dacttime"        , VAR(dacttime)        , MRC_VAR_DOUBLE        },
  { "max_time"        , VAR(max_time)        , PARAM_FLOAT(0.0)      },
  
  { "domain"          , VAR(domain)          , MRC_VAR_OBJ(mrc_domain)        },
  { "fld"             , VAR(fld)             , MRC_VAR_OBJ(mrc_fld)           },
  { "crds"            , VAR(crds)            , MRC_VAR_OBJ(ggcm_mhd_crds)     },
  { "step"            , VAR(step)            , MRC_VAR_OBJ(ggcm_mhd_step)     },
  { "diag"            , VAR(diag)            , MRC_VAR_OBJ(ggcm_mhd_diag)     },
  { "bnd"             , VAR(bnd)             , MRC_VAR_OBJ(ggcm_mhd_bnd)      },
  { "bnd1"            , VAR(bnd1)            , MRC_VAR_OBJ(ggcm_mhd_bnd)      },
  { "ic"              , VAR(ic)              , MRC_VAR_OBJ(ggcm_mhd_ic)       },

  { "do_badval_checks", VAR(do_badval_checks), PARAM_BOOL(true)             },
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
  // ggcm_mhd_fill_ghosts(mhd, x, 0, mhd->time);  // is this needed?
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
  
  mhd->time = time;
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
  
  mhd->time = mrc_ts_time(ts);
  mhd->dt = mrc_ts_dt(ts);
  mhd->istep = mrc_ts_step_number(ts);

  ggcm_mhd_step_run(mhd->step, x);

  mrc_ts_set_dt(ts, mhd->dt);
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
  mrc_ts_set_context(ts, ggcm_mhd_to_mrc_obj(mhd));

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

  mrc_ts_set_dt(ts, 1e-6);
  mrc_ts_set_solution(ts, mrc_fld_to_mrc_obj(mhd->fld));
  mrc_ts_set_rhs_function(ts, ts_ggcm_mhd_step_calc_rhs, mhd);
  mrc_ts_set_step_function(ts, ts_ggcm_mhd_step_run, mhd);
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

