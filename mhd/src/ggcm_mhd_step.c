
#include "ggcm_mhd_step_private.h"

#include "ggcm_mhd_defs.h"
#include "ggcm_mhd_defs_extra.h"
#include "ggcm_mhd_private.h"
#include "ggcm_mhd_diag_private.h"
#include "ggcm_mhd_diag_item_private.h"

#include <mrc_io.h>
#include <mrc_profile.h>
#include <mrc_bits.h>
#include <math.h>
#include <assert.h>
#include <mpi.h>

// ======================================================================
// ggcm_mhd_step class

// ----------------------------------------------------------------------
// ggcm_mhd_step_get_dt

double
ggcm_mhd_step_get_dt(struct ggcm_mhd_step *step, struct mrc_fld *x)
{
  struct ggcm_mhd *mhd = step->mhd;

  struct ggcm_mhd_step_ops *ops = ggcm_mhd_step_ops(step);
  if (!ops->get_dt) {
    return mhd->dt;
  }

  if (!step->legacy_dt_handling) {
    double dtn = ops->get_dt(step, x);

    if (mhd->istep > 0 && dtn < mhd->par.dtmin) {
      mpi_printf(ggcm_mhd_comm(mhd), "!!! dt < dtmin. Dying now!\n");
      mpi_printf(ggcm_mhd_comm(mhd), "!!! dt = %g (was %g), dtmin = %g\n",
		 dtn, mhd->dt, mhd->par.dtmin);
      ggcm_mhd_wrongful_death(mhd, mhd->fld, -1);
    }
    
    mhd->dt = dtn;
  } else { // legacy_dt_handling
    if (step->do_nwst) {
      step->dtn = ops->get_dt(step, x);
      // yes, dtn isn't set to mhd->dt until the end of the step... this
      // is what the fortran code did
    }
  }

  return mhd->dt;
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_calc_rhs

void
ggcm_mhd_step_calc_rhs(struct ggcm_mhd_step *step, struct mrc_fld *rhs,
		       struct mrc_fld *x)
{
  struct ggcm_mhd_step_ops *ops = ggcm_mhd_step_ops(step);
  assert(ops && ops->calc_rhs);
  ops->calc_rhs(step, rhs, x);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_get_e_ec

void
ggcm_mhd_step_get_e_ec(struct ggcm_mhd_step *step, struct mrc_fld *E,
                        struct mrc_fld *x)
{
  struct ggcm_mhd_step_ops *ops = ggcm_mhd_step_ops(step);
  assert(ops && ops->get_e_ec);
  ops->get_e_ec(step, E, x);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_run

void
ggcm_mhd_step_run(struct ggcm_mhd_step *step, struct mrc_fld *x)
{
  struct ggcm_mhd_step_ops *ops = ggcm_mhd_step_ops(step);
  static int pr;
  if (!pr) {
    pr = prof_register("ggcm_mhd_step_run", 0, 0, 0.);
  }

  if (step->debug_dump) {
    struct ggcm_mhd *mhd = step->mhd;
    static struct ggcm_mhd_diag *diag;
    static int cnt;
    if (!diag) {
      diag = ggcm_mhd_diag_create(ggcm_mhd_comm(mhd));
      ggcm_mhd_diag_set_type(diag, "c");
      ggcm_mhd_diag_set_name(diag, "ggcm_mhd_debug");
      ggcm_mhd_diag_set_param_obj(diag, "mhd", mhd);
      ggcm_mhd_diag_set_param_string(diag, "fields", "rr1:rv1:uu1:b1:rr:v:pp:b:divb");
      ggcm_mhd_diag_set_from_options(diag);
      // overwrite run name (needs to be after set_from_options)
      ggcm_mhd_diag_set_param_string(diag, "run", "dbg");
      ggcm_mhd_diag_setup(diag);
      ggcm_mhd_diag_view(diag);
    }
    ggcm_mhd_fill_ghosts(mhd, mhd->fld, 0, mhd->time);
    ggcm_mhd_diag_run_now(diag, mhd->fld, DIAG_TYPE_3D, cnt++);
  }

  prof_start(pr);
  assert(ops && ops->run);
  ops->run(step, x);
  prof_stop(pr);

  if (step->legacy_dt_handling && ops->get_dt) {
    struct ggcm_mhd *mhd = step->mhd;
    float dtn = fminf(1., step->dtn);

    if (dtn > 1.02 * mhd->dt || dtn < mhd->dt / 1.01) {
      mpi_printf(ggcm_mhd_comm(mhd), "switched dt %g <- %g\n", dtn, mhd->dt);
      
      if (mhd->istep > 0 &&
          (dtn < 0.5 * mhd->dt || dtn > 2.0 * mhd->dt)) {            
        mpi_printf(ggcm_mhd_comm(mhd), "!!! dt changed by > a factor of 2. "
                   "Dying now!\n");
        ggcm_mhd_wrongful_death(mhd, mhd->fld, 2);
      }
      
      mhd->dt = dtn;
    }
  }

  // FIXME, this should be done by mrc_ts
  struct ggcm_mhd *mhd = step->mhd;
  if ((mhd->istep % step->profile_every) == 0) {
    prof_print_mpi(ggcm_mhd_comm(mhd));
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_run_predcorr
//
// library-type function to be used by ggcm_mhd_step subclasses that
// implement the OpenGGCM predictor-corrector scheme

void
ggcm_mhd_step_run_predcorr(struct ggcm_mhd_step *step, struct mrc_fld *x)
{
  static int PR_push;
  if (!PR_push) {
    PR_push = prof_register("ggcm_mhd_step_run_predcorr", 1., 0, 0);
  }

  prof_start(PR_push);

  struct ggcm_mhd_step_ops *ops = ggcm_mhd_step_ops(step);
  struct ggcm_mhd *mhd = step->mhd;

  float dtn;
  if (step->do_nwst) {
    assert(ops && ops->newstep);
    ops->newstep(step, &dtn);
    // yes, dtn isn't set to mhd->dt until the end of the step... this
    // is what the fortran code did
  }

  ggcm_mhd_fill_ghosts(mhd, x, _RR1, mhd->time);
  assert(ops && ops->pred);
  ops->pred(step);

  ggcm_mhd_fill_ghosts(mhd, x, _RR2, mhd->time + mhd->bndt);
  assert(ops && ops->corr);
  ops->corr(step);

  if (step->do_nwst) {
    if (step->legacy_dt_handling) {
      dtn = fminf(1., dtn); // FIXME, only kept for compatibility
    }

    if (dtn > 1.02 * mhd->dt || dtn < mhd->dt / 1.01) {
      mpi_printf(ggcm_mhd_comm(mhd), "switched dt %g <- %g\n", dtn, mhd->dt);
      
      // FIXME: determining when to die on a bad dt should be generalized, since
      //        there's another hiccup if refining dt for actual AMR      
      bool first_step = mhd->istep <= 1;
      bool last_step = mhd->time + dtn > (1.0 - 1e-5) * mhd->max_time;

      if (!first_step && !last_step &&
          (dtn < 0.5 * mhd->dt || dtn > 2.0 * mhd->dt)) {            
        mpi_printf(ggcm_mhd_comm(mhd), "!!! dt changed by > a factor of 2. "
                   "Dying now!\n");
        ggcm_mhd_wrongful_death(mhd, mhd->fld, 2);
      }
      
      mhd->dt = dtn;
    }
  }

  prof_stop(PR_push);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_has_calc_rhs

bool
ggcm_mhd_step_has_calc_rhs(struct ggcm_mhd_step *step)
{
  struct ggcm_mhd_step_ops *ops = ggcm_mhd_step_ops(step);
  assert(ops);
  return ops->calc_rhs;
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_setup_flds

void
ggcm_mhd_step_setup_flds(struct ggcm_mhd_step *step)
{
  struct ggcm_mhd_step_ops *ops = ggcm_mhd_step_ops(step);
  assert(ops && ops->setup_flds);
  return ops->setup_flds(step);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_get_1d_fld
//

struct mrc_fld *
ggcm_mhd_step_get_1d_fld(struct ggcm_mhd_step *step, int nr_comps)
{
  int c = nr_comps - 1;
  struct mrc_fld_cache *cache = &step->cache_1d[c];
  if (cache->n > 0) {
    return cache->flds[--cache->n];
  }

  struct mrc_fld *fld = step->mhd->fld;

  struct mrc_fld *f = mrc_fld_create(ggcm_mhd_step_comm(step));
  int size = 0;
  for (int d = 0; d < 3; d++) {
    size = MAX(size, mrc_fld_dims(fld)[d + 1]); // FIXME assumes aos
  }

  mrc_fld_set_type(f, mrc_fld_type(fld));
  mrc_fld_set_param_int_array(f, "dims", 2, (int []) { nr_comps, size });
  mrc_fld_set_param_int_array(f, "sw"  , 2, (int []) { 0, fld->_nr_ghosts });
  mrc_fld_setup(f);

  return f;
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_put_1d_fld

void
ggcm_mhd_step_put_1d_fld(struct ggcm_mhd_step *step, struct mrc_fld *f)
{
  int c = mrc_fld_dims(f)[0] - 1;
  struct mrc_fld_cache *cache = &step->cache_1d[c];
  assert(cache->n < MRC_FLD_CACHE_SIZE);
  cache->flds[cache->n++] = f;
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_destroy

static void
_ggcm_mhd_step_destroy(struct ggcm_mhd_step *step)
{
  for (int c = 0; c < MRC_FLD_CACHE_COMPS; c++) {
    struct mrc_fld_cache *cache = &step->cache_1d[c];
    while (cache->n > 0) {
      mrc_fld_destroy(cache->flds[--cache->n]);
    }
  }
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_diag_item_zmask_run

static void
ggcm_mhd_step_diag_item_zmask_run(struct ggcm_mhd_step *step,
				  struct ggcm_mhd_diag_item *item,
				  struct mrc_io *io, struct mrc_fld *f,
				  int diag_type, float plane)
{
  struct ggcm_mhd_step_ops *ops = ggcm_mhd_step_ops(step);
  assert(ops && ops->diag_item_zmask_run);
  ops->diag_item_zmask_run(step, item, io, f, diag_type, plane);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_diag_item_rmask_run

static void
ggcm_mhd_step_diag_item_rmask_run(struct ggcm_mhd_step *step,
				  struct ggcm_mhd_diag_item *item,
				  struct mrc_io *io, struct mrc_fld *f,
				  int diag_type, float plane)
{
  struct ggcm_mhd_step_ops *ops = ggcm_mhd_step_ops(step);
  assert(ops && ops->diag_item_rmask_run);
  ops->diag_item_rmask_run(step, item, io, f, diag_type, plane);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step_init

static void
ggcm_mhd_step_init()
{
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_step, &ggcm_mhd_step_c3_float_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_step, &ggcm_mhd_step_c3_double_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_step, &ggcm_mhd_step_cweno_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_step, &ggcm_mhd_step_c_float_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_step, &ggcm_mhd_step_c_double_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_step, &ggcm_mhd_step_c2_float_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_step, &ggcm_mhd_step_c2_double_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_step, &ggcm_mhd_step_vlct_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_step, &ggcm_mhd_step_vl_ops);
  mrc_class_register_subclass(&mrc_class_ggcm_mhd_step, &ggcm_mhd_step_mhdcc_double_ops);
}

// ----------------------------------------------------------------------
// ggcm_mhd_step description

#define VAR(x) (void *)offsetof(struct ggcm_mhd_step, x)
static struct param ggcm_mhd_step_descr[] = {
  { "mhd"               , VAR(mhd)               , PARAM_OBJ(ggcm_mhd)     },
  // This is a bit hacky, do_nwst may normally be set in the timeloop
  // to determine whether to run newstep() the next timestep or not,
  // but this allows to set it to "always on" easily for test runs
  { "do_nwst"           , VAR(do_nwst)           , PARAM_BOOL(false)       },
  { "debug_dump"        , VAR(debug_dump)        , PARAM_BOOL(false)       },
  { "profile_every"     , VAR(profile_every)     , PARAM_INT(10)           },
  // FIXME, is there really a need to keep this around?
  // It limits the (normalized) to stay <= 1, and it should also wrap the
  // "find new dt first -- but do this timestep still with old dt"
  // weirdness
  { "legacy_dt_handling", VAR(legacy_dt_handling), PARAM_BOOL(true)        },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// ggcm_mhd_step class description

struct mrc_class_ggcm_mhd_step mrc_class_ggcm_mhd_step = {
  .name             = "ggcm_mhd_step",
  .size             = sizeof(struct ggcm_mhd_step),
  .param_descr      = ggcm_mhd_step_descr,
  .init             = ggcm_mhd_step_init,
  .destroy          = _ggcm_mhd_step_destroy,
};

/////////////////////////////////////////////////////////////////////////
// diag items to go with specific ggcm_mhd_step subclasses

// ======================================================================
// ggcm_mhd_diag_item subclass "zmask"

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item_zmask_run

static void
ggcm_mhd_diag_item_zmask_run(struct ggcm_mhd_diag_item *item,
			   struct mrc_io *io, struct mrc_fld *f,
			   int diag_type, float plane)
{
  struct ggcm_mhd_step *step = item->diag->mhd->step;
  ggcm_mhd_step_diag_item_zmask_run(step, item, io, f, diag_type, plane);
}

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item subclass "zmask"

struct ggcm_mhd_diag_item_ops ggcm_mhd_diag_item_ops_zmask = {
  .name             = "zmask",
  .run              = ggcm_mhd_diag_item_zmask_run,
};

// ======================================================================
// ggcm_mhd_diag_item subclass "rmask"

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item_rmask_run

static void
ggcm_mhd_diag_item_rmask_run(struct ggcm_mhd_diag_item *item,
			     struct mrc_io *io, struct mrc_fld *f,
			     int diag_type, float plane)
{
  struct ggcm_mhd_step *step = item->diag->mhd->step;
  ggcm_mhd_step_diag_item_rmask_run(step, item, io, f, diag_type, plane);
}

// ----------------------------------------------------------------------
// ggcm_mhd_diag_item subclass "rmask"

struct ggcm_mhd_diag_item_ops ggcm_mhd_diag_item_ops_rmask = {
  .name             = "rmask",
  .run              = ggcm_mhd_diag_item_rmask_run,
};

