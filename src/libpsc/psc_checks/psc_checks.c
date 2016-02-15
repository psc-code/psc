
#include "psc_checks_private.h"

// ----------------------------------------------------------------------
// psc_checks_continuity_before_particle_push

void
psc_checks_continuity_before_particle_push(struct psc_checks *checks, struct psc *psc)
{
  struct psc_checks_ops *ops = psc_checks_ops(checks);

  assert(ops && ops->continuity_before_particle_push);
  ops->continuity_before_particle_push(checks, psc);
}

// ----------------------------------------------------------------------
// psc_checks_continuity_after_particle_push

void
psc_checks_continuity_after_particle_push(struct psc_checks *checks, struct psc *psc)
{
  struct psc_checks_ops *ops = psc_checks_ops(checks);

  assert(ops && ops->continuity_after_particle_push);
  ops->continuity_after_particle_push(checks, psc);
}

// ----------------------------------------------------------------------
// psc_checks_gauss

void
psc_checks_gauss(struct psc_checks *checks, struct psc *psc)
{
  struct psc_checks_ops *ops = psc_checks_ops(checks);

  assert(ops && ops->gauss);
  ops->gauss(checks, psc);
}

// ----------------------------------------------------------------------
// psc_checks_init

static void
psc_checks_init()
{
  mrc_class_register_subclass(&mrc_class_psc_checks, &psc_checks_1st_double_ops);
  mrc_class_register_subclass(&mrc_class_psc_checks, &psc_checks_1st_single_ops);
}

// ----------------------------------------------------------------------
// psc_checks_descr

#define VAR(x) (void *)offsetof(struct psc_checks, x)

static struct param psc_checks_descr[] = {
  { "continuity_every_step" , VAR(continuity_every_step) , PARAM_INT(-1)       },
  { "continuity_threshold"  , VAR(continuity_threshold)  , PARAM_DOUBLE(1e-14) },
  { "continuity_verbose"    , VAR(continuity_verbose)    , PARAM_BOOL(false)   },
  { "continuity_dump_always", VAR(continuity_dump_always), PARAM_BOOL(false)   },

  { "gauss_every_step"      , VAR(gauss_every_step)      , PARAM_INT(-1)       },
  { "gauss_threshold"       , VAR(gauss_threshold)       , PARAM_DOUBLE(1e-14) },
  { "gauss_verbose"         , VAR(gauss_verbose)         , PARAM_BOOL(false)   },
  { "gauss_dump_always"     , VAR(gauss_dump_always)     , PARAM_BOOL(false)   },

  { "rho_m"                 , VAR(rho_m)                 , MRC_VAR_OBJ(psc_mfields) },
  { "rho_p"                 , VAR(rho_p)                 , MRC_VAR_OBJ(psc_mfields) },
  {},
};

#undef VAR

// ----------------------------------------------------------------------
// psc_checks class

struct mrc_class_psc_checks mrc_class_psc_checks = {
  .name             = "psc_checks",
  .size             = sizeof(struct psc_checks),
  .param_descr      = psc_checks_descr,
  .init             = psc_checks_init,
};


