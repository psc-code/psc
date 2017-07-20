
#include "psc_target_private.h"

// ======================================================================
// psc_target

// ----------------------------------------------------------------------
// psc_target_is_inside

bool
psc_target_is_inside(struct psc_target *target, double x[3])
{
  struct psc_target_ops *ops = psc_target_ops(target);

  assert(ops && ops->is_inside);
  return ops->is_inside(target, x);
}

// ----------------------------------------------------------------------
// psc_target_init_npt

void
psc_target_init_npt(struct psc_target *target, int pop, double x[3],
		    struct psc_particle_npt *npt)
{
  struct psc_target_ops *ops = psc_target_ops(target);

  assert(ops && ops->init_npt);
  return ops->init_npt(target, pop, x, npt);
}

// ----------------------------------------------------------------------
// psc_target_init

extern struct psc_target_ops psc_target_ops_slab;

static void
psc_target_init()
{
  mrc_class_register_subclass(&mrc_class_psc_target, &psc_target_ops_slab);
}

// ----------------------------------------------------------------------
// psc_target class

struct mrc_class_psc_target mrc_class_psc_target = {
  .name             = "psc_target",
  .size             = sizeof(struct psc_target),
  .init             = psc_target_init,
};

// ======================================================================
// psc_target subclass "slab"

struct psc_target_slab {
  // params
  double yl;
  double yh;
  double zl;
  double zh;
  double n;
  double Te;
  double Ti;
  int kind_ion;
  int kind_electron;
};

#define psc_target_slab(target) mrc_to_subobj(target, struct psc_target_slab)

#define VAR(x) (void *)offsetof(struct psc_target_slab, x)
static struct param psc_target_slab_descr[] _mrc_unused = {
  { "yl"           , VAR(yl)           , PARAM_DOUBLE(0.)       },
  { "yh"           , VAR(yh)           , PARAM_DOUBLE(0.)       },
  { "zl"           , VAR(zl)           , PARAM_DOUBLE(0.)       },
  { "zh"           , VAR(zh)           , PARAM_DOUBLE(0.)       },
  { "n"            , VAR(n)            , PARAM_DOUBLE(1.)       },
  { "Te"           , VAR(Te)           , PARAM_DOUBLE(.001)     },
  { "Ti"           , VAR(Ti)           , PARAM_DOUBLE(.001)     },
  { "kind_ion"     , VAR(kind_ion)     , PARAM_INT(-1)          },
  { "kind_electron", VAR(kind_electron), PARAM_INT(-1)          },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// psc_target_slab_is_inside

static bool
psc_target_slab_is_inside(struct psc_target *target, double x[3])
{
  struct psc_target_slab *sub = psc_target_slab(target);
  
  return (x[1] >= sub->yl && x[1] <= sub->yh &&
	  x[2] >= sub->zl && x[2] <= sub->zh);
}

// ----------------------------------------------------------------------
// psc_target_slab_init_npt

static void
psc_target_slab_init_npt(struct psc_target *target, int pop, double x[3],
			 struct psc_particle_npt *npt)
{
  struct psc_target_slab *sub = psc_target_slab(target);

  if (!psc_target_slab_is_inside(target, x)) {
    npt->n = 0;
    return;
  }

  if (pop == sub->kind_ion) {
    npt->n    = sub->n;
    npt->T[0] = sub->Ti;
    npt->T[1] = sub->Ti;
    npt->T[2] = sub->Ti;
  } else if (pop == sub->kind_electron) {
    npt->n    = sub->n;
    npt->T[0] = sub->Te;
    npt->T[1] = sub->Te;
    npt->T[2] = sub->Te;
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// psc_target "slab"

struct psc_target_ops psc_target_ops_slab = {
  .name                = "slab",
  .size                = sizeof(struct psc_target_slab),
  .param_descr         = psc_target_slab_descr,
  .is_inside           = psc_target_slab_is_inside,
  .init_npt            = psc_target_slab_init_npt,
};
