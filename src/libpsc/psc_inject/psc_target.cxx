
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
// psc_target class

struct mrc_class_psc_target_ : mrc_class_psc_target {
  mrc_class_psc_target_() {
    name             = "psc_target";
    size             = sizeof(struct psc_target);
  }
} mrc_class_psc_target;

