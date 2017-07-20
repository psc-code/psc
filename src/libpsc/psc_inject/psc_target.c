
#include "psc_target_private.h"

// ======================================================================
// psc_target

// ----------------------------------------------------------------------
// psc_target_is_inside

bool
psc_target_is_inside(struct psc_target *target, double x[3])
{
  return (x[1] >= target->yl && x[1] <= target->yh &&
	  x[2] >= target->zl && x[2] <= target->zh);
}

// ----------------------------------------------------------------------
// psc_target_init_npt

void
psc_target_init_npt(struct psc_target *target, int pop, double x[3],
		    struct psc_particle_npt *npt)
{
  if (!psc_target_is_inside(target, x)) {
    npt->n = 0;
    return;
  }

  if (pop == target->kind_ion) {
    npt->n    = target->n;
    npt->T[0] = target->Ti;
    npt->T[1] = target->Ti;
    npt->T[2] = target->Ti;
  } else if (pop == target->kind_electron) {
    npt->n    = target->n;
    npt->T[0] = target->Te;
    npt->T[1] = target->Te;
    npt->T[2] = target->Te;
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// psc_target class

struct mrc_class_psc_target mrc_class_psc_target = {
  .name             = "psc_target",
  .size             = sizeof(struct psc_target),
  .param_descr      = psc_target_descr,
};

