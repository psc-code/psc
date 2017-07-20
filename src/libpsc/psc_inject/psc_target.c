
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
// psc_target class

struct mrc_class_psc_target mrc_class_psc_target = {
  .name             = "psc_target",
  .size             = sizeof(struct psc_target),
  .param_descr      = psc_target_descr,
};

