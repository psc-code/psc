
#include "psc_pulse_private.h"

#include <math.h>

static double
psc_pulse_none_field(struct psc_pulse *pulse,
		     double xx, double yy, double zz, double tt)
{
  return 0.;
}

struct psc_pulse_ops psc_pulse_none_ops = {
  .name        = "none",
  .field_p     = psc_pulse_none_field,
  .field_s     = psc_pulse_none_field,
};
