
#include "psc_method_private.h"

#include <psc_push_particles.h>

// ======================================================================
// psc_method "default"

// ----------------------------------------------------------------------
// psc_method_default_initialize

static void
psc_method_default_initialize(struct psc_method *method, struct psc *psc)
{
  psc_push_particles_stagger(psc->push_particles, psc->particles, psc->flds);

  // initial output / stats
  psc_output(psc);
  psc_stats_log(psc);
  psc_print_profiling(psc);
}

// ----------------------------------------------------------------------
// psc_method "default"

struct psc_method_ops psc_method_ops_default = {
  .name                = "default",
  .initialize          = psc_method_default_initialize,
};
