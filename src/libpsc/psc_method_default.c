
#include "psc_method_private.h"

#include <psc_push_particles.h>
#include <psc_diag.h>
#include <psc_output_fields_collection.h>
#include <psc_output_particles.h>

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
// psc_method_default_output

void
psc_method_default_output(struct psc_method *method, struct psc *psc)
{
  psc_diag_run(psc->diag, psc);
  psc_output_fields_collection_run(psc->output_fields_collection, psc->flds, psc->particles);
  psc_output_particles_run(psc->output_particles, psc->particles);
}

// ----------------------------------------------------------------------
// psc_method "default"

struct psc_method_ops psc_method_ops_default = {
  .name                = "default",
  .initialize          = psc_method_default_initialize,
};
