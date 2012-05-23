
#include "psc.h"
#include "psc_particles_as_single.h"
#include "psc_fields_as_c.h"

// ======================================================================
// !!! These moments are shifted to (n+.5) * dt, rather than n * dt,
// since the "single" particles are shifted that way.
//
// node-centered moments -- probably mostly useful for testing
// charge continuity

#include "psc_output_fields_item_moments_1st_nc.c"

// ======================================================================
// psc_output_fields_item: subclass "n_1st_nc_single"

struct psc_output_fields_item_ops psc_output_fields_item_n_1st_nc_single_ops = {
  .name               = "n_1st_nc_single",
  .get_component_name = n_get_component_name,
  .get_nr_components  = n_get_nr_components,
  .run                = n_run,
  .flags              = POFI_ADD_GHOSTS,
};

