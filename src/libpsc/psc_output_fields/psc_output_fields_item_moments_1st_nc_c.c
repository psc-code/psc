
#include "psc.h"
#include "psc_particles_as_c.h"
#include "psc_fields_as_c.h"

// ======================================================================
// node-centered moments -- probably mostly useful for testing
// charge continuity

#include "psc_output_fields_item_moments_1st_nc.c"

// ======================================================================
// psc_output_fields_item: subclass "n_1st_nc_c"

struct psc_output_fields_item_ops psc_output_fields_item_n_1st_nc_c_ops = {
  .name               = "n_1st_nc_c",
  .get_component_name = n_get_component_name,
  .get_nr_components  = n_get_nr_components,
  .run                = n_run,
};

