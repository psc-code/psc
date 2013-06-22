
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
  .nr_comp            = 1,
  .fld_names          = { "n_nc" },
  .run                = n_run,
  .flags              = POFI_ADD_GHOSTS | POFI_BY_KIND,
};

