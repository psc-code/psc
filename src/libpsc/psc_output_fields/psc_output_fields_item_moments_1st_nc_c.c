
#include "psc.h"
#include "psc_particles_as_c.h"
#include "psc_fields_as_c.h"

// ======================================================================
// node-centered moments -- probably mostly useful for testing
// charge continuity
//
// NOTE: This is at time t^n, which may not be what you want for checking
// continuity / Gauss's Law, since that should hold at t^{n+.5}

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

// ======================================================================
// psc_output_fields_item: subclass "rho_1st_nc_c"

struct psc_output_fields_item_ops psc_output_fields_item_rho_1st_nc_c_ops = {
  .name               = "rho_1st_nc_c",
  .nr_comp            = 1,
  .fld_names          = { "rho_nc" },
  .run_patches        = rho_run_patches,
  .flags              = POFI_ADD_GHOSTS,
};

