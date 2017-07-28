
#include "psc.h"
#include "psc_particles_as_double.h"
#include "psc_fields_as_c.h"

// ======================================================================
// !!! These moments are shifted to (n+.5) * dt, rather than n * dt,
// since the "double" particles are shifted that way.
//
// node-centered moments -- probably mostly useful for testing
// charge continuity

#include "psc_output_fields_item_moments_2nd_nc.c"

// ======================================================================
// psc_output_fields_item: subclass "n_2nd_nc_double"

struct psc_output_fields_item_ops psc_output_fields_item_n_2nd_nc_double_ops = {
  .name               = "n_2nd_nc_double",
  .nr_comp            = 1,
  .fld_names          = { "n_nc" },
  .run_all            = n_run_all,
  .flags              = POFI_ADD_GHOSTS | POFI_BY_KIND,
};

// ======================================================================
// psc_output_fields_item: subclass "rho_2nd_nc_double"

struct psc_output_fields_item_ops psc_output_fields_item_rho_2nd_nc_double_ops = {
  .name               = "rho_2nd_nc_double",
  .nr_comp            = 1,
  .fld_names          = { "rho_nc" },
  .run_all            = rho_run_all,
  .flags              = POFI_ADD_GHOSTS,
};

