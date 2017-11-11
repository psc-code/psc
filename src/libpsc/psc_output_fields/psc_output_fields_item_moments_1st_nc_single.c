
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
  .nr_comp            = 1,
  .fld_names          = { "n_nc" },
  .run_all             = n_run_all,
  .flags              = POFI_ADD_GHOSTS | POFI_BY_KIND,
};

// ======================================================================
// psc_output_fields_item: subclass "rho_1st_nc_single"

struct psc_output_fields_item_ops psc_output_fields_item_rho_1st_nc_single_ops = {
  .name               = "rho_1st_nc_single",
  .nr_comp            = 1,
  .fld_names          = { "rho_nc" },
  .run_all            = rho_run_all,
  .flags              = POFI_ADD_GHOSTS,
};

// ======================================================================
// psc_output_fields_item: subclass "v_1st_nc_single"

struct psc_output_fields_item_ops psc_output_fields_item_v_1st_nc_single_ops = {
  .name               = "v_1st_nc_single",
  .nr_comp            = 3,
  .fld_names          = { "vx_nc", "vy_nc", "vz_nc" },
  .run_all             = v_run_all,
  .flags              = POFI_ADD_GHOSTS | POFI_BY_KIND,
};

