
#include "psc.h"
#include "psc_particles_as_double.h"
#include "psc_fields_as_c.h"

// ======================================================================
// !!! These moments are shifted to (n+.5) * dt, rather than n * dt,
// since the "double" particles are shifted that way.
//
// node-centered moments -- probably mostly useful for testing
// charge continuity

#include "psc_output_fields_item_moments_2nd_nc.cxx"

// ======================================================================
// psc_output_fields_item: subclass "n_2nd_nc_double"

struct psc_output_fields_item_ops_n : psc_output_fields_item_ops {
  psc_output_fields_item_ops_n() {
    name               = "n_2nd_nc_double";
    nr_comp            = 1;
    fld_names[0]       = "n_nc";
    run_all            = n_run_all;
    flags              = POFI_ADD_GHOSTS | POFI_BY_KIND;
  }
} psc_output_fields_item_n_2nd_nc_double_ops;

// ======================================================================
// psc_output_fields_item: subclass "rho_2nd_nc_double"

struct psc_output_fields_item_ops_rho : psc_output_fields_item_ops {
  psc_output_fields_item_ops_rho() {
    name               = "rho_2nd_nc_double";
    nr_comp            = 1;
    fld_names[0]       = "rho_nc";
    run_all            = rho_run_all;
    flags              = POFI_ADD_GHOSTS;
  }
} psc_output_fields_item_rho_2nd_nc_double_ops;

