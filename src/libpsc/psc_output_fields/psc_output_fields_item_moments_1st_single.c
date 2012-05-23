
#include "psc.h"
#include "psc_particles_as_single.h"
#include "psc_fields_as_c.h"

// ======================================================================
// !!! These moments are shifted to (n+.5) * dt, rather than n * dt,
// since the "single" particles are shifted that way.

#include "psc_output_fields_item_moments_1st.c"

// ======================================================================
// psc_output_fields_item: subclass "n_1st_single"

struct psc_output_fields_item_ops psc_output_fields_item_n_1st_single_ops = {
  .name               = "n_1st_single",
  .get_component_name = n_1st_get_component_name,
  .get_nr_components  = n_1st_get_nr_components,
  .run                = n_1st_run,
  .flags              = POFI_ADD_GHOSTS,
};

// ======================================================================
// psc_output_fields_item: subclass "v_1st_single"

struct psc_output_fields_item_ops psc_output_fields_item_v_1st_single_ops = {
  .name               = "v_1st_single",
  .get_component_name = v_1st_get_component_name,
  .get_nr_components  = v_1st_get_nr_components,
  .run                = v_1st_run,
  .flags              = POFI_ADD_GHOSTS,
};

// ======================================================================
// psc_output_fields_item: subclass "vv_1st_single"

struct psc_output_fields_item_ops psc_output_fields_item_vv_1st_single_ops = {
  .name               = "vv_1st_single",
  .get_component_name = vv_1st_get_component_name,
  .get_nr_components  = vv_1st_get_nr_components,
  .run                = vv_1st_run,
  .flags              = POFI_ADD_GHOSTS,
};

