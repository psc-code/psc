
#include "vpic_mfields.h"

#include <cassert>

extern vpic_simulation *simulation;

// ======================================================================
// vpic_mfields

// ----------------------------------------------------------------------
// vpic_mfields_create

struct vpic_mfields *
vpic_mfields_create()
{
  return new vpic_mfields;
}

// ----------------------------------------------------------------------
// vpic_mfields_ctor_from_simulation_fields

void
vpic_mfields_ctor_from_simulation_fields(struct vpic_mfields *vmflds)
{
  vmflds->field_array = simulation->field_array;

  // Accessing the data as a C array relies on fields_t to not change
  assert(sizeof(vmflds->field_array->f[0]) / sizeof(float) == VPIC_MFIELDS_N_COMP);
}

// ----------------------------------------------------------------------
// vpic_mfields_ctor_from_simulation_hydro

void
vpic_mfields_ctor_from_simulation_hydro(struct vpic_mfields *vmflds)
{
  vmflds->hydro_array = simulation->hydro_array;

  // Accessing the data as a C array relies on fields_t to not change
  assert(sizeof(vmflds->hydro_array->h[0]) / sizeof(float) == VPIC_HYDRO_N_COMP);
}

// ----------------------------------------------------------------------
// vpic_mfields_get_data

float *vpic_mfields_get_data(struct vpic_mfields *vmflds, int *ib, int *im)
{
  const int B = 1; // VPIC always uses one ghost cell (on c.c. grid)

  if (vmflds->field_array) {
    grid_t *g = vmflds->field_array->g;
    im[0] = g->nx + 2*B;
    im[1] = g->ny + 2*B;
    im[2] = g->nz + 2*B;
    ib[0] = -B;
    ib[1] = -B;
    ib[2] = -B;
    return &vmflds->field_array->f[0].ex;
  } else if (vmflds->hydro_array) {
    grid_t *g = vmflds->hydro_array->g;
    im[0] = g->nx + 2*B;
    im[1] = g->ny + 2*B;
    im[2] = g->nz + 2*B;
    ib[0] = -B;
    ib[1] = -B;
    ib[2] = -B;
    return &vmflds->hydro_array->h[0].jx;
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// C wrappers

double vpic_mfields_synchronize_tang_e_norm_b(struct vpic_mfields *vmflds)
{
  return vmflds->synchronize_tang_e_norm_b();
}

void vpic_mfields_compute_div_b_err(struct vpic_mfields *vmflds)
{
  vmflds->compute_div_b_err();
}

double vpic_mfields_compute_rms_div_b_err(struct vpic_mfields *vmflds)
{
  return vmflds->compute_rms_div_b_err();
}

void vpic_mfields_clean_div_b(struct vpic_mfields *vmflds)
{
  vmflds->clean_div_b();
}

void vpic_mfields_compute_div_e_err(struct vpic_mfields *vmflds)
{
  vmflds->compute_div_e_err();
}

double vpic_mfields_compute_rms_div_e_err(struct vpic_mfields *vmflds)
{
  return vmflds->compute_rms_div_e_err();
}

void vpic_mfields_clean_div_e(struct vpic_mfields *vmflds)
{
  vmflds->clean_div_e();
}

void vpic_mfields_clear_rhof(struct vpic_mfields *vmflds)
{
  vmflds->clear_rhof();
}

void vpic_mfields_accumulate_rho_p(struct vpic_mfields *vmflds,
				   struct vpic_mparticles *vmprts)
{
  vmflds->accumulate_rho_p(vmprts);
}

void vpic_mfields_synchronize_rho(struct vpic_mfields *vmflds)
{
  vmflds->synchronize_rho();
}

void vpic_mfields_compute_rhob(struct vpic_mfields *vmflds)
{
  vmflds->compute_rhob();
}

void vpic_mfields_compute_curl_b(struct vpic_mfields *vmflds)
{
  vmflds->compute_curl_b();
}
