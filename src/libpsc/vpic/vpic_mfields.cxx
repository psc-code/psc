
#include "simulation.h"
#include "vpic_mparticles.h"

#include <cassert>

// ======================================================================
// vpic_mfields

// ----------------------------------------------------------------------
// vpic_mfields_new_hydro_array

HydroArray* vpic_mfields_new_hydro_array(Simulation *sim)
{
  HydroArray* vmflds = static_cast<HydroArray*>(sim->hydro_array_);

  // Accessing the data as a C array relies on hydro_array_t to not change
  assert(sizeof(vmflds->h[0]) / sizeof(float) == VPIC_HYDRO_N_COMP);
  return vmflds;
}

// ----------------------------------------------------------------------
// vpic_mfields_hydro_get_data

float *vpic_mfields_hydro_get_data(struct HydroArray *vmflds, int *ib, int *im)
{
  return vmflds->getData(ib, im);
}

// ----------------------------------------------------------------------
// vpic_mfields_new_fields_array

FieldArray* vpic_mfields_new_fields_array(Simulation *sim)
{
  FieldArray* vmflds = static_cast<FieldArray*>(sim->field_array_);

  // Accessing the data as a C array relies on fields_array_t to not change
  assert(sizeof(vmflds->f[0]) / sizeof(float) == VPIC_MFIELDS_N_COMP);
  return vmflds;
}

// ----------------------------------------------------------------------
// vpic_mfields_get_data

float *vpic_mfields_get_data(struct FieldArray *vmflds, int *ib, int *im)
{
  return vmflds->getData(ib, im);
}

// ----------------------------------------------------------------------
// C wrappers

double vpic_mfields_synchronize_tang_e_norm_b(struct FieldArray *vmflds)
{
  return vmflds->synchronize_tang_e_norm_b();
}

void vpic_mfields_compute_div_b_err(struct FieldArray *vmflds)
{
  vmflds->compute_div_b_err();
}

double vpic_mfields_compute_rms_div_b_err(struct FieldArray *vmflds)
{
  return vmflds->compute_rms_div_b_err();
}

void vpic_mfields_clean_div_b(struct FieldArray *vmflds)
{
  vmflds->clean_div_b();
}

void vpic_mfields_compute_div_e_err(struct FieldArray *vmflds)
{
  vmflds->compute_div_e_err();
}

double vpic_mfields_compute_rms_div_e_err(struct FieldArray *vmflds)
{
  return vmflds->compute_rms_div_e_err();
}

void vpic_mfields_clean_div_e(struct FieldArray *vmflds)
{
  vmflds->clean_div_e();
}

void vpic_mfields_clear_rhof(struct FieldArray *vmflds)
{
  vmflds->clear_rhof();
}

void vpic_mfields_synchronize_rho(struct FieldArray *vmflds)
{
  vmflds->synchronize_rho();
}

void vpic_mfields_compute_rhob(struct FieldArray *vmflds)
{
  vmflds->compute_rhob();
}

void vpic_mfields_compute_curl_b(struct FieldArray *vmflds)
{
  vmflds->compute_curl_b();
}

void vpic_mfields_accumulate_rho_p(struct FieldArray *vmflds,
				   struct vpic_mparticles *vmprts)
{
  species_t *sp;
  LIST_FOR_EACH(sp, vmprts->sl_)
    TIC accumulate_rho_p(vmflds, sp); TOC( accumulate_rho_p, 1);
}

