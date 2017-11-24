
#include "vpic_iface.h"

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

float *vpic_mfields_hydro_get_data(HydroArray *vmflds, int *ib, int *im)
{
  return vmflds->getData(ib, im);
}

// ----------------------------------------------------------------------
// vpic_mfields_new_fields_array

FieldArray* vpic_mfields_new_fields_array(Simulation *sim)
{
  FieldArray* vmflds = sim->field_array_;

  // Accessing the data as a C array relies on fields_array_t to not change
  assert(sizeof(vmflds->f[0]) / sizeof(float) == VPIC_MFIELDS_N_COMP);
  return vmflds;
}

// ----------------------------------------------------------------------
// vpic_mfields_get_data

float *vpic_mfields_get_data(FieldArray *vmflds, int *ib, int *im)
{
  return vmflds->getData(ib, im);
}

// ----------------------------------------------------------------------
// C wrappers

double vpic_mfields_synchronize_tang_e_norm_b(FieldArray *vmflds)
{
  return vmflds->synchronize_tang_e_norm_b();
}

void vpic_mfields_compute_div_b_err(Simulation* sim, FieldArray* vmflds)
{
  TIC sim->compute_div_b_err(*vmflds); TOC(compute_div_b_err, 1);
}

double vpic_mfields_compute_rms_div_b_err(FieldArray *vmflds)
{
  return vmflds->compute_rms_div_b_err();
}

void vpic_mfields_clean_div_b(FieldArray *vmflds)
{
  vmflds->clean_div_b();
}

void vpic_mfields_compute_div_e_err(FieldArray *vmflds)
{
  vmflds->compute_div_e_err();
}

double vpic_mfields_compute_rms_div_e_err(FieldArray *vmflds)
{
  return vmflds->compute_rms_div_e_err();
}

void vpic_mfields_clean_div_e(FieldArray *vmflds)
{
  vmflds->clean_div_e();
}

void vpic_mfields_clear_rhof(Simulation* sim, FieldArray* vmflds)
{
  TIC sim->clear_rhof(*vmflds); TOC(clear_jf, 1);
}

void vpic_mfields_synchronize_rho(FieldArray *vmflds)
{
  vmflds->synchronize_rho();
}

void vpic_mfields_compute_rhob(FieldArray *vmflds)
{
  vmflds->compute_rhob();
}

void vpic_mfields_compute_curl_b(FieldArray *vmflds)
{
  vmflds->compute_curl_b();
}

void vpic_mfields_accumulate_rho_p(FieldArray *vmflds, Particles *vmprts)
{
  for (Particles::Iter sp = vmprts->begin(); sp != vmprts->end(); ++sp) {
    TIC accumulate_rho_p(vmflds, &*sp); TOC( accumulate_rho_p, 1);
  }
}

