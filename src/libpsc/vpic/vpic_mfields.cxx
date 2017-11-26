
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

void vpic_mfields_accumulate_rho_p(FieldArray *vmflds, Particles *vmprts)
{
  for (Particles::Iter sp = vmprts->begin(); sp != vmprts->end(); ++sp) {
    TIC accumulate_rho_p(vmflds, &*sp); TOC( accumulate_rho_p, 1);
  }
}

