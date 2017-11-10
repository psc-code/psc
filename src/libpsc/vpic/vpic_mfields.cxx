
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
// vpic_mfields_ctor_from_simulation

void
vpic_mfields_ctor_from_simulation(struct vpic_mfields *vmflds)
{
  vmflds->field_array = simulation->field_array;

  // Accessing the data as a C array relies on fields_t to not change
  assert(sizeof(vmflds->field_array->f[0]) / sizeof(float) == VPIC_MFIELDS_N_COMP);
}

// ----------------------------------------------------------------------
// vpic_mfields_get_data

float *vpic_mfields_get_data(struct vpic_mfields *vmflds, int *ib, int *im)
{
  if (!vmflds->field_array) {
    MHERE;
    return NULL;
  }
  const int B = 1; // VPIC always uses one ghost cell (on c.c. grid)
  grid_t *g = vmflds->field_array->g;
  im[0] = g->nx + 2*B;
  im[1] = g->ny + 2*B;
  im[2] = g->nz + 2*B;
  ib[0] = -B;
  ib[1] = -B;
  ib[2] = -B;
  
  return &vmflds->field_array->f[0].ex;
}

