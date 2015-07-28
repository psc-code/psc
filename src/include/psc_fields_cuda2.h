
#ifndef PSC_FIELDS_CUDA2_H
#define PSC_FIELDS_CUDA2_H

#include "psc_fields_private.h"

typedef float fields_cuda2_real_t;

#define MPI_FIELDS_CUDA2_REAL MPI_FLOAT

struct psc_fields_cuda2 {
  fields_cuda2_real_t *d_flds;
};

#define psc_fields_cuda2(pf) mrc_to_subobj(pf, struct psc_fields_cuda2)

// ----------------------------------------------------------------------

struct psc_mfields_cuda2 {
  fields_cuda2_real_t *d_flds;
};

#define psc_mfields_cuda2(pf) mrc_to_subobj(pf, struct psc_mfields_cuda2)

#endif
