
#ifndef PSC_FIELDS_CUDA_H
#define PSC_FIELDS_CUDA_H

#include "psc_fields_private.h"

#define FTYPE FTYPE_CUDA
#include "psc_fields_common.h"
#undef FTYPE

// ----------------------------------------------------------------------

struct psc_mfields_cuda {
  struct cuda_mfields *cmflds;
  struct cuda_mfields_bnd *cbnd;
};

#define psc_mfields_cuda(pf) mrc_to_subobj(pf, struct psc_mfields_cuda)

#endif
