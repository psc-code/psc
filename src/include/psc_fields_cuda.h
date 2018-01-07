
#ifndef PSC_FIELDS_CUDA_H
#define PSC_FIELDS_CUDA_H

#include "fields3d.hxx"

// ----------------------------------------------------------------------

struct psc_mfields_cuda {
  struct cuda_mfields *cmflds;
};

#define psc_mfields_cuda(pf) mrc_to_subobj(pf, struct psc_mfields_cuda)

#endif
