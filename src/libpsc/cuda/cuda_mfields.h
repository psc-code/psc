
#ifndef CUDA_MFIELDS_H
#define CUDA_MFIELDS_H

#include "cuda_iface.h"

#ifdef __cplusplus
#define EXTERN_C extern "C"
#else
#define EXTERN_C
#endif

// FIXME, better call it cuda_mfields_real_t
typedef float fields_cuda_real_t;

// ----------------------------------------------------------------------
// cuda_mfields

struct cuda_mfields {
  fields_cuda_real_t *d_flds;
  int ib[3], im[3]; // FIXME, should be called off, ldims
  int n_patches;
  int n_fields;
  int n_cells_per_patch;
  int n_cells;
  struct psc_fields_cuda_bnd **bnd;
};

#endif
