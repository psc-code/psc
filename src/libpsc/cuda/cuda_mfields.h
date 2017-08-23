
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
// cuda_mfields_bnd_patch

struct cuda_mfields_bnd_patch {
  fields_cuda_real_t *arr_off;
  int im[3];
  int ib[3];
  fields_cuda_real_t *arr;
};

// ----------------------------------------------------------------------
// cuda_mfields

struct cuda_mfields {
  fields_cuda_real_t *d_flds;
  int ib[3], im[3]; // FIXME, should be called off, ldims
  int n_patches;
  int n_fields;
  int n_cells_per_patch;
  int n_cells;
  fields_cuda_real_t **d_flds_by_patch;

  // bnd related
  struct cuda_mfields_bnd_patch *bnd_by_patch;
  fields_cuda_real_t *d_bnd_buf;
  fields_cuda_real_t *h_bnd_buf;
  int *h_nei_patch;
  int *d_nei_patch;
};

#endif
