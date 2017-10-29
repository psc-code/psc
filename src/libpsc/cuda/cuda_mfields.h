
#ifndef CUDA_MFIELDS_H
#define CUDA_MFIELDS_H

#include "cuda_iface.h"
#include "cuda_iface_bnd.h"

// ----------------------------------------------------------------------
// cuda_mfields_bnd

struct cuda_mfields_bnd {
  struct cuda_mfields_bnd_patch *bnd_by_patch;
  fields_cuda_real_t *d_bnd_buf;
  fields_cuda_real_t *h_bnd_buf;
  int *h_nei_patch;
  int *d_nei_patch;
  int *h_map_out; // maps thread id to a particular offset for ghosts in the flds array 
  int *d_map_out;
  int nr_map_out; // number of entries in the map
  int *h_map_in; // maps thread id to a particular offset for ghosts in the flds array 
  int *d_map_in;
  int nr_map_in; // number of entries in the map
  int n_patches;
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
  int ldims[3];                   // number of cells per direction in each patch
  float dx[3];                    // cell size (in actual length units)
};

#endif
