
#ifndef CUDA_MFIELDS_H
#define CUDA_MFIELDS_H

#include "cuda_iface.h"
#include "cuda_iface_bnd.h"

#define MAX_BND_FIELDS (17)
#define MAX_BND_COMPONENTS (3)

// ----------------------------------------------------------------------
// cuda_mfields_bnd_map
//
// the maps differs depending on the number of components in the field
// (incorporates specific strides)

struct cuda_mfields_bnd_map {
  int *h_map_out; // maps thread id to a particular offset for ghosts in the flds array 
  int *d_map_out;
  int nr_map_out; // number of entries in the map
  int *h_map_in; // maps thread id to a particular offset for ghosts in the flds array 
  int *d_map_in;
  int nr_map_in; // number of entries in the map
};

// ----------------------------------------------------------------------
// cuda_mfields_bnd

struct cuda_mfields_bnd {
  int n_patches;
  int im[3];
  int ib[3];
  struct cuda_mfields_bnd_patch *bnd_by_patch;
  fields_cuda_real_t *d_bnd_buf;
  fields_cuda_real_t *h_bnd_buf;
  int *h_nei_patch;
  int *d_nei_patch;
  struct cuda_mfields_bnd_map map[MAX_BND_FIELDS];
};

// ----------------------------------------------------------------------
// cuda_mfields

struct cuda_mfields
{
  cuda_mfields(Grid_t& grid, mrc_json_t json);
  cuda_mfields(const cuda_mfields&) = delete;
  ~cuda_mfields();
  
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
