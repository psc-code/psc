
#include "cuda_mfields.h"
#include "cuda_bits.h"

#include <cstdio>
#include <cassert>

// ----------------------------------------------------------------------
// cuda_mfields_create

struct cuda_mfields *
cuda_mfields_create()
{
  struct cuda_mfields *cmflds = 
    (struct cuda_mfields *) calloc(1, sizeof(*cmflds));

  return cmflds;
}

// ----------------------------------------------------------------------
// cuda_mfields_destroy

void
cuda_mfields_destroy(struct cuda_mfields *cmflds)
{
  free(cmflds);
}

// ----------------------------------------------------------------------
// cuda_mfields_ctor

void
cuda_mfields_ctor(struct cuda_mfields *cmflds, mrc_json_t json)
{
  cudaError_t ierr;

  mrc_json_t json_info = mrc_json_get_object_entry(json, "info");
  
  cmflds->n_patches = mrc_json_get_object_entry_integer(json_info, "n_patches");
  cmflds->n_fields = mrc_json_get_object_entry_integer(json_info, "n_fields");
  mrc_json_get_object_entry_int3(json_info, "ib", cmflds->ib);
  mrc_json_get_object_entry_int3(json_info, "im", cmflds->im);
  mrc_json_get_object_entry_int3(json_info, "ldims", cmflds->ldims);
  double dx[3];
  mrc_json_get_object_entry_double3(json_info, "dx", dx);
  for (int d = 0; d < 3; d++) {
    cmflds->dx[d] = dx[d];
  }

  cmflds->n_cells_per_patch = cmflds->im[0] * cmflds->im[1] * cmflds->im[2];
  cmflds->n_cells = cmflds->n_patches * cmflds->n_cells_per_patch;

  ierr = cudaMalloc((void **) &cmflds->d_flds,
		    cmflds->n_fields * cmflds->n_cells * sizeof(*cmflds->d_flds)); cudaCheck(ierr);

  cmflds->d_flds_by_patch = new fields_cuda_real_t *[cmflds->n_patches];
  for (int p = 0; p < cmflds->n_patches; p++) {
    cmflds->d_flds_by_patch[p] = cmflds->d_flds + p * cmflds->n_fields * cmflds->n_cells_per_patch;
  }
}

// ----------------------------------------------------------------------
// cuda_mfields_dtor

void
cuda_mfields_dtor(struct cuda_mfields *cmflds)
{
  cudaError_t ierr;

  ierr = cudaFree(cmflds->d_flds); cudaCheck(ierr);
  
  delete[] cmflds->d_flds_by_patch;
}

// ----------------------------------------------------------------------
// cuda_mfields_get_host_fields

fields_single_t
cuda_mfields_get_host_fields(struct cuda_mfields *cmflds)
{
  return fields_single_t_ctor(cmflds->ib, cmflds->im, cmflds->n_fields);
}

// ----------------------------------------------------------------------
// cuda_mfields_copy_to_device

void
cuda_mfields_copy_to_device(struct cuda_mfields *cmflds, int p, fields_single_t h_flds, int mb, int me)
{
  cudaError_t ierr;
  
  if (mb == me) {
    return;
  }
  assert(mb < me);

  unsigned int size = cmflds->n_cells_per_patch;
  ierr = cudaMemcpy(cmflds->d_flds_by_patch[p] + mb * size,
		    h_flds.data + mb * size,
		    (me - mb) * size * sizeof(float),
		    cudaMemcpyHostToDevice); cudaCheck(ierr);
}

// ----------------------------------------------------------------------
// cuda_mfields_copy_from_device

void
cuda_mfields_copy_from_device(struct cuda_mfields *cmflds, int p, fields_single_t h_flds, int mb, int me)
{
  cudaError_t ierr;

  if (mb == me) {
    return;
  }
  assert(mb < me);

  unsigned int size = cmflds->n_cells_per_patch;
  ierr = cudaMemcpy(h_flds.data + mb * size,
		    cmflds->d_flds_by_patch[p] + mb * size,
		    (me - mb) * size * sizeof(float),
		    cudaMemcpyDeviceToHost); cudaCheck(ierr);
}

