
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
cuda_mfields_ctor(struct cuda_mfields *cmflds, int ib[3], int im[3],
		  int n_fields, int n_patches)
{
  cudaError_t ierr;
  
  cmflds->n_patches = n_patches;
  cmflds->n_fields = n_fields;
  
  for (int d = 0; d < 3; d++) {
    cmflds->im[d] = im[d];
    cmflds->ib[d] = ib[d];
  }

  cmflds->n_cells_per_patch = im[0] * im[1] * im[2];
  cmflds->n_cells = n_patches * cmflds->n_cells_per_patch;

  ierr = cudaMalloc((void **) &cmflds->d_flds,
		    n_fields * cmflds->n_cells * sizeof(*cmflds->d_flds)); cudaCheck(ierr);

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

