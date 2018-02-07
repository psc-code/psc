
#ifndef CUDA_MFIELDS_CONST_H
#define CUDA_MFIELDS_CONST_H

#include <cuda_bits.h>

#include <cstdio>

// ----------------------------------------------------------------------
// cuda_mfields_const
//
// cuda_mfields parameters in CUDA constant memory
// (to avoid having to pass them in)

struct cuda_mfields_const {
  int ib[3];
  int im[3];
  int n_cells_per_patch;
};

__constant__ __device__ struct cuda_mfields_const d_cmflds_const;

static void
cuda_mfields_const_set(struct cuda_mfields *cmflds)
{
  struct cuda_mfields_const c;
  for (int d = 0; d < 3; d++) {
    c.ib[d] = cmflds->ib[d];
    c.im[d] = cmflds->im[d];
  }
  c.n_cells_per_patch = cmflds->n_cells_per_patch;
  /* assert(c.im[0] == 1); */
  /* assert(c.ib[0] == 0 && c.ib[1] == -2 && c.ib[2] == -2); */

  cudaError_t ierr = cudaMemcpyToSymbol(d_cmflds_const, &c, sizeof(c)); cudaCheck(ierr);
}

#endif

