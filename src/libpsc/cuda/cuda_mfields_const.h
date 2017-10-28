
#ifndef CUDA_MFIELDS_CONST_H
#define CUDA_MFIELDS_CONST_H

// ----------------------------------------------------------------------
// cuda_mfields_const
//
// cuda_mfields parameters in CUDA constant memory
// (to avoid having to pass them in)

__constant__ __device__ struct cuda_mfields_const d_cmflds_const;

static void
cuda_mfields_const_set(struct cuda_mfields *cmflds)
{
  struct cuda_mfields_const c;
  for (int d = 0; d < 3; d++) {
    c.ib[d] = cmflds->ib[d];
    c.im[d] = cmflds->im[d];
  }

  // FIXME, should also use cudaError_t
  int ierr = cudaMemcpyToSymbol(d_cmflds_const, &c, sizeof(c)); assert(ierr == 0); //cudaCheck(ierr);
}

#endif

