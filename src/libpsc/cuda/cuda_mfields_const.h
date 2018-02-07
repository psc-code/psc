
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

// ======================================================================
// DFields

struct DFields
{
  using real_t = float;
  
  __host__ __device__ DFields(real_t* d_flds, int im[3])
    : d_flds_(d_flds),
      im_{ im[0], im[1], im[2] }
  {}
  
  __device__ real_t  operator()(int m, int i, int j, int k) const { return d_flds_[index(m, i,j,k)]; }
  __device__ real_t& operator()(int m, int i, int j, int k)       { return d_flds_[index(m, i,j,k)]; }

  __host__ real_t *d_flds() { return d_flds_; }

private:
  __device__ int index(int m, int i, int j, int k) const
  {
    return ((((m)
	      *im_[2] + (k + 2))
	     *im_[1] + (j + 2))
	    *1 + (0));
  }

private:
  real_t *d_flds_;
  int im_[3];
};

// ======================================================================
// DMFields

struct DMFields
{
  using real_t = float;
  
  __host__ DMFields(cuda_mfields *cmflds)
    : d_flds_(cmflds->d_flds.data().get()),
      stride_(cmflds->n_cells_per_patch * cmflds->n_fields),
      im_{ cmflds->im[0], cmflds->im[1], cmflds->im[2] }
  {}
  
  __host__ __device__ DFields operator[](int p)
  {
    return DFields(d_flds_ + p * stride_, im_); }

private:
  real_t *d_flds_;
  uint stride_;
  int im_[3];
};

#endif

