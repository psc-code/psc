
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

// ----------------------------------------------------------------------
// D_F3
//
// field access on the device

#define D_F3_OFF_YZ(m, jy,jz)						\
  ((((m)								\
     *d_cmflds_const.im[2] + ((jz)+2))					\
    *d_cmflds_const.im[1] + ((jy)+2))					\
   *1 + (0))

#define D_F3(d_flds, m, jx,jy,jz) ((d_flds)[D_F3_OFF_YZ(m, jy,jz)])

// ======================================================================
// DFields

struct DFields
{
  using real_t = float;
  
  __device__ DFields(real_t* d_flds)
    : d_flds_(d_flds)
  {}
  
  __device__ real_t  operator()(int m, int i, int j, int k) const { return D_F3(d_flds_, m, i,j,k); }
  __device__ real_t& operator()(int m, int i, int j, int k)       { return D_F3(d_flds_, m, i,j,k); }

  real_t *d_flds_;
  uint stride_;
};

// ======================================================================
// DMFields

struct DMFields
{
  using real_t = float;
  
  __host__ DMFields(cuda_mfields *cmflds)
    : d_flds_(cmflds->d_flds.data().get()),
      stride_(cmflds->n_cells_per_patch * cmflds->n_fields)
  {}
  
  __device__ real_t  operator()(int m, int i, int j, int k, int p) const { return D_F3(d_flds_ + p * stride_, m, i,j,k); }
  __device__ real_t& operator()(int m, int i, int j, int k, int p)       { return D_F3(d_flds_ + p * stride_, m, i,j,k); }

  __device__ DFields operator[](int p) { return DFields(d_flds_ + p * stride_); }
  
  real_t *d_flds_;
  uint stride_;
};

#endif

