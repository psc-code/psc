
#include "gtest/gtest.h"

#include "cuda_mfields.h"
#include "cuda_bits.h"

// ----------------------------------------------------------------------
// cuda_mfields_const
//
// cuda_mfields parameters in CUDA constant memory
// (to avoid having to pass them in)

struct _cuda_mfields_const {
  int ib[3];
  int im[3];
  int n_cells_per_patch;
};

__constant__ __device__ struct _cuda_mfields_const _d_cmflds_const;

void _cuda_mfields_const_set(struct cuda_mfields *cmflds)
{
  struct _cuda_mfields_const c;
  for (int d = 0; d < 3; d++) {
    c.ib[d] = cmflds->ib[d];
    c.im[d] = cmflds->im[d];
  }
  c.n_cells_per_patch = cmflds->n_cells_per_patch;
  /* assert(c.im[0] == 1); */
  /* assert(c.ib[0] == 0 && c.ib[1] == -2 && c.ib[2] == -2); */

  cudaError_t ierr = cudaMemcpyToSymbol(_d_cmflds_const, &c, sizeof(c)); cudaCheck(ierr);
}

// ======================================================================
// _DFields

struct _DFields
{
  using real_t = float;
  
  __host__ __device__ _DFields(real_t* d_flds)
    : d_flds_(d_flds)
  {}
  
  __device__ real_t  operator()(int m, int i, int j, int k) const { return d_flds_[index(m, i,j,k)]; }
  __device__ real_t& operator()(int m, int i, int j, int k)       { return d_flds_[index(m, i,j,k)]; }

  __host__ real_t *d_flds() { return d_flds_; }

private:
  __device__ int index(int m, int i, int j, int k) const
  {
    return ((((m)
	      *_d_cmflds_const.im[2] + (k + 2))
	     *_d_cmflds_const.im[1] + (j + 2))
	    *1 + (0));
  }

private:
  real_t *d_flds_;
};

// ======================================================================
// _DMFields

struct _DMFields
{
  using real_t = float;
  
  __host__ _DMFields(cuda_mfields *cmflds)
    : d_flds_(cmflds->data()),
      stride_(cmflds->n_cells_per_patch * cmflds->n_fields)
  {}
  
  __host__ __device__ _DFields operator[](int p) { return _DFields(d_flds_ + p * stride_); }

private:
  real_t *d_flds_;
  uint stride_;
};

// ======================================================================
// test1

__global__ void test1(_DMFields MF)
{
  int j = threadIdx.x + blockDim.x * blockIdx.x;
  int k = threadIdx.y + blockDim.y * blockIdx.y;
  int p = threadIdx.z;

  auto F = MF[p];

  F(1, 0,j,k) += 5.;
}

// ======================================================================
// test2

__global__ void test2(DMFields MF)
{
  int j = threadIdx.x + blockDim.x * blockIdx.x;
  int k = threadIdx.y + blockDim.y * blockIdx.y;
  int p = threadIdx.z;

  auto F = MF[p];

  F(1, 0,j,k) += 5.;
}


// ======================================================================
// Accel test

TEST(Sample, SampleCase)
{
}

