
#undef _GLIBCXX_USE_INT128

#include "cuda_mparticles.h"

#include <thrust/functional.h>
#include <thrust/transform_scan.h>
#include <thrust/count.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

#include "psc_cuda.h"
#include "particles_cuda.h"

#include <b40c/radixsort_scanscatter_kernel4.h>

using namespace b40c_thrust;

typedef unsigned int K;
typedef unsigned int V;

// static const int RADIX_BITS = 4;

struct count_if_equal : public thrust::unary_function<unsigned int, unsigned int> {
  const unsigned int value;

  __device__ __host__ count_if_equal(unsigned int _value) : value(_value) { }

  __device__ __host__ unsigned int operator()(unsigned int value_in) {
    return value_in == value;
  }
};

#if 0

EXTERN_C int
cuda_exclusive_scan_2(struct psc_particles *prts, unsigned int *_d_vals,
		      unsigned int *_d_sums, int n_prts)
{
  struct cuda_mparticles *cmprts = psc_mparticles_cuda(prts->mprts)->cmprts;
  thrust::device_ptr<unsigned int> d_vals(_d_vals);
  thrust::device_ptr<unsigned int> d_sums(_d_sums);
  
  count_if_equal unary_op(cmprts->n_blocks_per_patch);
  thrust::transform_exclusive_scan(d_vals, d_vals + n_prts, d_sums, unary_op,
				   0, thrust::plus<unsigned int>());

  // OPT, don't mv to host
  int sum = d_sums[n_prts - 1] + (d_vals[n_prts - 1] == cmprts->n_blocks_per_patch);
  return sum;
}

EXTERN_C int
_cuda_exclusive_scan_2(struct psc_particles *prts, unsigned int *d_bidx,
		       unsigned int *d_sums, int n_prts)
{
  struct cuda_mparticles *cmprts = psc_mparticles_cuda(prts->mprts)->cmprts;
  unsigned int *bidx = new unsigned int[n_prts];
  unsigned int *sums = new unsigned int[n_prts];
  check(cudaMemcpy(bidx, d_bidx, n_prts * sizeof(*bidx),
		   cudaMemcpyDeviceToHost));

  unsigned int sum = 0;
  for (int i = 0; i < n_prts; i++) {
    sums[i] = sum;
    sum += (bidx[i] == cmprts->n_blocks_per_patch ? 1 : 0);
  }

  check(cudaMemcpy(d_sums, sums, n_prts * sizeof(*d_sums),
		   cudaMemcpyHostToDevice));
  delete[] sums;
  delete[] bidx;
  return sum;
}

#endif

