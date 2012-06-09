
#include "psc_cuda.h"
#include "particles_cuda.h"

#include <thrust/functional.h>
#include <thrust/transform_scan.h>
#include <thrust/count.h>

#define NBLOCKS_X 1
#define NBLOCKS_Y 8
#define NBLOCKS_Z 8

struct count_if_equal : public thrust::unary_function<unsigned int, unsigned int> {
  const unsigned int value;

  __device__ __host__ count_if_equal(unsigned int _value) : value(_value) { }

  __device__ __host__ unsigned int operator()(unsigned int value_in) {
    return value_in == value;
  }
};

EXTERN_C int
cuda_exclusive_scan_2(struct psc_particles *prts, unsigned int *_d_vals,
		      unsigned int *_d_sums)
{
  struct psc_particles_cuda *cuda = psc_particles_cuda(prts);
  thrust::device_ptr<unsigned int> d_vals(_d_vals);
  thrust::device_ptr<unsigned int> d_sums(_d_sums);

  count_if_equal unary_op(cuda->nr_blocks);
  thrust::transform_exclusive_scan(d_vals, d_vals + prts->n_part, d_sums, unary_op,
				   0, thrust::plus<unsigned int>());

  // OPT, don't mv to host
  int sum = d_sums[prts->n_part - 1] + (d_vals[prts->n_part - 1] == cuda->nr_blocks);
  return sum;
}

EXTERN_C int
_cuda_exclusive_scan_2(struct psc_particles *prts, unsigned int *d_bidx,
		       unsigned int *d_sums)
{
  struct psc_particles_cuda *cuda = psc_particles_cuda(prts);
  unsigned int *bidx = new unsigned int[prts->n_part];
  unsigned int *sums = new unsigned int[prts->n_part];
  check(cudaMemcpy(bidx, d_bidx, prts->n_part * sizeof(*bidx),
		   cudaMemcpyDeviceToHost));

  unsigned int sum = 0;
  for (int i = 0; i < prts->n_part; i++) {
    sums[i] = sum;
    sum += (bidx[i] == cuda->nr_blocks ? 1 : 0);
  }

  check(cudaMemcpy(d_sums, sums, prts->n_part * sizeof(*d_sums),
		   cudaMemcpyHostToDevice));
  delete[] sums;
  delete[] bidx;
  return sum;
}

#if 0
void
cuda_mprts_scan_send_buf(struct psc_mparticles *mprts)
{
  for (int p = 0; p < mprts->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
    struct psc_particles_cuda *cuda = psc_particles_cuda(prts);
    cuda->bnd_n_send = cuda_exclusive_scan_2(prts, cuda->h_dev->bidx, cuda->h_dev->sums);
  }
}
#endif

void
cuda_mprts_find_n_send(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);

  cuda_mprts_spine_reduce(mprts);

  thrust::device_ptr<unsigned int> d_spine_cnts(mprts_cuda->d_bnd_spine_cnts);
  thrust::device_ptr<unsigned int> d_spine_sums(mprts_cuda->d_bnd_spine_sums);
  unsigned int off = 0;
  for (int p = 0; p < mprts->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
    struct psc_particles_cuda *cuda = psc_particles_cuda(prts);
    int nr_total_blocks = cuda->nr_blocks * mprts->nr_patches;

    unsigned int n_send = d_spine_sums[nr_total_blocks * 10 + (p + 1) * cuda->nr_blocks];
    cuda->bnd_n_send = n_send - off;
    off = n_send;
  }
  mprts_cuda->nr_prts_send = off;
}

void
cuda_mprts_scan_send_buf_total(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);
  int nr_blocks = 0;

  int nr_send = 0;
  for (int p = 0; p < mprts->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
    struct psc_particles_cuda *cuda = psc_particles_cuda(prts);
    nr_send += cuda->bnd_n_send;
    nr_blocks = cuda->nr_blocks;
  }
  assert(nr_blocks > 0);
  
  thrust::device_ptr<unsigned int> d_vals(mprts_cuda->d_bidx);
  thrust::device_ptr<unsigned int> d_sums(mprts_cuda->d_sums);

  count_if_equal unary_op(mprts->nr_patches * nr_blocks);
  thrust::transform_exclusive_scan(d_vals, d_vals + mprts_cuda->nr_prts, d_sums, unary_op,
				   0, thrust::plus<unsigned int>());
  
  // OPT, don't mv to host
  int sum = d_sums[mprts_cuda->nr_prts - 1] + 
    (d_vals[mprts_cuda->nr_prts - 1] == mprts->nr_patches * nr_blocks);
  assert(sum == nr_send);
}
