
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

static const int RADIX_BITS = 4;

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

void
cuda_mprts_find_n_send(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);

  int nr_total_blocks = mprts_cuda->nr_total_blocks;

  thrust::device_ptr<unsigned int> d_spine_sums(mprts_cuda->d_bnd_spine_sums);
  thrust::host_vector<unsigned int> h_spine_sums(nr_total_blocks + 1);

  thrust::copy(d_spine_sums + nr_total_blocks * 10,
	       d_spine_sums + nr_total_blocks * 11 + 1,
	       h_spine_sums.begin());

  unsigned int off = 0;
  for (int p = 0; p < mprts->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
    struct psc_particles_cuda *cuda = psc_particles_cuda(prts);

    unsigned int n_send = h_spine_sums[(p + 1) * mprts_cuda->nr_blocks];
    cuda->bnd_n_send = n_send - off;
    off = n_send;
  }
  mprts_cuda->nr_prts_send = off;
}

// ======================================================================
// cuda_mprts_reorder_send_by_id

static void __global__
mprts_reorder_send_by_id(unsigned int nr_prts_send, unsigned int *d_xchg_ids,
			 float4 *d_xi4, float4 *d_pxi4,
			 float4 *d_xchg_xi4, float4 *d_xchg_pxi4)
{
  int n = threadIdx.x + THREADS_PER_BLOCK * blockIdx.x;
  if (n >= nr_prts_send) {
    return;
  }

  unsigned int id = d_xchg_ids[n];
  d_xchg_xi4[n]  = d_xi4[id];
  d_xchg_pxi4[n] = d_pxi4[id];
}


void
cuda_mprts_reorder_send_by_id(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);
  struct cuda_mparticles *cmprts = mprts_cuda->cmprts;
  assert(cmprts);
  
  int dimGrid = (mprts_cuda->nr_prts_send + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;

  mprts_reorder_send_by_id<<<dimGrid, THREADS_PER_BLOCK>>>
    (mprts_cuda->nr_prts_send, mprts_cuda->d_ids + mprts_cuda->nr_prts - mprts_cuda->nr_prts_send,
     cmprts->d_xi4, cmprts->d_pxi4,
     cmprts->d_xi4 + mprts_cuda->nr_prts, cmprts->d_pxi4 + mprts_cuda->nr_prts);
}

void
cuda_mprts_reorder_send_by_id_gold(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);
  struct cuda_mparticles *cmprts = mprts_cuda->cmprts;
  assert(cmprts);
  
  thrust::device_ptr<unsigned int> d_ids(mprts_cuda->d_ids);
  thrust::device_ptr<float4> d_xi4(cmprts->d_xi4);
  thrust::device_ptr<float4> d_pxi4(cmprts->d_pxi4);
  thrust::host_vector<unsigned int> h_ids(d_ids, d_ids + mprts_cuda->nr_prts);
  thrust::host_vector<float4> h_xi4(d_xi4, d_xi4 + mprts_cuda->nr_prts + mprts_cuda->nr_prts_send);
  thrust::host_vector<float4> h_pxi4(d_pxi4, d_pxi4 + mprts_cuda->nr_prts + mprts_cuda->nr_prts_send);
  
  for (int n = 0; n < mprts_cuda->nr_prts_send; n++) {
    unsigned int id = h_ids[mprts_cuda->nr_prts - mprts_cuda->nr_prts_send + n];
    h_xi4[mprts_cuda->nr_prts + n]  = h_xi4[id];
    h_pxi4[mprts_cuda->nr_prts + n] = h_pxi4[id];
  }

  thrust::copy(h_xi4.begin(), h_xi4.end(), d_xi4);
  thrust::copy(h_pxi4.begin(), h_pxi4.end(), d_pxi4);
}

// ======================================================================
// cuda_mprts_scan_send_buf_total

void
cuda_mprts_scan_send_buf_total_gold(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);

  unsigned int nr_total_blocks = mprts_cuda->nr_total_blocks;

  thrust::device_ptr<unsigned int> d_bidx(mprts_cuda->d_bidx);
  thrust::device_ptr<unsigned int> d_sums(mprts_cuda->d_sums);
  thrust::device_ptr<unsigned int> d_off(mprts_cuda->d_off);
  thrust::device_ptr<unsigned int> d_spine_sums(mprts_cuda->d_bnd_spine_sums);
  thrust::host_vector<unsigned int> h_off(d_off, d_off + nr_total_blocks + 1);
  thrust::host_vector<unsigned int> h_bidx(d_bidx, d_bidx + mprts_cuda->nr_prts);
  thrust::host_vector<unsigned int> h_sums(d_sums, d_sums + mprts_cuda->nr_prts);
  
  for (unsigned int bid = 0; bid < nr_total_blocks; bid++) {
    unsigned int sum = d_spine_sums[nr_total_blocks * 10 + bid];
    for (int n = h_off[bid]; n < h_off[bid+1]; n++) {
      if (h_bidx[n] == CUDA_BND_S_OOB) {
	h_sums[n] = sum;
	sum++;
      }
    }
  }

  thrust::copy(h_sums.begin(), h_sums.end(), d_sums);

  cuda_mprts_reorder_send_buf_total(mprts);
}

void
cuda_mprts_scan_send_buf_total(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);

  unsigned int nr_total_blocks = mprts_cuda->nr_total_blocks;
  int *b_mx = mprts_cuda->b_mx;

  // OPT, we could do this from the beginning and adapt find_n_send()
  thrust::device_ptr<unsigned int> d_spine_cnts(mprts_cuda->d_bnd_spine_cnts);
  thrust::device_ptr<unsigned int> d_spine_sums(mprts_cuda->d_bnd_spine_sums);
  thrust::exclusive_scan(d_spine_cnts + nr_total_blocks * 10,
			 d_spine_cnts + nr_total_blocks * 11 + 1,
			 d_spine_sums + nr_total_blocks * 10,
			 mprts_cuda->nr_prts - mprts_cuda->nr_prts_send);
  // OPT, we could somehow not fill in ids for not oob at all
  // this should make sure at least those within bounds don't screw anything up
  thrust::fill(d_spine_sums, d_spine_sums + nr_total_blocks * 10, 0);

  if (b_mx[0] == 1 && b_mx[1] == 8 && b_mx[2] == 8) {
    ScanScatterDigits4<K, V, 0, RADIX_BITS, 0,
		       NopFunctor<K>,
		       NopFunctor<K>,
		       8, 8> 
      <<<nr_total_blocks, B40C_RADIXSORT_THREADS>>>
      (mprts_cuda->d_bnd_spine_sums, mprts_cuda->d_bidx,
       mprts_cuda->d_ids, mprts_cuda->d_off, nr_total_blocks);
  } else if (b_mx[0] == 1 && b_mx[1] == 16 && b_mx[2] == 16) {
    ScanScatterDigits4<K, V, 0, RADIX_BITS, 0,
		       NopFunctor<K>,
		       NopFunctor<K>,
		       16, 16> 
      <<<nr_total_blocks, B40C_RADIXSORT_THREADS>>>
      (mprts_cuda->d_bnd_spine_sums, mprts_cuda->d_bidx,
       mprts_cuda->d_ids, mprts_cuda->d_off, nr_total_blocks);
  } else if (b_mx[0] == 1 && b_mx[1] == 32 && b_mx[2] == 32) {
    ScanScatterDigits4<K, V, 0, RADIX_BITS, 0,
		       NopFunctor<K>,
		       NopFunctor<K>,
		       32, 32> 
      <<<nr_total_blocks, B40C_RADIXSORT_THREADS>>>
      (mprts_cuda->d_bnd_spine_sums, mprts_cuda->d_bidx,
       mprts_cuda->d_ids, mprts_cuda->d_off, nr_total_blocks);
  } else if (b_mx[0] == 1 && b_mx[1] == 64 && b_mx[2] == 64) {
    ScanScatterDigits4<K, V, 0, RADIX_BITS, 0,
		       NopFunctor<K>,
		       NopFunctor<K>,
		       64, 64> 
      <<<nr_total_blocks, B40C_RADIXSORT_THREADS>>>
      (mprts_cuda->d_bnd_spine_sums, mprts_cuda->d_bidx,
       mprts_cuda->d_ids, mprts_cuda->d_off, nr_total_blocks);
  } else if (b_mx[0] == 1 && b_mx[1] == 128 && b_mx[2] == 128) {
    ScanScatterDigits4<K, V, 0, RADIX_BITS, 0,
                       NopFunctor<K>,
                       NopFunctor<K>,
                       128, 128>
      <<<nr_total_blocks, B40C_RADIXSORT_THREADS>>>
      (mprts_cuda->d_bnd_spine_sums, mprts_cuda->d_bidx,
       mprts_cuda->d_ids, mprts_cuda->d_off, nr_total_blocks);
  } else {
    mprintf("no support for b_mx %d x %d x %d!\n", b_mx[0], b_mx[1], b_mx[2]);
    assert(0);
  }
  cuda_sync_if_enabled();

  cuda_mprts_reorder_send_by_id(mprts);
}

