
#include "cuda_mparticles.h"
#include "cuda_bits.h"

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/scan.h>

#include <b40c/radixsort_scanscatter_kernel4.h>

using namespace b40c_thrust;

typedef unsigned int K;
typedef unsigned int V;

static const int RADIX_BITS = 4;

#include <cstdio>
#include <cassert>

#define THREADS_PER_BLOCK 256

// ----------------------------------------------------------------------
// cuda_mparticles_reorder_send_by_id

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
cuda_mparticles_reorder_send_by_id(struct cuda_mparticles *cmprts)
{
  if (cmprts->bnd.n_prts_send == 0) {
    return;
  }

  int dimGrid = (cmprts->bnd.n_prts_send + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;

  mprts_reorder_send_by_id<<<dimGrid, THREADS_PER_BLOCK>>>
    (cmprts->bnd.n_prts_send, cmprts->d_id + cmprts->n_prts - cmprts->bnd.n_prts_send,
     cmprts->d_xi4, cmprts->d_pxi4,
     cmprts->d_xi4 + cmprts->n_prts, cmprts->d_pxi4 + cmprts->n_prts);
  cuda_sync_if_enabled();
}

// ----------------------------------------------------------------------
// cuda_mparticles_reorder_send_by_id_gold

void
cuda_mparticles_reorder_send_by_id_gold(struct cuda_mparticles *cmprts)
{
  thrust::device_ptr<unsigned int> d_id(cmprts->d_id);
  thrust::device_ptr<float4> d_xi4(cmprts->d_xi4);
  thrust::device_ptr<float4> d_pxi4(cmprts->d_pxi4);
  thrust::host_vector<unsigned int> h_id(d_id, d_id + cmprts->n_prts);
  thrust::host_vector<float4> h_xi4(d_xi4, d_xi4 + cmprts->n_prts + cmprts->bnd.n_prts_send);
  thrust::host_vector<float4> h_pxi4(d_pxi4, d_pxi4 + cmprts->n_prts + cmprts->bnd.n_prts_send);
  
  for (int n = 0; n < cmprts->bnd.n_prts_send; n++) {
    unsigned int id = h_id[cmprts->n_prts - cmprts->bnd.n_prts_send + n];
    h_xi4[cmprts->n_prts + n]  = h_xi4[id];
    h_pxi4[cmprts->n_prts + n] = h_pxi4[id];
  }

  thrust::copy(h_xi4.begin(), h_xi4.end(), d_xi4);
  thrust::copy(h_pxi4.begin(), h_pxi4.end(), d_pxi4);
}

// ----------------------------------------------------------------------
// cuda_mparticles_reorder_send_buf_total

__global__ static void
mprts_reorder_send_buf_total(int nr_prts, int nr_total_blocks,
			     unsigned int *d_bidx, unsigned int *d_sums,
			     float4 *d_xi4, float4 *d_pxi4,
			     float4 *d_xchg_xi4, float4 *d_xchg_pxi4)
{
  int i = threadIdx.x + THREADS_PER_BLOCK * blockIdx.x;
  if (i >= nr_prts)
    return;

  if (d_bidx[i] == CUDA_BND_S_OOB) {
    int j = d_sums[i];
    d_xchg_xi4[j]  = d_xi4[i];
    d_xchg_pxi4[j] = d_pxi4[i];
  }
}

void
cuda_mparticles_reorder_send_buf_total(struct cuda_mparticles *cmprts)
{
  if (cmprts->n_patches == 0)
    return;

  float4 *xchg_xi4 = cmprts->d_xi4 + cmprts->n_prts;
  float4 *xchg_pxi4 = cmprts->d_pxi4 + cmprts->n_prts;
  assert(cmprts->n_prts + cmprts->bnd.n_prts_send < cmprts->n_alloced);
  
  dim3 dimBlock(THREADS_PER_BLOCK, 1);
  dim3 dimGrid((cmprts->n_prts + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK, 1);
  
  mprts_reorder_send_buf_total<<<dimGrid, dimBlock>>>(cmprts->n_prts, cmprts->n_blocks,
						      cmprts->d_bidx, cmprts->bnd.d_sums,
						      cmprts->d_xi4, cmprts->d_pxi4,
						      xchg_xi4, xchg_pxi4);
  cuda_sync_if_enabled();
}

// ----------------------------------------------------------------------
// cuda_mparticles_scan_send_buf_total

void
cuda_mparticles_scan_send_buf_total(struct cuda_mparticles *cmprts)
{
  unsigned int n_blocks = cmprts->n_blocks;
  int *b_mx = cmprts->b_mx;

  // OPT, we could do this from the beginning and adapt find_n_send()
  thrust::device_ptr<unsigned int> d_spine_cnts(cmprts->bnd.d_bnd_spine_cnts);
  thrust::device_ptr<unsigned int> d_spine_sums(cmprts->bnd.d_bnd_spine_sums);
  thrust::exclusive_scan(d_spine_cnts + n_blocks * 10,
			 d_spine_cnts + n_blocks * 11 + 1,
			 d_spine_sums + n_blocks * 10,
			 cmprts->n_prts - cmprts->bnd.n_prts_send);
  // OPT, we could somehow not fill in ids for not oob at all
  // this should make sure at least those within bounds don't screw anything up
  thrust::fill(d_spine_sums, d_spine_sums + n_blocks * 10, 0);

  if (b_mx[0] == 1 && b_mx[1] == 8 && b_mx[2] == 8) {
    ScanScatterDigits4<K, V, 0, RADIX_BITS, 0,
		       NopFunctor<K>,
		       NopFunctor<K>,
		       8, 8> 
      <<<n_blocks, B40C_RADIXSORT_THREADS>>>
      (cmprts->bnd.d_bnd_spine_sums, cmprts->d_bidx,
       cmprts->d_id, cmprts->d_off, n_blocks);
  } else if (b_mx[0] == 1 && b_mx[1] == 16 && b_mx[2] == 16) {
    ScanScatterDigits4<K, V, 0, RADIX_BITS, 0,
		       NopFunctor<K>,
		       NopFunctor<K>,
		       16, 16> 
      <<<n_blocks, B40C_RADIXSORT_THREADS>>>
      (cmprts->bnd.d_bnd_spine_sums, cmprts->d_bidx,
       cmprts->d_id, cmprts->d_off, n_blocks);
  } else if (b_mx[0] == 1 && b_mx[1] == 32 && b_mx[2] == 32) {
    ScanScatterDigits4<K, V, 0, RADIX_BITS, 0,
		       NopFunctor<K>,
		       NopFunctor<K>,
		       32, 32> 
      <<<n_blocks, B40C_RADIXSORT_THREADS>>>
      (cmprts->bnd.d_bnd_spine_sums, cmprts->d_bidx,
       cmprts->d_id, cmprts->d_off, n_blocks);
  } else if (b_mx[0] == 1 && b_mx[1] == 64 && b_mx[2] == 64) {
    ScanScatterDigits4<K, V, 0, RADIX_BITS, 0,
		       NopFunctor<K>,
		       NopFunctor<K>,
		       64, 64> 
      <<<n_blocks, B40C_RADIXSORT_THREADS>>>
      (cmprts->bnd.d_bnd_spine_sums, cmprts->d_bidx,
       cmprts->d_id, cmprts->d_off, n_blocks);
  } else if (b_mx[0] == 1 && b_mx[1] == 128 && b_mx[2] == 128) {
    ScanScatterDigits4<K, V, 0, RADIX_BITS, 0,
                       NopFunctor<K>,
                       NopFunctor<K>,
                       128, 128>
      <<<n_blocks, B40C_RADIXSORT_THREADS>>>
      (cmprts->bnd.d_bnd_spine_sums, cmprts->d_bidx,
       cmprts->d_id, cmprts->d_off, n_blocks);
  } else {
    printf("no support for b_mx %d x %d x %d!\n", b_mx[0], b_mx[1], b_mx[2]);
    assert(0);
  }
  cuda_sync_if_enabled();

  cuda_mparticles_reorder_send_by_id(cmprts);
}

// ----------------------------------------------------------------------
// cuda_mprts_scan_send_buf_total_gold

void
cuda_mparticles_scan_send_buf_total_gold(struct cuda_mparticles *cmprts)
{
  unsigned int n_blocks = cmprts->n_blocks;

  thrust::device_ptr<unsigned int> d_bidx(cmprts->d_bidx);
  thrust::device_ptr<unsigned int> d_sums(cmprts->bnd.d_sums);
  thrust::device_ptr<unsigned int> d_off(cmprts->d_off);
  thrust::device_ptr<unsigned int> d_spine_sums(cmprts->bnd.d_bnd_spine_sums);
  thrust::host_vector<unsigned int> h_off(d_off, d_off + n_blocks + 1);
  thrust::host_vector<unsigned int> h_bidx(d_bidx, d_bidx + cmprts->n_prts);
  thrust::host_vector<unsigned int> h_sums(d_sums, d_sums + cmprts->n_prts);
  
  for (unsigned int bid = 0; bid < n_blocks; bid++) {
    unsigned int sum = d_spine_sums[n_blocks * 10 + bid];
    for (int n = h_off[bid]; n < h_off[bid+1]; n++) {
      if (h_bidx[n] == CUDA_BND_S_OOB) {
	h_sums[n] = sum;
	sum++;
      }
    }
  }

  thrust::copy(h_sums.begin(), h_sums.end(), d_sums);

  cuda_mparticles_reorder_send_buf_total(cmprts);
}

