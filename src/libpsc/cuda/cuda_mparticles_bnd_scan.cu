
#include "cuda_particles_bnd.h"
#include "cuda_mparticles.h"
#include "cuda_bits.h"

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/scan.h>

#include <b40c/radixsort_scanscatter_kernel4.h>

using namespace b40c_thrust;

typedef uint K;
typedef uint V;

static const int RADIX_BITS = 4;

#include <cstdio>
#include <cassert>

#define THREADS_PER_BLOCK 256

// ----------------------------------------------------------------------
// k_reorder_send_by_id

static void __global__
k_reorder_send_by_id(uint nr_prts_send, uint *d_xchg_ids,
			 float4 *d_xi4, float4 *d_pxi4,
			 float4 *d_xchg_xi4, float4 *d_xchg_pxi4)
{
  int n = threadIdx.x + THREADS_PER_BLOCK * blockIdx.x;
  if (n >= nr_prts_send) {
    return;
  }

  uint id = d_xchg_ids[n];
  d_xchg_xi4[n]  = d_xi4[id];
  d_xchg_pxi4[n] = d_pxi4[id];
}

// ----------------------------------------------------------------------
// reorder_send_by_id

void cuda_particles_bnd::reorder_send_by_id(struct cuda_mparticles *cmprts)
{
  if (cmprts->n_prts_send == 0) {
    return;
  }

  int dimGrid = (cmprts->n_prts_send + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;

  k_reorder_send_by_id<<<dimGrid, THREADS_PER_BLOCK>>>
    (cmprts->n_prts_send, cmprts->d_id.data().get() + cmprts->n_prts - cmprts->n_prts_send,
     cmprts->d_xi4.data().get(), cmprts->d_pxi4.data().get(),
     cmprts->d_xi4.data().get() + cmprts->n_prts, cmprts->d_pxi4.data().get() + cmprts->n_prts);
  cuda_sync_if_enabled();
}

// ----------------------------------------------------------------------
// reorder_send_by_id_gold

void cuda_particles_bnd::reorder_send_by_id_gold(cuda_mparticles *cmprts)
{
  uint n_prts_send = cmprts->n_prts_send;
  thrust::host_vector<uint> h_id(cmprts->d_id.data(), cmprts->d_id.data() + cmprts->n_prts);
  thrust::host_vector<float4> h_xi4(cmprts->d_xi4.data(), cmprts->d_xi4.data() + cmprts->n_prts + n_prts_send);
  thrust::host_vector<float4> h_pxi4(cmprts->d_pxi4.data(), cmprts->d_pxi4.data() + cmprts->n_prts + n_prts_send);
  
  for (int n = 0; n < n_prts_send; n++) {
    uint id = h_id[cmprts->n_prts - n_prts_send + n];
    h_xi4[cmprts->n_prts + n]  = h_xi4[id];
    h_pxi4[cmprts->n_prts + n] = h_pxi4[id];
  }

  thrust::copy(h_xi4.begin(), h_xi4.end(), cmprts->d_xi4.data());
  thrust::copy(h_pxi4.begin(), h_pxi4.end(), cmprts->d_pxi4.data());
}

// ----------------------------------------------------------------------
// k_send_buf_total

__global__ static void
k_reorder_send_buf_total(int nr_prts, int nr_total_blocks,
			     uint *d_bidx, uint *d_sums,
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

// ----------------------------------------------------------------------
// reorder_send_buf_total

void cuda_mparticles_bnd::reorder_send_buf_total(cuda_mparticles *cmprts)
{
  if (cmprts->n_patches == 0)
    return;

  float4 *xchg_xi4 = cmprts->d_xi4.data().get() + cmprts->n_prts;
  float4 *xchg_pxi4 = cmprts->d_pxi4.data().get() + cmprts->n_prts;
  assert(cmprts->n_prts + n_prts_send < cmprts->n_alloced);
  
  dim3 dimBlock(THREADS_PER_BLOCK, 1);
  dim3 dimGrid((cmprts->n_prts + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK, 1);
  
  k_reorder_send_buf_total<<<dimGrid, dimBlock>>>(cmprts->n_prts, cmprts->n_blocks,
						  cmprts->d_bidx.data().get(), d_sums.data().get(),
						  cmprts->d_xi4.data().get(), cmprts->d_pxi4.data().get(),
						  xchg_xi4, xchg_pxi4);
  cuda_sync_if_enabled();
}

// ----------------------------------------------------------------------
// scan_send_buf_total

void cuda_particles_bnd::scan_send_buf_total(struct cuda_mparticles *cmprts)
{
  thrust::device_vector<uint>& d_spine_cnts = cmprts->d_spine_cnts;
  thrust::device_vector<uint>& d_spine_sums = cmprts->d_spine_sums;

  uint n_blocks = cmprts->n_blocks;
  int *b_mx = cmprts->indexer.b_mx_;

  // OPT, we could do this from the beginning and adapt find_n_send()
  thrust::exclusive_scan(d_spine_cnts.data() + n_blocks * 10,
			 d_spine_cnts.data() + n_blocks * 11 + 1,
			 d_spine_sums.data() + n_blocks * 10,
			 cmprts->n_prts - cmprts->n_prts_send);
  // OPT, we could somehow not fill in ids for not oob at all
  // this should make sure at least those within bounds don't screw anything up
  thrust::fill(d_spine_sums.data(), d_spine_sums.data() + n_blocks * 10, 0);

  if (b_mx[0] == 1 && b_mx[1] == 4 && b_mx[2] == 4) {
    ScanScatterDigits4<K, V, 0, RADIX_BITS, 0,
		       NopFunctor<K>,
		       NopFunctor<K>,
		       4, 4> 
      <<<n_blocks, B40C_RADIXSORT_THREADS>>>
      (d_spine_sums.data().get(), cmprts->d_bidx.data().get(), cmprts->d_id.data().get(), cmprts->d_off.data().get(), n_blocks);
  } else if (b_mx[0] == 1 && b_mx[1] == 8 && b_mx[2] == 8) {
    ScanScatterDigits4<K, V, 0, RADIX_BITS, 0,
		       NopFunctor<K>,
		       NopFunctor<K>,
		       8, 8> 
      <<<n_blocks, B40C_RADIXSORT_THREADS>>>
      (d_spine_sums.data().get(), cmprts->d_bidx.data().get(), cmprts->d_id.data().get(), cmprts->d_off.data().get(), n_blocks);
  } else if (b_mx[0] == 1 && b_mx[1] == 16 && b_mx[2] == 16) {
    ScanScatterDigits4<K, V, 0, RADIX_BITS, 0,
		       NopFunctor<K>,
		       NopFunctor<K>,
		       16, 16> 
      <<<n_blocks, B40C_RADIXSORT_THREADS>>>
      (d_spine_sums.data().get(), cmprts->d_bidx.data().get(), cmprts->d_id.data().get(), cmprts->d_off.data().get(), n_blocks);
  } else if (b_mx[0] == 1 && b_mx[1] == 32 && b_mx[2] == 32) {
    ScanScatterDigits4<K, V, 0, RADIX_BITS, 0,
		       NopFunctor<K>,
		       NopFunctor<K>,
		       32, 32> 
      <<<n_blocks, B40C_RADIXSORT_THREADS>>>
      (d_spine_sums.data().get(), cmprts->d_bidx.data().get(), cmprts->d_id.data().get(), cmprts->d_off.data().get(), n_blocks);
  } else if (b_mx[0] == 1 && b_mx[1] == 64 && b_mx[2] == 64) {
    ScanScatterDigits4<K, V, 0, RADIX_BITS, 0,
		       NopFunctor<K>,
		       NopFunctor<K>,
		       64, 64> 
      <<<n_blocks, B40C_RADIXSORT_THREADS>>>
      (d_spine_sums.data().get(), cmprts->d_bidx.data().get(), cmprts->d_id.data().get(), cmprts->d_off.data().get(), n_blocks);
  } else if (b_mx[0] == 1 && b_mx[1] == 128 && b_mx[2] == 128) {
    ScanScatterDigits4<K, V, 0, RADIX_BITS, 0,
                       NopFunctor<K>,
                       NopFunctor<K>,
                       128, 128>
      <<<n_blocks, B40C_RADIXSORT_THREADS>>>
      (d_spine_sums.data().get(), cmprts->d_bidx.data().get(), cmprts->d_id.data().get(), cmprts->d_off.data().get(), n_blocks);
  } else {
    printf("no support for b_mx %d x %d x %d!\n", b_mx[0], b_mx[1], b_mx[2]);
    assert(0);
  }
  cuda_sync_if_enabled();

  reorder_send_by_id(cmprts);
}

// ----------------------------------------------------------------------
// scan_send_buf_total_gold

void cuda_particles_bnd::scan_send_buf_total_gold(cuda_mparticles *cmprts)
{
  thrust::device_vector<uint>& d_spine_sums = cmprts->d_spine_sums;

  uint n_blocks = cmprts->n_blocks;

  thrust::host_vector<uint> h_off(cmprts->d_off);
  thrust::host_vector<uint> h_bidx(cmprts->d_bidx.data(), cmprts->d_bidx.data() + cmprts->n_prts);
  thrust::host_vector<uint> h_sums(cmprts->d_sums.data(), cmprts->d_sums.data() + cmprts->n_prts);
  
  for (uint bid = 0; bid < n_blocks; bid++) {
    uint sum = d_spine_sums[n_blocks * 10 + bid];
    for (int n = h_off[bid]; n < h_off[bid+1]; n++) {
      if (h_bidx[n] == CUDA_BND_S_OOB) {
	h_sums[n] = sum;
	sum++;
      }
    }
  }

  thrust::copy(h_sums.begin(), h_sums.end(), cmprts->d_sums.begin());

  cmprts->reorder_send_buf_total(cmprts);
}

