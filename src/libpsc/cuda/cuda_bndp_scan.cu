
#include "cuda_bndp.h"
#include "cuda_mparticles.h"
#include "cuda_bits.h"

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/scan.h>

#include "b40c/radixsort_scanscatter_kernel4.h"

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
//
// copies particles to be sent to contiguous area following
// the existing n_prts particles
//
// in: d_id[n_prts - n_prts_send:n_prts_send[
// in: d_xi4, d_pxi4[0:n_prts[
// out: d_xi4, d_pxi4[n_prts:n_prts_send[

template<typename CudaMparticles, typename DIM>
void cuda_bndp<CudaMparticles, DIM>::reorder_send_by_id(CudaMparticles* cmprts, uint n_prts_send)
{
  cmprts->d_xi4.resize(cmprts->n_prts + n_prts_send);
  cmprts->d_pxi4.resize(cmprts->n_prts + n_prts_send);
  
  if (n_prts_send == 0) {
    return;
  }

  int dimGrid = (n_prts_send + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;

  k_reorder_send_by_id<<<dimGrid, THREADS_PER_BLOCK>>>
    (n_prts_send, cmprts->by_block_.d_id.data().get() + cmprts->n_prts - n_prts_send,
     cmprts->d_xi4.data().get(), cmprts->d_pxi4.data().get(),
     cmprts->d_xi4.data().get() + cmprts->n_prts, cmprts->d_pxi4.data().get() + cmprts->n_prts);
  cuda_sync_if_enabled();
}

// ----------------------------------------------------------------------
// reorder_send_by_id_gold

template<typename CudaMparticles, typename DIM>
void cuda_bndp<CudaMparticles, DIM>::reorder_send_by_id_gold(CudaMparticles *cmprts, uint n_prts_send)
{
  thrust::host_vector<uint> h_id(cmprts->by_block_.d_id.data(), cmprts->by_block_.d_id.data() + cmprts->n_prts);
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
// k_reorder_send_buf_total

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

template<typename CudaMparticles, typename DIM>
void cuda_bndp<CudaMparticles, DIM>::reorder_send_buf_total(CudaMparticles *cmprts, uint n_prts_send)
{
  if (n_patches == 0) {
    return;
  }

  cmprts->resize(cmprts->n_prts + n_prts_send);

  float4 *xchg_xi4 = cmprts->d_xi4.data().get() + cmprts->n_prts;
  float4 *xchg_pxi4 = cmprts->d_pxi4.data().get() + cmprts->n_prts;
  
  dim3 dimBlock(THREADS_PER_BLOCK, 1);
  dim3 dimGrid((cmprts->n_prts + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK, 1);
  
  k_reorder_send_buf_total<<<dimGrid, dimBlock>>>(cmprts->n_prts, cmprts->n_blocks,
						  cmprts->by_block_.d_idx.data().get(), d_sums.data().get(),
						  cmprts->d_xi4.data().get(), cmprts->d_pxi4.data().get(),
						  xchg_xi4, xchg_pxi4);
  cuda_sync_if_enabled();
}

// ----------------------------------------------------------------------
// scan_send_buf_total
//
// in: d_spine_cnts
// out: d_spine_sums, d_id[n_prts - n_prts_send: n_prts[
//

template<typename CudaMparticles, typename DIM>
void cuda_bndp<CudaMparticles, DIM>::scan_send_buf_total(CudaMparticles* cmprts, uint n_prts_send)
{
  // OPT, we could do this from the beginning and adapt find_n_send()
  thrust::exclusive_scan(d_spine_cnts.data() + n_blocks * 10,
			 d_spine_cnts.data() + n_blocks * 11 + 1,
			 d_spine_sums.data() + n_blocks * 10,
			 cmprts->n_prts - n_prts_send);
  // OPT, we could somehow not fill in ids for not oob at all
  // this should make sure at least those within bounds don't screw anything up
  thrust::fill(d_spine_sums.data(), d_spine_sums.data() + n_blocks * 10, 0);

  Int3 mx = b_mx();
  if (mx[0] == 1 && mx[1] == 4 && mx[2] == 4) {
    ScanScatterDigits4<K, V, 0, RADIX_BITS, 0,
		       NopFunctor<K>,
		       NopFunctor<K>,
		       4, 4> 
      <<<n_blocks, B40C_RADIXSORT_THREADS>>>
      (d_spine_sums.data().get(), cmprts->by_block_.d_idx.data().get(), cmprts->by_block_.d_id.data().get(), cmprts->by_block_.d_off.data().get(), n_blocks);
  } else if (mx[0] == 1 && mx[1] == 8 && mx[2] == 8) {
    ScanScatterDigits4<K, V, 0, RADIX_BITS, 0,
		       NopFunctor<K>,
		       NopFunctor<K>,
		       8, 8> 
      <<<n_blocks, B40C_RADIXSORT_THREADS>>>
      (d_spine_sums.data().get(), cmprts->by_block_.d_idx.data().get(), cmprts->by_block_.d_id.data().get(), cmprts->by_block_.d_off.data().get(), n_blocks);
  } else if (mx[0] == 1 && mx[1] == 16 && mx[2] == 16) {
    ScanScatterDigits4<K, V, 0, RADIX_BITS, 0,
		       NopFunctor<K>,
		       NopFunctor<K>,
		       16, 16> 
      <<<n_blocks, B40C_RADIXSORT_THREADS>>>
      (d_spine_sums.data().get(), cmprts->by_block_.d_idx.data().get(), cmprts->by_block_.d_id.data().get(), cmprts->by_block_.d_off.data().get(), n_blocks);
  } else if (mx[0] == 1 && mx[1] == 32 && mx[2] == 32) {
    ScanScatterDigits4<K, V, 0, RADIX_BITS, 0,
		       NopFunctor<K>,
		       NopFunctor<K>,
		       32, 32> 
      <<<n_blocks, B40C_RADIXSORT_THREADS>>>
      (d_spine_sums.data().get(), cmprts->by_block_.d_idx.data().get(), cmprts->by_block_.d_id.data().get(), cmprts->by_block_.d_off.data().get(), n_blocks);
  } else if (mx[0] == 1 && mx[1] == 64 && mx[2] == 64) {
    ScanScatterDigits4<K, V, 0, RADIX_BITS, 0,
		       NopFunctor<K>,
		       NopFunctor<K>,
		       64, 64> 
      <<<n_blocks, B40C_RADIXSORT_THREADS>>>
      (d_spine_sums.data().get(), cmprts->by_block_.d_idx.data().get(), cmprts->by_block_.d_id.data().get(), cmprts->by_block_.d_off.data().get(), n_blocks);
  } else if (mx[0] == 1 && mx[1] == 128 && mx[2] == 128) {
    ScanScatterDigits4<K, V, 0, RADIX_BITS, 0,
                       NopFunctor<K>,
                       NopFunctor<K>,
                       128, 128>
      <<<n_blocks, B40C_RADIXSORT_THREADS>>>
      (d_spine_sums.data().get(), cmprts->by_block_.d_idx.data().get(), cmprts->by_block_.d_id.data().get(), cmprts->by_block_.d_off.data().get(), n_blocks);
  } else {
    printf("no support for b_mx %d x %d x %d!\n", mx[0], mx[1], mx[2]);
    assert(0);
  }
  cuda_sync_if_enabled();

  reorder_send_by_id(cmprts, n_prts_send);
}

// ----------------------------------------------------------------------
// scan_send_buf_total_gold

template<typename CudaMparticles, typename DIM>
void cuda_bndp<CudaMparticles, DIM>::scan_send_buf_total_gold(CudaMparticles *cmprts, uint n_prts_send)
{
  thrust::host_vector<uint> h_off(cmprts->by_block_.d_off);
  thrust::host_vector<uint> h_bidx(cmprts->by_block_.d_idx.data(), cmprts->by_block_.d_idx.data() + cmprts->n_prts);
  thrust::host_vector<uint> h_sums(cmprts->n_prts);
  
  for (uint bid = 0; bid < n_blocks; bid++) {
    uint sum = d_spine_sums[n_blocks * 10 + bid];
    for (int n = h_off[bid]; n < h_off[bid+1]; n++) {
      if (h_bidx[n] == CUDA_BND_S_OOB) {
	h_sums[n] = sum;
	sum++;
      }
    }
  }

  d_sums.resize(cmprts->n_prts);
  thrust::copy(h_sums.begin(), h_sums.end(), d_sums.begin());

  reorder_send_buf_total(cmprts, n_prts_send);
}

template struct cuda_bndp<cuda_mparticles<BS144>, dim_yz>;