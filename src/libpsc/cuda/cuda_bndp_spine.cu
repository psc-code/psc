
#include "cuda_bndp.h"
#include "cuda_mparticles.h"
#include "cuda_bits.h"

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/scan.h>

#include "b40c/radixsort_reduction_kernel.h"
#include "b40c/radixsort_scanscatter_kernel3.h"

#include <mrc_profile.h>

using namespace b40c_thrust;

typedef uint K;
typedef uint V;

static const int RADIX_BITS = 4;

#define THREADS_PER_BLOCK 256

// ----------------------------------------------------------------------
// spine_reduce

template<typename CudaMparticles, typename DIM>
void cuda_bndp<CudaMparticles, DIM>::spine_reduce(CudaMparticles *cmprts)
{
  // OPT?
  thrust::fill(d_spine_cnts.data(), d_spine_cnts.data() + 1 + n_blocks * (CUDA_BND_STRIDE + 1), 0);

  const int threads = B40C_RADIXSORT_THREADS;
  Int3 mx = b_mx();
  if (mx[0] == 1 && mx[1] == 2 && mx[2] == 2) {
    RakingReduction3x<K, V, 0, RADIX_BITS, 0,
		      NopFunctor<K>, 2, 2> <<<n_blocks, threads>>>
      (d_spine_cnts.data().get(), cmprts->by_block_.d_idx.data().get(), cmprts->by_block_.d_off.data().get(), n_blocks);
  } else if (mx[0] == 1 && mx[1] == 4 && mx[2] == 4) {
    RakingReduction3x<K, V, 0, RADIX_BITS, 0,
		      NopFunctor<K>, 4, 4> <<<n_blocks, threads>>>
      (d_spine_cnts.data().get(), cmprts->by_block_.d_idx.data().get(), cmprts->by_block_.d_off.data().get(), n_blocks);
  } else if (mx[0] == 1 && mx[1] == 8 && mx[2] == 8) {
    RakingReduction3x<K, V, 0, RADIX_BITS, 0,
		      NopFunctor<K>, 8, 8> <<<n_blocks, threads>>>
      (d_spine_cnts.data().get(), cmprts->by_block_.d_idx.data().get(), cmprts->by_block_.d_off.data().get(), n_blocks);
  } else if (mx[0] == 1 && mx[1] == 16 && mx[2] == 16) {
    RakingReduction3x<K, V, 0, RADIX_BITS, 0,
		      NopFunctor<K>, 16, 16> <<<n_blocks, threads>>>
      (d_spine_cnts.data().get(), cmprts->by_block_.d_idx.data().get(), cmprts->by_block_.d_off.data().get(), n_blocks);
  } else if (mx[0] == 1 && mx[1] == 32 && mx[2] == 32) {
    RakingReduction3x<K, V, 0, RADIX_BITS, 0,
		      NopFunctor<K>, 32, 32> <<<n_blocks, threads>>>
      (d_spine_cnts.data().get(), cmprts->by_block_.d_idx.data().get(), cmprts->by_block_.d_off.data().get(), n_blocks);
  } else if (mx[0] == 1 && mx[1] == 64 && mx[2] == 64) {
    RakingReduction3x<K, V, 0, RADIX_BITS, 0,
		      NopFunctor<K>, 64, 64> <<<n_blocks, threads>>>
      (d_spine_cnts.data().get(), cmprts->by_block_.d_idx.data().get(), cmprts->by_block_.d_off.data().get(), n_blocks);
  } else if (mx[0] == 1 && mx[1] == 128 && mx[2] == 128) {
    RakingReduction3x<K, V, 0, RADIX_BITS, 0,
                      NopFunctor<K>, 128, 128> <<<n_blocks, threads>>>
      (d_spine_cnts.data().get(), cmprts->by_block_.d_idx.data().get(), cmprts->by_block_.d_off.data().get(), n_blocks);
  } else {
    printf("no support for b_mx %d x %d x %d!\n", mx[0], mx[1], mx[2]);
    assert(0);
  }
  cuda_sync_if_enabled();

  thrust::exclusive_scan(d_spine_cnts.data() + n_blocks * 10,
			 d_spine_cnts.data() + n_blocks * 10 + n_blocks + 1,
			 d_spine_sums.data() + n_blocks * 10);
}

// ----------------------------------------------------------------------
// cuda_mprts_spine_reduce_gold

template<typename CudaMparticles, typename DIM>
void cuda_bndp<CudaMparticles, DIM>::spine_reduce_gold(CudaMparticles *cmprts)
{
  thrust::fill(d_spine_cnts.data(), d_spine_cnts.data() + 1 + n_blocks * (CUDA_BND_STRIDE + 1), 0);

  thrust::host_vector<uint> h_bidx(cmprts->by_block_.d_idx.data(), cmprts->by_block_.d_idx.data() + cmprts->n_prts);
  thrust::host_vector<uint> h_off(cmprts->by_block_.d_off);
  thrust::host_vector<uint> h_spine_cnts(d_spine_cnts.data(), d_spine_cnts.data() + 1 + n_blocks * (CUDA_BND_STRIDE + 1));

  Int3 mx = b_mx();
  for (int p = 0; p < n_patches; p++) {
    for (int b = 0; b < n_blocks_per_patch; b++) {
      uint bid = b + p * n_blocks_per_patch;
      for (int n = h_off[bid]; n < h_off[bid+1]; n++) {
	uint key = h_bidx[n];
	if (key < 9) {
	  int dy = key % 3;
	  int dz = key / 3;
	  int by = b % mx[1];
	  int bz = b / mx[1];
	  uint bby = by + 1 - dy;
	  uint bbz = bz + 1 - dz;
	  uint bb = bbz * mx[1] + bby;
	  if (bby < mx[1] && bbz < mx[2]) {
	    h_spine_cnts[(bb + p * n_blocks_per_patch) * 10 + key]++;
	  } else {
	    assert(0);
	  }
	} else if (key == CUDA_BND_S_OOB) {
	  h_spine_cnts[mx[1]*mx[2]*n_patches * 10 + bid]++;
	}
      }
    }
  }  

  thrust::copy(h_spine_cnts.begin(), h_spine_cnts.end(), d_spine_cnts.begin());
  thrust::exclusive_scan(d_spine_cnts.data() + n_blocks * 10,
			 d_spine_cnts.data() + n_blocks * 10 + n_blocks + 1,
			 d_spine_sums.data() + n_blocks * 10);
}

// ----------------------------------------------------------------------
// k_count_received

__global__ static void
k_count_received(int nr_total_blocks, uint *d_n_recv_by_block, uint *d_spine_cnts)
{
  int bid = threadIdx.x + THREADS_PER_BLOCK * blockIdx.x;

  if (bid < nr_total_blocks) {
    d_spine_cnts[bid * 10 + CUDA_BND_S_NEW] = d_n_recv_by_block[bid];
  }
}

// ----------------------------------------------------------------------
// count_received

template<typename CudaMparticles, typename DIM>
void cuda_bndp<CudaMparticles, DIM>::count_received(CudaMparticles *cmprts)
{
  k_count_received<<<n_blocks, THREADS_PER_BLOCK>>>
    (n_blocks, d_spine_cnts.data().get() + 10 * n_blocks, d_spine_cnts.data().get());
}

// ----------------------------------------------------------------------
// count_received_gold

template<typename CudaMparticles, typename DIM>
void cuda_bndp<CudaMparticles, DIM>::count_received_gold(CudaMparticles *cmprts)
{
  thrust::host_vector<uint> h_spine_cnts(1 + n_blocks * (10 + 1));

  thrust::copy(d_spine_cnts.data(), d_spine_cnts.data() + 1 + n_blocks * (10 + 1), h_spine_cnts.begin());

  for (int bid = 0; bid < n_blocks; bid++) {
    h_spine_cnts[bid * 10 + CUDA_BND_S_NEW] = h_spine_cnts[10 * n_blocks + bid];
  }

  thrust::copy(h_spine_cnts.begin(), h_spine_cnts.end(), d_spine_cnts.begin());
}

#if 0
void cuda_bndp::count_received_v1(CudaMparticles *cmprts)
{
  thrust::device_ptr<uint> d_bidx(cmprts->d_bidx);
  thrust::device_ptr<uint> d_spine_cnts(d_bnd_spine_cnts);

  thrust::host_vector<uint> h_bidx(cmprts->n_prts);
  thrust::host_vector<uint> h_spine_cnts(1 + n_blocks * (10 + 1));

  thrust::copy(d_bidx, d_bidx + cmprts->n_prts, h_bidx.begin());
  thrust::copy(d_spine_cnts, d_spine_cnts + 1 + n_blocks * (10 + 1), h_spine_cnts.begin());
  for (int n = cmprts->n_prts - n_prts_recv; n < cmprts->n_prts; n++) {
    assert(h_bidx[n] < n_blocks);
    h_spine_cnts[h_bidx[n] * 10 + CUDA_BND_S_NEW]++;
  }
  thrust::copy(h_spine_cnts.begin(), h_spine_cnts.end(), d_spine_cnts.begin());
}
#endif

// ----------------------------------------------------------------------
// k_scan_scatter_received

static void __global__
k_scan_scatter_received(uint nr_recv, uint nr_prts_prev,
			    uint *d_spine_sums, uint *d_bnd_off,
			    uint *d_bidx, uint *d_ids)
{
  int n0 = threadIdx.x + THREADS_PER_BLOCK * blockIdx.x;
  if (n0 >= nr_recv) {
    return;
  }

  int n = n0 + nr_prts_prev;

  int nn = d_spine_sums[d_bidx[n] * 10 + CUDA_BND_S_NEW] + d_bnd_off[n0];
  d_ids[nn] = n;
}

// ----------------------------------------------------------------------
// scan_scatter_received

template<typename CudaMparticles, typename DIM>
void cuda_bndp<CudaMparticles, DIM>::scan_scatter_received(CudaMparticles *cmprts, uint n_prts_recv)
{
  if (n_prts_recv == 0) {
    return;
  }
  
  uint n_prts_prev = cmprts->n_prts - n_prts_recv;

  int dimGrid = (n_prts_recv + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;

  k_scan_scatter_received<<<dimGrid, THREADS_PER_BLOCK>>>
    (n_prts_recv, n_prts_prev, d_spine_sums.data().get(), d_bnd_off.data().get(),
     cmprts->by_block_.d_idx.data().get(), cmprts->by_block_.d_id.data().get());
  cuda_sync_if_enabled();
}

// ----------------------------------------------------------------------
// scan_scatter_received_gold

template<typename CudaMparticles, typename DIM>
void cuda_bndp<CudaMparticles, DIM>::scan_scatter_received_gold(CudaMparticles *cmprts, uint n_prts_recv)
{
  thrust::host_vector<uint> h_bidx(cmprts->n_prts);
  thrust::host_vector<uint> h_bnd_off(n_prts_recv);
  thrust::host_vector<uint> h_id(cmprts->n_prts);
  thrust::host_vector<uint> h_spine_sums(1 + n_blocks * (10 + 1));

  thrust::copy(d_spine_sums.data(), d_spine_sums.data() + n_blocks * 11, h_spine_sums.begin());
  thrust::copy(cmprts->by_block_.d_idx.data(), cmprts->by_block_.d_idx.data() + cmprts->n_prts, h_bidx.begin());

  uint n_prts_prev = cmprts->n_prts - n_prts_recv;
  thrust::copy(d_bnd_off.begin(), d_bnd_off.end(), h_bnd_off.begin());
  for (int n0 = 0; n0 < n_prts_recv; n0++) {
    int n = n0 + n_prts_prev;
    int nn = h_spine_sums[h_bidx[n] * 10 + CUDA_BND_S_NEW] + h_bnd_off[n0];
    h_id[nn] = n;
  }
  thrust::copy(h_id.begin(), h_id.end(), cmprts->by_block_.d_id.begin());
}

// ----------------------------------------------------------------------
// sort_pairs_device

template<typename CudaMparticles, typename DIM>
void cuda_bndp<CudaMparticles, DIM>::sort_pairs_device(CudaMparticles *cmprts, uint n_prts_recv)
{
  static int pr_A, pr_B, pr_C, pr_D;
  if (!pr_B) {
    pr_A = prof_register("xchg_cnt_recvd", 1., 0, 0);
    pr_B = prof_register("xchg_top_scan", 1., 0, 0);
    pr_C = prof_register("xchg_ss_recvd", 1., 0, 0);
    pr_D = prof_register("xchg_bottom_scan", 1., 0, 0);
  }

  prof_start(pr_A);
  count_received(cmprts);
  prof_stop(pr_A);

  prof_start(pr_B);
  // FIXME why isn't 10 + 0 enough?
  thrust::exclusive_scan(d_spine_cnts.data(), d_spine_cnts.data() + 1 + n_blocks * (10 + 1), d_spine_sums.data());
  prof_stop(pr_B);

  prof_start(pr_C);
  scan_scatter_received(cmprts, n_prts_recv);
  prof_stop(pr_C);

  prof_start(pr_D);
  Int3 mx = b_mx();
  if (mx[0] == 1 && mx[1] == 4 && mx[2] == 4) {
    ScanScatterDigits3x<K, V, 0, RADIX_BITS, 0,
			NopFunctor<K>,
			NopFunctor<K>,
			4, 4> 
      <<<n_blocks, B40C_RADIXSORT_THREADS>>>
      (d_spine_sums.data().get(), cmprts->by_block_.d_idx.data().get(), cmprts->by_block_.d_id.data().get(), cmprts->by_block_.d_off.data().get(), n_blocks);
  } else if (mx[0] == 1 && mx[1] == 8 && mx[2] == 8) {
    ScanScatterDigits3x<K, V, 0, RADIX_BITS, 0,
			NopFunctor<K>,
			NopFunctor<K>,
			8, 8> 
      <<<n_blocks, B40C_RADIXSORT_THREADS>>>
      (d_spine_sums.data().get(), cmprts->by_block_.d_idx.data().get(), cmprts->by_block_.d_id.data().get(), cmprts->by_block_.d_off.data().get(), n_blocks);
  } else if (mx[0] == 1 && mx[1] == 16 && mx[2] == 16) {
    ScanScatterDigits3x<K, V, 0, RADIX_BITS, 0,
			NopFunctor<K>,
			NopFunctor<K>,
			16, 16> 
      <<<n_blocks, B40C_RADIXSORT_THREADS>>>
      (d_spine_sums.data().get(), cmprts->by_block_.d_idx.data().get(), cmprts->by_block_.d_id.data().get(), cmprts->by_block_.d_off.data().get(), n_blocks);
  } else if (mx[0] == 1 && mx[1] == 32 && mx[2] == 32) {
    ScanScatterDigits3x<K, V, 0, RADIX_BITS, 0,
			NopFunctor<K>,
			NopFunctor<K>,
			32, 32> 
      <<<n_blocks, B40C_RADIXSORT_THREADS>>>
      (d_spine_sums.data().get(), cmprts->by_block_.d_idx.data().get(), cmprts->by_block_.d_id.data().get(), cmprts->by_block_.d_off.data().get(), n_blocks);
  } else if (mx[0] == 1 && mx[1] == 64 && mx[2] == 64) {
    ScanScatterDigits3x<K, V, 0, RADIX_BITS, 0,
			NopFunctor<K>,
			NopFunctor<K>,
			64, 64> 
      <<<n_blocks, B40C_RADIXSORT_THREADS>>>
      (d_spine_sums.data().get(), cmprts->by_block_.d_idx.data().get(), cmprts->by_block_.d_id.data().get(), cmprts->by_block_.d_off.data().get(), n_blocks);
  } else if (mx[0] == 1 && mx[1] == 128 && mx[2] == 128) {
    ScanScatterDigits3x<K, V, 0, RADIX_BITS, 0,
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
  prof_stop(pr_D);

  // d_ids now contains the indices to reorder by
}

template<typename CudaMparticles, typename DIM>
void cuda_bndp<CudaMparticles, DIM>::sort_pairs_gold(CudaMparticles *cmprts, uint n_prts_recv)
{
  thrust::host_vector<uint> h_bidx(cmprts->by_block_.d_idx.data(), cmprts->by_block_.d_idx.data() + cmprts->n_prts);
  thrust::host_vector<uint> h_id(cmprts->n_prts);
  thrust::host_vector<uint> h_off(cmprts->by_block_.d_off);
  thrust::host_vector<uint> h_spine_cnts(d_spine_cnts.data(), d_spine_cnts.data() + 1 + n_blocks * (10 + 1));

  thrust::host_vector<uint> h_spine_sums(1 + n_blocks * (10 + 1));

  for (int n = cmprts->n_prts - n_prts_recv; n < cmprts->n_prts; n++) {
    assert(h_bidx[n] < n_blocks);
    h_spine_cnts[h_bidx[n] * 10 + CUDA_BND_S_NEW]++;
  }

  thrust::exclusive_scan(h_spine_cnts.begin(), h_spine_cnts.end(), h_spine_sums.begin());
  thrust::copy(h_spine_sums.begin(), h_spine_sums.end(), d_spine_sums.begin());

  Int3 mx = b_mx();
  for (int bid = 0; bid < n_blocks; bid++) {
    int b = bid % n_blocks_per_patch;
    int p = bid / n_blocks_per_patch;
    for (int n = h_off[bid]; n < h_off[bid+1]; n++) {
      uint key = h_bidx[n];
      if (key < 9) {
	int dy = key % 3;
	int dz = key / 3;
	int by = b % mx[1];
	int bz = b / mx[1];
	uint bby = by + 1 - dy;
	uint bbz = bz + 1 - dz;
	assert(bby < mx[1] && bbz < mx[2]);
	uint bb = bbz * mx[1] + bby;
	int nn = h_spine_sums[(bb + p * n_blocks_per_patch) * 10 + key]++;
	h_id[nn] = n;
      } else { // OOB
	assert(0);
      }
    }
  }
  for (int n = cmprts->n_prts - n_prts_recv; n < cmprts->n_prts; n++) {
      int nn = h_spine_sums[h_bidx[n] * 10 + CUDA_BND_S_NEW]++;
      h_id[nn] = n;
  }

  thrust::copy(h_id.begin(), h_id.end(), cmprts->by_block_.d_id.begin());
  // d_ids now contains the indices to reorder by
}

template struct cuda_bndp<cuda_mparticles<BS144>, dim_yz>;