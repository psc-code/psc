
#undef _GLIBCXX_USE_INT128

#include "cuda_mparticles.h"

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/scan.h>

#include "psc_cuda.h"
#include "particles_cuda.h"
#include <mrc_profile.h>

#include <b40c/radixsort_reduction_kernel.h>
#include <b40c/radixsort_scanscatter_kernel3.h>

using namespace b40c_thrust;

typedef unsigned int K;
typedef unsigned int V;

static const int RADIX_BITS = 4;

// layout of the spine
//     lt             self             rb        # from left-top .. self .. right-bottom 
//     0   1   2   3   4   5   6   7   8   NEW
// b0 |   |   |   |   |   |   |   |   |   |   |
// b1 |   |   |   |   |   |   |   |   |   |   |
// b2 |   |   |   |   |   |   |   |   |   |   |
// ...
// bn |   |   |   |   |   |   |   |   |   |   |

//    |   |   |   |   |   |   |   |   |   |   |   |   | ... |   | # oob
//     b0  b1  b2  b3                                        bn

// ======================================================================
// cuda_mprts_bidx_to_key

template<int NBLOCKS_X, int NBLOCKS_Y, int NBLOCKS_Z>
__global__ static void
mprts_bidx_to_key(int nr_total_blocks, unsigned int *d_off, unsigned int *d_bidx)
{
  int tid = threadIdx.x;
  int bid = blockIdx.x;

  int block_begin = d_off[bid];
  int block_end   = d_off[bid + 1];

  for (int n = block_begin + tid; n < block_end; n += THREADS_PER_BLOCK) {
    unsigned int bidx = d_bidx[n];
    if (bidx == nr_total_blocks) {
      bidx = CUDA_BND_S_OOB;
    } else {
      int b_diff = bid - bidx + NBLOCKS_Y + 1;
      int d1 = b_diff % NBLOCKS_Y;
      int d2 = b_diff / NBLOCKS_Y;
      bidx = d2 * 3 + d1;
    }
    d_bidx[n] = bidx;
  }
}

void
cuda_mprts_bidx_to_key(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);
  struct cuda_mparticles *cmprts = mprts_cuda->cmprts;

  unsigned int n_blocks = cmprts->n_blocks;
  
  int *b_mx = mprts_cuda->b_mx;
  if (b_mx[0] == 1 && b_mx[1] == 8 && b_mx[2] == 8) {
    mprts_bidx_to_key<1, 8, 8><<<n_blocks, THREADS_PER_BLOCK>>>
      (n_blocks, mprts_cuda->d_off, cmprts->d_bidx);
  } else {
    mprintf("no support for b_mx %d x %d x %d!\n", b_mx[0], b_mx[1], b_mx[2]);
    assert(0);
  }
}

void
cuda_mprts_bidx_to_key_gold(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);
  struct cuda_mparticles *cmprts = mprts_cuda->cmprts;

  unsigned int n_blocks = cmprts->n_blocks;
  unsigned int n_blocks_per_patch = cmprts->n_blocks_per_patch;

  thrust::device_ptr<unsigned int> d_bidx(cmprts->d_bidx);
  thrust::device_ptr<unsigned int> d_off(mprts_cuda->d_off);

  thrust::host_vector<unsigned int> h_bidx(d_bidx, d_bidx + cmprts->n_prts);
  thrust::host_vector<unsigned int> h_off(d_off, d_off + n_blocks + 1);

  int *b_mx = mprts_cuda->b_mx;

  for (int bid = 0; bid < n_blocks; bid++) {
    int p = bid / n_blocks_per_patch;
    for (int n = h_off[bid]; n < h_off[bid+1]; n++) {
      assert((h_bidx[n] >= p * n_blocks_per_patch && h_bidx[n] < (p+1) * n_blocks_per_patch) ||
	     (h_bidx[n] == n_blocks));
      int bidx;
      if (h_bidx[n] == n_blocks) {
	bidx = CUDA_BND_S_OOB;
      } else {
	int b_diff = bid - h_bidx[n] + b_mx[1] + 1;
	int d1 = b_diff % b_mx[1];
	int d2 = b_diff / b_mx[1];
	bidx = d2 * 3 + d1;
      }
      
      h_bidx[n] = bidx;
    }
  }

  thrust::copy(h_bidx.begin(), h_bidx.end(), d_bidx);
}

// ======================================================================
// cuda_mprts_count_received

__global__ static void
mprts_count_received(int nr_total_blocks, unsigned int *d_alt_bidx, unsigned int *d_spine_cnts)
{
  int bid = threadIdx.x + THREADS_PER_BLOCK * blockIdx.x;

  if (bid < nr_total_blocks) {
    d_spine_cnts[bid * 10 + CUDA_BND_S_NEW] = d_alt_bidx[bid];
  }
}

void
cuda_mprts_count_received(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);
  struct cuda_mparticles *cmprts = mprts_cuda->cmprts;

  unsigned int n_blocks = cmprts->n_blocks;
  
  mprts_count_received<<<n_blocks, THREADS_PER_BLOCK>>>
    (n_blocks, mprts_cuda->d_bnd_spine_cnts + 10 * n_blocks, mprts_cuda->d_bnd_spine_cnts);
}

void
cuda_mprts_count_received_gold(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);
  struct cuda_mparticles *cmprts = mprts_cuda->cmprts;

  int n_blocks = cmprts->n_blocks;

  thrust::device_ptr<unsigned int> d_spine_cnts(mprts_cuda->d_bnd_spine_cnts);

  thrust::host_vector<unsigned int> h_spine_cnts(1 + n_blocks * (10 + 1));

  thrust::copy(d_spine_cnts, d_spine_cnts + 1 + n_blocks * (10 + 1), h_spine_cnts.begin());

  for (int bid = 0; bid < n_blocks; bid++) {
    h_spine_cnts[bid * 10 + CUDA_BND_S_NEW] = h_spine_cnts[10 * n_blocks + bid];
  }

  thrust::copy(h_spine_cnts.begin(), h_spine_cnts.end(), d_spine_cnts);
}

void
cuda_mprts_count_received_v1(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);
  struct cuda_mparticles *cmprts = mprts_cuda->cmprts;

  int n_blocks = cmprts->n_blocks;

  thrust::device_ptr<unsigned int> d_bidx(cmprts->d_bidx);
  thrust::device_ptr<unsigned int> d_spine_cnts(mprts_cuda->d_bnd_spine_cnts);

  thrust::host_vector<unsigned int> h_bidx(cmprts->n_prts);
  thrust::host_vector<unsigned int> h_spine_cnts(1 + n_blocks * (10 + 1));

  thrust::copy(d_bidx, d_bidx + cmprts->n_prts, h_bidx.begin());
  thrust::copy(d_spine_cnts, d_spine_cnts + 1 + n_blocks * (10 + 1), h_spine_cnts.begin());
  for (int n = cmprts->n_prts - mprts_cuda->nr_prts_recv; n < cmprts->n_prts; n++) {
    assert(h_bidx[n] < n_blocks);
    h_spine_cnts[h_bidx[n] * 10 + CUDA_BND_S_NEW]++;
  }
  thrust::copy(h_spine_cnts.begin(), h_spine_cnts.end(), d_spine_cnts);
}

// ======================================================================
// cuda_mprts_scan_scatter_received

static void __global__
mprts_scan_scatter_received(unsigned int nr_recv, unsigned int nr_prts_prev,
			    unsigned int *d_spine_sums, unsigned int *d_alt_bidx,
			    unsigned int *d_bidx, unsigned int *d_ids)
{
  int n = threadIdx.x + THREADS_PER_BLOCK * blockIdx.x;
  if (n >= nr_recv) {
    return;
  }

  n += nr_prts_prev;

  int nn = d_spine_sums[d_bidx[n] * 10 + CUDA_BND_S_NEW] + d_alt_bidx[n];
  d_ids[nn] = n;
}

void
cuda_mprts_scan_scatter_received(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);
  struct cuda_mparticles *cmprts = mprts_cuda->cmprts;

  int nr_recv = mprts_cuda->nr_prts_recv;
  int nr_prts_prev = cmprts->n_prts - nr_recv;

  int dimGrid = (nr_recv + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;

  mprts_scan_scatter_received<<<dimGrid, THREADS_PER_BLOCK>>>
    (nr_recv, nr_prts_prev, mprts_cuda->d_bnd_spine_sums, mprts_cuda->d_alt_bidx,
     cmprts->d_bidx, cmprts->d_id);
}

void
cuda_mprts_scan_scatter_received_gold(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);
  struct cuda_mparticles *cmprts = mprts_cuda->cmprts;

  unsigned int n_blocks = cmprts->n_blocks;

  thrust::device_ptr<unsigned int> d_bidx(cmprts->d_bidx);
  thrust::device_ptr<unsigned int> d_alt_bidx(mprts_cuda->d_alt_bidx);
  thrust::device_ptr<unsigned int> d_id(cmprts->d_id);
  thrust::device_ptr<unsigned int> d_spine_sums(mprts_cuda->d_bnd_spine_sums);

  thrust::host_vector<unsigned int> h_bidx(cmprts->n_prts);
  thrust::host_vector<unsigned int> h_alt_bidx(cmprts->n_prts);
  thrust::host_vector<unsigned int> h_id(cmprts->n_prts);
  thrust::host_vector<unsigned int> h_spine_sums(1 + n_blocks * (10 + 1));

  thrust::copy(d_spine_sums, d_spine_sums + n_blocks * 11, h_spine_sums.begin());
  thrust::copy(d_bidx, d_bidx + cmprts->n_prts, h_bidx.begin());
  thrust::copy(d_alt_bidx, d_alt_bidx + cmprts->n_prts, h_alt_bidx.begin());
  for (int n = cmprts->n_prts - mprts_cuda->nr_prts_recv; n < cmprts->n_prts; n++) {
    int nn = h_spine_sums[h_bidx[n] * 10 + CUDA_BND_S_NEW] + h_alt_bidx[n];
    h_id[nn] = n;
  }
  thrust::copy(h_id.begin(), h_id.end(), d_id);
}

// ======================================================================
// cuda_mprts_spine_reduce

void
cuda_mprts_spine_reduce(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);
  struct cuda_mparticles *cmprts = mprts_cuda->cmprts;

  unsigned int n_blocks = cmprts->n_blocks;
  int *b_mx = mprts_cuda->b_mx;

  thrust::device_ptr<unsigned int> d_spine_cnts(mprts_cuda->d_bnd_spine_cnts);
  thrust::device_ptr<unsigned int> d_spine_sums(mprts_cuda->d_bnd_spine_sums);
  thrust::device_ptr<unsigned int> d_bidx(cmprts->d_bidx);
  thrust::device_ptr<unsigned int> d_off(mprts_cuda->d_off);

  // OPT?
  thrust::fill(d_spine_cnts, d_spine_cnts + 1 + n_blocks * (CUDA_BND_STRIDE + 1), 0);

  const int threads = B40C_RADIXSORT_THREADS;
  if (b_mx[0] == 1 && b_mx[1] == 2 && b_mx[2] == 2) {
    RakingReduction3x<K, V, 0, RADIX_BITS, 0,
		      NopFunctor<K>, 2, 2> <<<n_blocks, threads>>>
      (mprts_cuda->d_bnd_spine_cnts, cmprts->d_bidx, mprts_cuda->d_off, n_blocks);
  } else if (b_mx[0] == 1 && b_mx[1] == 4 && b_mx[2] == 4) {
    RakingReduction3x<K, V, 0, RADIX_BITS, 0,
		      NopFunctor<K>, 4, 4> <<<n_blocks, threads>>>
      (mprts_cuda->d_bnd_spine_cnts, cmprts->d_bidx, mprts_cuda->d_off, n_blocks);
  } else if (b_mx[0] == 1 && b_mx[1] == 8 && b_mx[2] == 8) {
    RakingReduction3x<K, V, 0, RADIX_BITS, 0,
		      NopFunctor<K>, 8, 8> <<<n_blocks, threads>>>
      (mprts_cuda->d_bnd_spine_cnts, cmprts->d_bidx, mprts_cuda->d_off, n_blocks);
  } else if (b_mx[0] == 1 && b_mx[1] == 16 && b_mx[2] == 16) {
    RakingReduction3x<K, V, 0, RADIX_BITS, 0,
		      NopFunctor<K>, 16, 16> <<<n_blocks, threads>>>
      (mprts_cuda->d_bnd_spine_cnts, cmprts->d_bidx, mprts_cuda->d_off, n_blocks);
  } else if (b_mx[0] == 1 && b_mx[1] == 32 && b_mx[2] == 32) {
    RakingReduction3x<K, V, 0, RADIX_BITS, 0,
		      NopFunctor<K>, 32, 32> <<<n_blocks, threads>>>
      (mprts_cuda->d_bnd_spine_cnts, cmprts->d_bidx, mprts_cuda->d_off, n_blocks);
  } else if (b_mx[0] == 1 && b_mx[1] == 64 && b_mx[2] == 64) {
    RakingReduction3x<K, V, 0, RADIX_BITS, 0,
		      NopFunctor<K>, 64, 64> <<<n_blocks, threads>>>
      (mprts_cuda->d_bnd_spine_cnts, cmprts->d_bidx, mprts_cuda->d_off, n_blocks);
  } else if (b_mx[0] == 1 && b_mx[1] == 128 && b_mx[2] == 128) {
    RakingReduction3x<K, V, 0, RADIX_BITS, 0,
                      NopFunctor<K>, 128, 128> <<<n_blocks, threads>>>
      (mprts_cuda->d_bnd_spine_cnts, cmprts->d_bidx, mprts_cuda->d_off, n_blocks);
  } else {
    mprintf("no support for b_mx %d x %d x %d!\n", b_mx[0], b_mx[1], b_mx[2]);
    assert(0);
  }
  cuda_sync_if_enabled();

  thrust::exclusive_scan(d_spine_cnts + n_blocks * 10,
			 d_spine_cnts + n_blocks * 10 + n_blocks + 1,
			 d_spine_sums + n_blocks * 10);
}

void
cuda_mprts_spine_reduce_gold(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);
  struct cuda_mparticles *cmprts = mprts_cuda->cmprts;

  unsigned int n_blocks = cmprts->n_blocks;
  unsigned int n_blocks_per_patch = cmprts->n_blocks_per_patch;
  int *b_mx = mprts_cuda->b_mx;

  thrust::device_ptr<unsigned int> d_spine_cnts(mprts_cuda->d_bnd_spine_cnts);
  thrust::device_ptr<unsigned int> d_spine_sums(mprts_cuda->d_bnd_spine_sums);
  thrust::device_ptr<unsigned int> d_bidx(cmprts->d_bidx);
  thrust::device_ptr<unsigned int> d_off(mprts_cuda->d_off);

  thrust::fill(d_spine_cnts, d_spine_cnts + 1 + n_blocks * (CUDA_BND_STRIDE + 1), 0);

  thrust::host_vector<unsigned int> h_bidx(d_bidx, d_bidx + cmprts->n_prts);
  thrust::host_vector<unsigned int> h_off(d_off, d_off + n_blocks + 1);
  thrust::host_vector<unsigned int> h_spine_cnts(d_spine_cnts, d_spine_cnts + 1 + n_blocks * (CUDA_BND_STRIDE + 1));

  
  for (int p = 0; p < mprts->nr_patches; p++) {
    for (int b = 0; b < n_blocks_per_patch; b++) {
      unsigned int bid = b + p * n_blocks_per_patch;
      for (int n = h_off[bid]; n < h_off[bid+1]; n++) {
	unsigned int key = h_bidx[n];
	if (key < 9) {
	  int dy = key % 3;
	  int dz = key / 3;
	  int by = b % b_mx[1];
	  int bz = b / b_mx[1];
	  unsigned int bby = by + 1 - dy;
	  unsigned int bbz = bz + 1 - dz;
	  unsigned int bb = bbz * b_mx[1] + bby;
	  if (bby < b_mx[1] && bbz < b_mx[2]) {
	    h_spine_cnts[(bb + p * n_blocks_per_patch) * 10 + key]++;
	  } else {
	    assert(0);
	  }
	} else if (key == CUDA_BND_S_OOB) {
	  h_spine_cnts[b_mx[1]*b_mx[2]*mprts->nr_patches * 10 + bid]++;
	}
      }
    }
  }  

  thrust::copy(h_spine_cnts.begin(), h_spine_cnts.end(), d_spine_cnts);
  thrust::exclusive_scan(d_spine_cnts + n_blocks * 10,
			 d_spine_cnts + n_blocks * 10 + n_blocks + 1,
			 d_spine_sums + n_blocks * 10);
}

// ======================================================================
// cuda_mprts_sort_pairs_device

void
cuda_mprts_sort_pairs_device(struct psc_mparticles *mprts)
{
  static int pr_A, pr_B, pr_C, pr_D;
  if (!pr_B) {
    pr_A = prof_register("xchg_cnt_recvd", 1., 0, 0);
    pr_B = prof_register("xchg_top_scan", 1., 0, 0);
    pr_C = prof_register("xchg_ss_recvd", 1., 0, 0);
    pr_D = prof_register("xchg_bottom_scan", 1., 0, 0);
  }

  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);
  struct cuda_mparticles *cmprts = mprts_cuda->cmprts;

  unsigned int n_blocks = cmprts->n_blocks;

  prof_start(pr_A);
  cuda_mprts_count_received(mprts);
  prof_stop(pr_A);

  prof_start(pr_B);
  // FIXME why isn't 10 + 0 enough?
  thrust::device_ptr<unsigned int> d_spine_cnts(mprts_cuda->d_bnd_spine_cnts);
  thrust::device_ptr<unsigned int> d_spine_sums(mprts_cuda->d_bnd_spine_sums);
  thrust::exclusive_scan(d_spine_cnts, d_spine_cnts + 1 + n_blocks * (10 + 1), d_spine_sums);
  prof_stop(pr_B);

  prof_start(pr_C);
  cuda_mprts_scan_scatter_received(mprts);
  prof_stop(pr_C);

  prof_start(pr_D);
  int *b_mx = mprts_cuda->b_mx;
  if (b_mx[0] == 1 && b_mx[1] == 8 && b_mx[2] == 8) {
    ScanScatterDigits3x<K, V, 0, RADIX_BITS, 0,
			NopFunctor<K>,
			NopFunctor<K>,
			8, 8> 
      <<<n_blocks, B40C_RADIXSORT_THREADS>>>
      (mprts_cuda->d_bnd_spine_sums, cmprts->d_bidx,
       cmprts->d_id, mprts_cuda->d_off, n_blocks);
  } else if (b_mx[0] == 1 && b_mx[1] == 16 && b_mx[2] == 16) {
    ScanScatterDigits3x<K, V, 0, RADIX_BITS, 0,
			NopFunctor<K>,
			NopFunctor<K>,
			16, 16> 
      <<<n_blocks, B40C_RADIXSORT_THREADS>>>
      (mprts_cuda->d_bnd_spine_sums, cmprts->d_bidx,
       cmprts->d_id, mprts_cuda->d_off, n_blocks);
  } else if (b_mx[0] == 1 && b_mx[1] == 32 && b_mx[2] == 32) {
    ScanScatterDigits3x<K, V, 0, RADIX_BITS, 0,
			NopFunctor<K>,
			NopFunctor<K>,
			32, 32> 
      <<<n_blocks, B40C_RADIXSORT_THREADS>>>
      (mprts_cuda->d_bnd_spine_sums, cmprts->d_bidx,
       cmprts->d_id, mprts_cuda->d_off, n_blocks);
  } else if (b_mx[0] == 1 && b_mx[1] == 64 && b_mx[2] == 64) {
    ScanScatterDigits3x<K, V, 0, RADIX_BITS, 0,
			NopFunctor<K>,
			NopFunctor<K>,
			64, 64> 
      <<<n_blocks, B40C_RADIXSORT_THREADS>>>
      (mprts_cuda->d_bnd_spine_sums, cmprts->d_bidx,
       cmprts->d_id, mprts_cuda->d_off, n_blocks);
  } else if (b_mx[0] == 1 && b_mx[1] == 128 && b_mx[2] == 128) {
    ScanScatterDigits3x<K, V, 0, RADIX_BITS, 0,
                        NopFunctor<K>,
                        NopFunctor<K>,
                        128, 128>
      <<<n_blocks, B40C_RADIXSORT_THREADS>>>
      (mprts_cuda->d_bnd_spine_sums, cmprts->d_bidx,
       cmprts->d_id, mprts_cuda->d_off, n_blocks);
  } else {
    mprintf("no support for b_mx %d x %d x %d!\n", b_mx[0], b_mx[1], b_mx[2]);
    assert(0);
  }
  cuda_sync_if_enabled();
  prof_stop(pr_D);

  // d_ids now contains the indices to reorder by
}

void
cuda_mprts_sort_pairs_gold(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);
  struct cuda_mparticles *cmprts = mprts_cuda->cmprts;

  unsigned int n_blocks_per_patch = cmprts->n_blocks_per_patch;
  unsigned int n_blocks = cmprts->n_blocks;
  int *b_mx = mprts_cuda->b_mx;

  thrust::device_ptr<unsigned int> d_bidx(cmprts->d_bidx);
  thrust::device_ptr<unsigned int> d_id(cmprts->d_id);
  thrust::device_ptr<unsigned int> d_off(mprts_cuda->d_off);
  thrust::device_ptr<unsigned int> d_spine_cnts(mprts_cuda->d_bnd_spine_cnts);
  thrust::device_ptr<unsigned int> d_spine_sums(mprts_cuda->d_bnd_spine_sums);

  thrust::host_vector<unsigned int> h_bidx(d_bidx, d_bidx + cmprts->n_prts);
  thrust::host_vector<unsigned int> h_id(cmprts->n_prts);
  thrust::host_vector<unsigned int> h_off(d_off, d_off + n_blocks + 1);
  thrust::host_vector<unsigned int> h_spine_cnts(d_spine_cnts, d_spine_cnts + 1 + n_blocks * (10 + 1));

  thrust::host_vector<unsigned int> h_spine_sums(1 + n_blocks * (10 + 1));

  for (int n = cmprts->n_prts - mprts_cuda->nr_prts_recv; n < cmprts->n_prts; n++) {
    assert(h_bidx[n] < n_blocks);
    h_spine_cnts[h_bidx[n] * 10 + CUDA_BND_S_NEW]++;
  }

  thrust::exclusive_scan(h_spine_cnts.begin(), h_spine_cnts.end(), h_spine_sums.begin());
  thrust::copy(h_spine_sums.begin(), h_spine_sums.end(), d_spine_sums);

  for (int bid = 0; bid < n_blocks; bid++) {
    int b = bid % n_blocks_per_patch;
    int p = bid / n_blocks_per_patch;
    for (int n = h_off[bid]; n < h_off[bid+1]; n++) {
      unsigned int key = h_bidx[n];
      if (key < 9) {
	int dy = key % 3;
	int dz = key / 3;
	int by = b % b_mx[1];
	int bz = b / b_mx[1];
	unsigned int bby = by + 1 - dy;
	unsigned int bbz = bz + 1 - dz;
	assert(bby < b_mx[1] && bbz < b_mx[2]);
	unsigned int bb = bbz * b_mx[1] + bby;
	int nn = h_spine_sums[(bb + p * n_blocks_per_patch) * 10 + key]++;
	h_id[nn] = n;
      } else { // OOB
	assert(0);
      }
    }
  }
  for (int n = cmprts->n_prts - mprts_cuda->nr_prts_recv; n < cmprts->n_prts; n++) {
      int nn = h_spine_sums[h_bidx[n] * 10 + CUDA_BND_S_NEW]++;
      h_id[nn] = n;
  }

  thrust::copy(h_id.begin(), h_id.end(), d_id);
  // d_ids now contains the indices to reorder by
}

// ======================================================================
// cuda_mprts_update_offsets

__global__ static void
mprts_update_offsets(int nr_total_blocks, unsigned int *d_off, unsigned int *d_spine_sums)
{
  int bid = threadIdx.x + THREADS_PER_BLOCK * blockIdx.x;
  
  if (bid <= nr_total_blocks) {
    d_off[bid] = d_spine_sums[bid * CUDA_BND_STRIDE + 0];
  }
}

void
cuda_mprts_update_offsets(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);
  struct cuda_mparticles *cmprts = mprts_cuda->cmprts;

  unsigned int n_blocks = cmprts->n_blocks;
  int dimGrid = (n_blocks + 1 + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;

  mprts_update_offsets<<<dimGrid, THREADS_PER_BLOCK>>>
    (n_blocks, mprts_cuda->d_off, mprts_cuda->d_bnd_spine_sums);
}

void
cuda_mprts_update_offsets_gold(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);
  struct cuda_mparticles *cmprts = mprts_cuda->cmprts;

  unsigned int n_blocks = cmprts->n_blocks;

  thrust::device_ptr<unsigned int> d_spine_sums(mprts_cuda->d_bnd_spine_sums);
  thrust::device_ptr<unsigned int> d_off(mprts_cuda->d_off);

  thrust::host_vector<unsigned int> h_spine_sums(d_spine_sums, d_spine_sums + 1 + n_blocks * (10 + 1));
  thrust::host_vector<unsigned int> h_off(n_blocks + 1);

  for (int bid = 0; bid <= n_blocks; bid++) {
    h_off[bid] = h_spine_sums[bid * 10];
  }

  thrust::copy(h_off.begin(), h_off.end(), d_off);
}

