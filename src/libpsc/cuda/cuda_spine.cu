
#include "psc_cuda.h"
#include "particles_cuda.h"

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/scan.h>

#include <b40c/radixsort_reduction_kernel.h>
#include <b40c/radixsort_spine_kernel.h>
#include <b40c/radixsort_scanscatter_kernel3.h>

using namespace b40c_thrust;

typedef unsigned int K;
typedef unsigned int V;

static const int RADIX_BITS = 4;

#define NBLOCKS_X 1
#define NBLOCKS_Y 8
#define NBLOCKS_Z 8

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

  unsigned int nr_total_blocks = mprts_cuda->nr_total_blocks;

  int dimGrid  = nr_total_blocks;
  
  mprts_bidx_to_key<<<dimGrid, THREADS_PER_BLOCK>>>
    (mprts_cuda->nr_total_blocks, mprts_cuda->d_off, mprts_cuda->d_bidx);
}

void
cuda_mprts_bidx_to_key_gold(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);

  unsigned int nr_total_blocks = mprts_cuda->nr_total_blocks;
  unsigned int nr_blocks = mprts_cuda->nr_blocks;

  thrust::device_ptr<unsigned int> d_bidx(mprts_cuda->d_bidx);
  thrust::device_ptr<unsigned int> d_off(mprts_cuda->d_off);

  thrust::host_vector<unsigned int> h_bidx(d_bidx, d_bidx + mprts_cuda->nr_prts);
  thrust::host_vector<unsigned int> h_off(d_off, d_off + nr_total_blocks + 1);

  for (int bid = 0; bid < nr_total_blocks; bid++) {
    int p = bid / nr_blocks;
    for (int n = h_off[bid]; n < h_off[bid+1]; n++) {
      assert((h_bidx[n] >= p * nr_blocks && h_bidx[n] < (p+1) * nr_blocks) ||
	     (h_bidx[n] == nr_total_blocks));
      int bidx;
      if (h_bidx[n] == mprts_cuda->nr_total_blocks) {
	bidx = CUDA_BND_S_OOB;
      } else {
	int b_diff = bid - h_bidx[n] + NBLOCKS_Y + 1;
	int d1 = b_diff % NBLOCKS_Y;
	int d2 = b_diff / NBLOCKS_Y;
	bidx = d2 * 3 + d1;
      }
      
      h_bidx[n] = bidx;
    }
  }

  thrust::copy(h_bidx.begin(), h_bidx.end(), d_bidx);
}

// ======================================================================
// cuda_mprts_spine_reduce

void
cuda_mprts_spine_reduce(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);

  unsigned int nr_total_blocks = mprts_cuda->nr_total_blocks;
  unsigned int nr_blocks = mprts_cuda->nr_blocks;
  assert(nr_blocks == NBLOCKS_Y * NBLOCKS_Z);

  thrust::device_ptr<unsigned int> d_spine_cnts(mprts_cuda->d_bnd_spine_cnts);
  thrust::device_ptr<unsigned int> d_spine_sums(mprts_cuda->d_bnd_spine_sums);
  thrust::device_ptr<unsigned int> d_bidx(mprts_cuda->d_bidx);
  thrust::device_ptr<unsigned int> d_off(mprts_cuda->d_off);

  // OPT?
  thrust::fill(d_spine_cnts, d_spine_cnts + 1 + nr_total_blocks * (CUDA_BND_STRIDE + 1), 0);

#if 0
  thrust::host_vector<unsigned int> h_bidx(d_bidx, d_bidx + mprts_cuda->nr_prts);
  thrust::host_vector<unsigned int> h_off(d_off, d_off + nr_total_blocks + 1);
  thrust::host_vector<unsigned int> h_spine_cnts(d_spine_cnts, d_spine_cnts + 1 + nr_total_blocks * (CUDA_BND_STRIDE + 1));
#endif

  const int threads = B40C_RADIXSORT_THREADS;
  RakingReduction3x<K, V, 0, RADIX_BITS, 0,
		    NopFunctor<K>, NBLOCKS_Y, NBLOCKS_Z> <<<nr_total_blocks, threads>>>
    (mprts_cuda->d_bnd_spine_cnts, mprts_cuda->d_bidx, mprts_cuda->d_off, nr_total_blocks);
  cuda_sync_if_enabled();

#if 0
  for (int p = 0; p < mprts->nr_patches; p++) {
    for (int b = 0; b < nr_blocks; b++) {
      unsigned int bid = b + p * nr_blocks;
      for (int n = h_off[bid]; n < h_off[bid+1]; n++) {
	unsigned int key = h_bidx[n];
	if (key < 9) {
	  int dy = key % 3;
	  int dz = key / 3;
	  int by = b % NBLOCKS_Y;
	  int bz = b / NBLOCKS_Y;
	  unsigned int bby = by + 1 - dy;
	  unsigned int bbz = bz + 1 - dz;
	  unsigned int bb = bbz * NBLOCKS_Y + bby;
	  if (bby < NBLOCKS_Y && bbz < NBLOCKS_Z) {
	    h_spine_cnts[(bb + p * nr_blocks) * 10 + key]++;
	  } else {
	    assert(0);
	  }
	} else if (key == CUDA_BND_S_OOB) {
	  h_spine_cnts[NBLOCKS_Y*NBLOCKS_Z*mprts->nr_patches * 10 + bid]++;
	}
      }
    }
  }  
  //thrust::copy(h_spine_cnts.begin(), h_spine_cnts.end(), d_spine_cnts);
  thrust::host_vector<unsigned int> h_spine_cnts2(nr_total_blocks * 10 + nr_total_blocks + 1);
  thrust::copy(d_spine_cnts, d_spine_cnts + nr_total_blocks * 11 + 1, h_spine_cnts2.begin());
  for (int i = 0; i < nr_total_blocks * 11 + 1; i++) {
    if (h_spine_cnts[i] != h_spine_cnts2[i]) {
      mprintf("i %d: %d // %d\n", i, h_spine_cnts[i], h_spine_cnts2[i]);
    }
  }
  assert(0);
#endif

  thrust::exclusive_scan(d_spine_cnts + nr_total_blocks * 10,
			 d_spine_cnts + nr_total_blocks * 10 + nr_total_blocks + 1,
			 d_spine_sums + nr_total_blocks * 10);
}

void
cuda_mprts_spine_reduce_gold(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);

  unsigned int nr_total_blocks = mprts_cuda->nr_total_blocks;
  unsigned int nr_blocks = mprts_cuda->nr_blocks;
  assert(nr_blocks == NBLOCKS_Y * NBLOCKS_Z);

  thrust::device_ptr<unsigned int> d_spine_cnts(mprts_cuda->d_bnd_spine_cnts);
  thrust::device_ptr<unsigned int> d_spine_sums(mprts_cuda->d_bnd_spine_sums);
  thrust::device_ptr<unsigned int> d_bidx(mprts_cuda->d_bidx);
  thrust::device_ptr<unsigned int> d_off(mprts_cuda->d_off);

  thrust::fill(d_spine_cnts, d_spine_cnts + 1 + nr_total_blocks * (CUDA_BND_STRIDE + 1), 0);

  thrust::host_vector<unsigned int> h_bidx(d_bidx, d_bidx + mprts_cuda->nr_prts);
  thrust::host_vector<unsigned int> h_off(d_off, d_off + nr_total_blocks + 1);
  thrust::host_vector<unsigned int> h_spine_cnts(d_spine_cnts, d_spine_cnts + 1 + nr_total_blocks * (CUDA_BND_STRIDE + 1));

  
  for (int p = 0; p < mprts->nr_patches; p++) {
    for (int b = 0; b < nr_blocks; b++) {
      unsigned int bid = b + p * nr_blocks;
      for (int n = h_off[bid]; n < h_off[bid+1]; n++) {
	unsigned int key = h_bidx[n];
	if (key < 9) {
	  int dy = key % 3;
	  int dz = key / 3;
	  int by = b % NBLOCKS_Y;
	  int bz = b / NBLOCKS_Y;
	  unsigned int bby = by + 1 - dy;
	  unsigned int bbz = bz + 1 - dz;
	  unsigned int bb = bbz * NBLOCKS_Y + bby;
	  if (bby < NBLOCKS_Y && bbz < NBLOCKS_Z) {
	    h_spine_cnts[(bb + p * nr_blocks) * 10 + key]++;
	  } else {
	    assert(0);
	  }
	} else if (key == CUDA_BND_S_OOB) {
	  h_spine_cnts[NBLOCKS_Y*NBLOCKS_Z*mprts->nr_patches * 10 + bid]++;
	}
      }
    }
  }  

  thrust::copy(h_spine_cnts.begin(), h_spine_cnts.end(), d_spine_cnts);
  thrust::exclusive_scan(d_spine_cnts + nr_total_blocks * 10,
			 d_spine_cnts + nr_total_blocks * 10 + nr_total_blocks + 1,
			 d_spine_sums + nr_total_blocks * 10);
}

void
cuda_mprts_sort_pairs_device(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);

  int nr_blocks = mprts_cuda->nr_blocks;
  int nr_total_blocks = mprts_cuda->nr_total_blocks;

  thrust::device_ptr<unsigned int> d_bidx(mprts_cuda->d_bidx);
  thrust::device_ptr<unsigned int> d_ids(mprts_cuda->d_ids);
  thrust::device_ptr<unsigned int> d_off(mprts_cuda->d_off);
  thrust::device_ptr<unsigned int> d_spine_cnts(mprts_cuda->d_bnd_spine_cnts);
  thrust::device_ptr<unsigned int> d_spine_sums(mprts_cuda->d_bnd_spine_sums);

  thrust::host_vector<unsigned int> h_bidx(d_bidx, d_bidx + mprts_cuda->nr_prts);
  thrust::host_vector<unsigned int> h_ids(mprts_cuda->nr_prts);
  thrust::host_vector<unsigned int> h_off(d_off, d_off + nr_total_blocks + 1);
  thrust::host_vector<unsigned int> h_spine_cnts(d_spine_cnts, d_spine_cnts + 1 + nr_total_blocks * (10 + 1));

  thrust::host_vector<unsigned int> h_spine_sums(1 + nr_total_blocks * (10 + 1));

  for (int n = mprts_cuda->nr_prts - mprts_cuda->nr_prts_recv; n < mprts_cuda->nr_prts; n++) {
    assert(h_bidx[n] < mprts_cuda->nr_total_blocks);
    h_spine_cnts[h_bidx[n] * 10 + CUDA_BND_S_NEW]++;
  }

  thrust::exclusive_scan(h_spine_cnts.begin(), h_spine_cnts.end(), h_spine_sums.begin());
  thrust::copy(h_spine_sums.begin(), h_spine_sums.end(), d_spine_sums);

  for (int bid = 0; bid < nr_total_blocks; bid++) {
    int b = bid % nr_blocks;
    int p = bid / nr_blocks;
    for (int n = h_off[bid]; n < h_off[bid+1]; n++) {
      unsigned int key = h_bidx[n];
      if (key < 9) {
	int dy = key % 3;
	int dz = key / 3;
	int by = b % NBLOCKS_Y;
	int bz = b / NBLOCKS_Y;
	unsigned int bby = by + 1 - dy;
	unsigned int bbz = bz + 1 - dz;
	assert(bby < NBLOCKS_Y && bbz < NBLOCKS_Z);
	unsigned int bb = bbz * NBLOCKS_Y + bby;
	int nn = h_spine_sums[(bb + p * nr_blocks) * 10 + key]++;
	h_ids[nn] = n;
      } else { // OOB
      }
    }
  }
  for (int n = mprts_cuda->nr_prts - mprts_cuda->nr_prts_recv; n < mprts_cuda->nr_prts; n++) {
      int nn = h_spine_sums[h_bidx[n] * 10 + CUDA_BND_S_NEW]++;
      h_ids[nn] = n;
  }

  thrust::copy(h_ids.begin(), h_ids.end(), d_ids);
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

  int nr_total_blocks = mprts_cuda->nr_total_blocks;
  int dimGrid = (nr_total_blocks + 1 + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;

  mprts_update_offsets<<<dimGrid, THREADS_PER_BLOCK>>>
    (mprts_cuda->nr_total_blocks, mprts_cuda->d_off, mprts_cuda->d_bnd_spine_sums);
}

void
cuda_mprts_update_offsets_gold(struct psc_mparticles *mprts)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);

  int nr_total_blocks = mprts_cuda->nr_total_blocks;

  thrust::device_ptr<unsigned int> d_spine_sums(mprts_cuda->d_bnd_spine_sums);
  thrust::device_ptr<unsigned int> d_off(mprts_cuda->d_off);

  thrust::host_vector<unsigned int> h_spine_sums(d_spine_sums, d_spine_sums + 1 + nr_total_blocks * (10 + 1));
  thrust::host_vector<unsigned int> h_off(nr_total_blocks + 1);

  for (int bid = 0; bid <= nr_total_blocks; bid++) {
    h_off[bid] = h_spine_sums[bid * 10];
  }

  thrust::copy(h_off.begin(), h_off.end(), d_off);
}

