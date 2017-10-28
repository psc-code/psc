
#undef _GLIBCXX_USE_INT128

#include "cuda_mparticles.h"

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/scan.h>

#include "psc_cuda.h"
#include <mrc_profile.h>

#include <b40c/radixsort_reduction_kernel.h>
#include <b40c/radixsort_scanscatter_kernel3.h>

using namespace b40c_thrust;

typedef unsigned int K;
typedef unsigned int V;

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
  struct cuda_mparticles *cmprts = psc_mparticles_cuda(mprts)->cmprts;

  unsigned int n_blocks = cmprts->n_blocks;
  
  int *b_mx = cmprts->b_mx;
  if (b_mx[0] == 1 && b_mx[1] == 8 && b_mx[2] == 8) {
    mprts_bidx_to_key<1, 8, 8><<<n_blocks, THREADS_PER_BLOCK>>>
      (n_blocks, cmprts->d_off, cmprts->d_bidx);
  } else {
    mprintf("no support for b_mx %d x %d x %d!\n", b_mx[0], b_mx[1], b_mx[2]);
    assert(0);
  }
  cuda_sync_if_enabled();
}

void
cuda_mprts_bidx_to_key_gold(struct psc_mparticles *mprts)
{
  struct cuda_mparticles *cmprts = psc_mparticles_cuda(mprts)->cmprts;

  unsigned int n_blocks = cmprts->n_blocks;
  unsigned int n_blocks_per_patch = cmprts->n_blocks_per_patch;

  thrust::device_ptr<unsigned int> d_bidx(cmprts->d_bidx);
  thrust::device_ptr<unsigned int> d_off(cmprts->d_off);

  thrust::host_vector<unsigned int> h_bidx(d_bidx, d_bidx + cmprts->n_prts);
  thrust::host_vector<unsigned int> h_off(d_off, d_off + n_blocks + 1);

  int *b_mx = cmprts->b_mx;

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


