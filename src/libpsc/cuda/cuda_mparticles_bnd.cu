
#include "cuda_mparticles.h"
#include "cuda_bits.h"

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/scan.h>

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

#include <cstdio>
#include <cassert>

// ----------------------------------------------------------------------
// cuda_mparticles_bnd_setup

void
cuda_mparticles_bnd_setup(struct cuda_mparticles *cmprts)
{
  cudaError_t ierr;

  cmprts->bnd.h_bnd_cnt = new unsigned int[cmprts->n_blocks];

  ierr = cudaMalloc((void **) &cmprts->bnd.d_bnd_spine_cnts,
		    (1 + cmprts->n_blocks * (CUDA_BND_STRIDE + 1)) * sizeof(unsigned int)); cudaCheck(ierr);
  ierr = cudaMalloc((void **) &cmprts->bnd.d_bnd_spine_sums,
		    (1 + cmprts->n_blocks * (CUDA_BND_STRIDE + 1)) * sizeof(unsigned int)); cudaCheck(ierr);
}  

// ----------------------------------------------------------------------
// cuda_mparticles_bnd_free_particle_mem

void
cuda_mparticles_bnd_free_particle_mem(struct cuda_mparticles *cmprts)
{
  cudaError_t ierr;

  ierr = cudaFree(cmprts->bnd.d_alt_bidx); cudaCheck(ierr);
  ierr = cudaFree(cmprts->bnd.d_sums); cudaCheck(ierr);
}

// ----------------------------------------------------------------------
// cuda_mparticles_bnd_destroy

void
cuda_mparticles_bnd_destroy(struct cuda_mparticles *cmprts)
{
  cudaError_t ierr;

  delete[] cmprts->bnd.h_bnd_cnt;

  ierr = cudaFree(cmprts->bnd.d_bnd_spine_cnts); cudaCheck(ierr);
  ierr = cudaFree(cmprts->bnd.d_bnd_spine_sums); cudaCheck(ierr);
}

// ----------------------------------------------------------------------
// cuda_mparticles_bnd_reserve_all

void
cuda_mparticles_bnd_reserve_all(struct cuda_mparticles *cmprts)
{
  cudaError_t ierr;

  int n_alloced = cmprts->n_alloced;
  ierr = cudaMalloc((void **) &cmprts->bnd.d_alt_bidx, n_alloced * sizeof(unsigned int)); cudaCheck(ierr);
  ierr = cudaMalloc((void **) &cmprts->bnd.d_sums, n_alloced * sizeof(unsigned int)); cudaCheck(ierr);
}

// ----------------------------------------------------------------------
// cuda_mparticles_spine_reduce

void
cuda_mparticles_spine_reduce(struct cuda_mparticles *cmprts)
{
  unsigned int n_blocks = cmprts->n_blocks;
  int *b_mx = cmprts->b_mx;

  thrust::device_ptr<unsigned int> d_spine_cnts(cmprts->bnd.d_bnd_spine_cnts);
  thrust::device_ptr<unsigned int> d_spine_sums(cmprts->bnd.d_bnd_spine_sums);
  thrust::device_ptr<unsigned int> d_bidx(cmprts->d_bidx);
  thrust::device_ptr<unsigned int> d_off(cmprts->d_off);

  // OPT?
  thrust::fill(d_spine_cnts, d_spine_cnts + 1 + n_blocks * (CUDA_BND_STRIDE + 1), 0);

  const int threads = B40C_RADIXSORT_THREADS;
  if (b_mx[0] == 1 && b_mx[1] == 2 && b_mx[2] == 2) {
    RakingReduction3x<K, V, 0, RADIX_BITS, 0,
		      NopFunctor<K>, 2, 2> <<<n_blocks, threads>>>
      (cmprts->bnd.d_bnd_spine_cnts, cmprts->d_bidx, cmprts->d_off, n_blocks);
  } else if (b_mx[0] == 1 && b_mx[1] == 4 && b_mx[2] == 4) {
    RakingReduction3x<K, V, 0, RADIX_BITS, 0,
		      NopFunctor<K>, 4, 4> <<<n_blocks, threads>>>
      (cmprts->bnd.d_bnd_spine_cnts, cmprts->d_bidx, cmprts->d_off, n_blocks);
  } else if (b_mx[0] == 1 && b_mx[1] == 8 && b_mx[2] == 8) {
    RakingReduction3x<K, V, 0, RADIX_BITS, 0,
		      NopFunctor<K>, 8, 8> <<<n_blocks, threads>>>
      (cmprts->bnd.d_bnd_spine_cnts, cmprts->d_bidx, cmprts->d_off, n_blocks);
  } else if (b_mx[0] == 1 && b_mx[1] == 16 && b_mx[2] == 16) {
    RakingReduction3x<K, V, 0, RADIX_BITS, 0,
		      NopFunctor<K>, 16, 16> <<<n_blocks, threads>>>
      (cmprts->bnd.d_bnd_spine_cnts, cmprts->d_bidx, cmprts->d_off, n_blocks);
  } else if (b_mx[0] == 1 && b_mx[1] == 32 && b_mx[2] == 32) {
    RakingReduction3x<K, V, 0, RADIX_BITS, 0,
		      NopFunctor<K>, 32, 32> <<<n_blocks, threads>>>
      (cmprts->bnd.d_bnd_spine_cnts, cmprts->d_bidx, cmprts->d_off, n_blocks);
  } else if (b_mx[0] == 1 && b_mx[1] == 64 && b_mx[2] == 64) {
    RakingReduction3x<K, V, 0, RADIX_BITS, 0,
		      NopFunctor<K>, 64, 64> <<<n_blocks, threads>>>
      (cmprts->bnd.d_bnd_spine_cnts, cmprts->d_bidx, cmprts->d_off, n_blocks);
  } else if (b_mx[0] == 1 && b_mx[1] == 128 && b_mx[2] == 128) {
    RakingReduction3x<K, V, 0, RADIX_BITS, 0,
                      NopFunctor<K>, 128, 128> <<<n_blocks, threads>>>
      (cmprts->bnd.d_bnd_spine_cnts, cmprts->d_bidx, cmprts->d_off, n_blocks);
  } else {
    printf("no support for b_mx %d x %d x %d!\n", b_mx[0], b_mx[1], b_mx[2]);
    assert(0);
  }
  cuda_sync_if_enabled();

  thrust::exclusive_scan(d_spine_cnts + n_blocks * 10,
			 d_spine_cnts + n_blocks * 10 + n_blocks + 1,
			 d_spine_sums + n_blocks * 10);
}

// ----------------------------------------------------------------------
// cuda_mprts_spine_reduce_gold

void
cuda_mparticles_spine_reduce_gold(struct cuda_mparticles *cmprts)
{
  unsigned int n_blocks = cmprts->n_blocks;
  unsigned int n_blocks_per_patch = cmprts->n_blocks_per_patch;
  int *b_mx = cmprts->b_mx;

  thrust::device_ptr<unsigned int> d_spine_cnts(cmprts->bnd.d_bnd_spine_cnts);
  thrust::device_ptr<unsigned int> d_spine_sums(cmprts->bnd.d_bnd_spine_sums);
  thrust::device_ptr<unsigned int> d_bidx(cmprts->d_bidx);
  thrust::device_ptr<unsigned int> d_off(cmprts->d_off);

  thrust::fill(d_spine_cnts, d_spine_cnts + 1 + n_blocks * (CUDA_BND_STRIDE + 1), 0);

  thrust::host_vector<unsigned int> h_bidx(d_bidx, d_bidx + cmprts->n_prts);
  thrust::host_vector<unsigned int> h_off(d_off, d_off + n_blocks + 1);
  thrust::host_vector<unsigned int> h_spine_cnts(d_spine_cnts, d_spine_cnts + 1 + n_blocks * (CUDA_BND_STRIDE + 1));

  
  for (int p = 0; p < cmprts->n_patches; p++) {
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
	  h_spine_cnts[b_mx[1]*b_mx[2]*cmprts->n_patches * 10 + bid]++;
	}
      }
    }
  }  

  thrust::copy(h_spine_cnts.begin(), h_spine_cnts.end(), d_spine_cnts);
  thrust::exclusive_scan(d_spine_cnts + n_blocks * 10,
			 d_spine_cnts + n_blocks * 10 + n_blocks + 1,
			 d_spine_sums + n_blocks * 10);
}

// ----------------------------------------------------------------------
// cuda_mparticles_find_n_send

void
cuda_mparticles_find_n_send(struct cuda_mparticles *cmprts)
{
  unsigned int n_blocks = cmprts->n_blocks;

  thrust::device_ptr<unsigned int> d_spine_sums(cmprts->bnd.d_bnd_spine_sums);
  thrust::host_vector<unsigned int> h_spine_sums(n_blocks + 1);

  thrust::copy(d_spine_sums + n_blocks * 10,
	       d_spine_sums + n_blocks * 11 + 1,
	       h_spine_sums.begin());

  unsigned int off = 0;
  for (int p = 0; p < cmprts->n_patches; p++) {
    unsigned int n_send = h_spine_sums[(p + 1) * cmprts->n_blocks_per_patch];
    cmprts->bnd.bpatch[p].n_send = n_send - off;
    off = n_send;
  }
  cmprts->bnd.n_prts_send = off;
}

// ----------------------------------------------------------------------
// cuda_mparticles_copy_from_dev

void
cuda_mparticles_copy_from_dev(struct cuda_mparticles *cmprts)
{
  cudaError_t ierr;

  if (cmprts->n_patches == 0) {
    return;
  }

  cmprts->bnd.h_bnd_xi4 = new float4[cmprts->bnd.n_prts_send];
  cmprts->bnd.h_bnd_pxi4 = new float4[cmprts->bnd.n_prts_send];

  assert(cmprts->n_prts + cmprts->bnd.n_prts_send < cmprts->n_alloced);

  ierr = cudaMemcpy(cmprts->bnd.h_bnd_xi4, cmprts->d_xi4 + cmprts->n_prts,
		    cmprts->bnd.n_prts_send * sizeof(float4), cudaMemcpyDeviceToHost); cudaCheck(ierr);
  ierr = cudaMemcpy(cmprts->bnd.h_bnd_pxi4, cmprts->d_pxi4 + cmprts->n_prts,
		    cmprts->bnd.n_prts_send * sizeof(float4), cudaMemcpyDeviceToHost); cudaCheck(ierr);
}

// ----------------------------------------------------------------------
// cuda_mparticles_convert_from_cuda

void
cuda_mparticles_convert_from_cuda(struct cuda_mparticles *cmprts)
{
  if (cmprts->n_patches == 0) {
    return;
  }

  float4 *bnd_xi4 = cmprts->bnd.h_bnd_xi4;
  float4 *bnd_pxi4 = cmprts->bnd.h_bnd_pxi4;
  for (int p = 0; p < cmprts->n_patches; p++) {
    for (int n = 0; n < cmprts->bnd.bpatch[p].n_send; n++) {
      particle_cuda_t *prt = &cmprts->bnd.bpatch[p].buf.m_data[n];
      prt->xi      = bnd_xi4[n].x;
      prt->yi      = bnd_xi4[n].y;
      prt->zi      = bnd_xi4[n].z;
      prt->kind    = cuda_float_as_int(bnd_xi4[n].w);
      prt->pxi     = bnd_pxi4[n].x;
      prt->pyi     = bnd_pxi4[n].y;
      prt->pzi     = bnd_pxi4[n].z;
      prt->qni_wni = bnd_pxi4[n].w;
    }
    bnd_xi4 += cmprts->bnd.bpatch[p].n_send;
    bnd_pxi4 += cmprts->bnd.bpatch[p].n_send;
  }
  delete[] cmprts->bnd.h_bnd_xi4;
  delete[] cmprts->bnd.h_bnd_pxi4;
}

// ----------------------------------------------------------------------
// cuda_mparticles_copy_to_dev

void
cuda_mparticles_copy_to_dev(struct cuda_mparticles *cmprts)
{
  cudaError_t ierr;

  float4 *d_xi4 = cmprts->d_xi4;
  float4 *d_pxi4 = cmprts->d_pxi4;

  unsigned int nr_recv = 0;
  for (int p = 0; p < cmprts->n_patches; p++) {
    nr_recv += cmprts->bnd.bpatch[p].n_recv;
  }
  assert(cmprts->n_prts + nr_recv <= cmprts->n_alloced);

  ierr = cudaMemcpy(d_xi4 + cmprts->n_prts, cmprts->bnd.h_bnd_xi4,
		    nr_recv * sizeof(*d_xi4), cudaMemcpyHostToDevice); cudaCheck(ierr);
  ierr = cudaMemcpy(d_pxi4 + cmprts->n_prts, cmprts->bnd.h_bnd_pxi4,
		    nr_recv * sizeof(*d_pxi4), cudaMemcpyHostToDevice); cudaCheck(ierr);

  free(cmprts->bnd.h_bnd_xi4);
  free(cmprts->bnd.h_bnd_pxi4);

  cmprts->bnd.n_prts_recv = nr_recv;
  cmprts->n_prts += nr_recv;
}

