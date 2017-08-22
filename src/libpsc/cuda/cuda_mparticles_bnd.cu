
#include "cuda_mparticles.h"
#include "cuda_bits.h"

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/scan.h>

#include <b40c/radixsort_reduction_kernel.h>
#include <b40c/radixsort_scanscatter_kernel3.h>

#include "psc_particle_buf_cuda.h"

#include <mrc_profile.h>

#define THREADS_PER_BLOCK 256

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
    psc_particle_cuda_buf_t *buf = &cmprts->bnd.bpatch[p].buf;
    unsigned int n_send = cmprts->bnd.bpatch[p].n_send;
    psc_particle_cuda_buf_ctor(buf);
    psc_particle_cuda_buf_reserve(buf, n_send);
    psc_particle_cuda_buf_resize(buf, n_send);

    for (int n = 0; n < n_send; n++) {
      particle_cuda_t *prt = psc_particle_cuda_buf_at_ptr(buf, n);
      prt->xi      = bnd_xi4[n].x;
      prt->yi      = bnd_xi4[n].y;
      prt->zi      = bnd_xi4[n].z;
      prt->kind    = cuda_float_as_int(bnd_xi4[n].w);
      prt->pxi     = bnd_pxi4[n].x;
      prt->pyi     = bnd_pxi4[n].y;
      prt->pzi     = bnd_pxi4[n].z;
      prt->qni_wni = bnd_pxi4[n].w;
    }
    bnd_xi4 += n_send;
    bnd_pxi4 += n_send;
  }

  delete[] cmprts->bnd.h_bnd_xi4;
  delete[] cmprts->bnd.h_bnd_pxi4;
}

// ----------------------------------------------------------------------
// cuda_mparticles_convert_to_cuda

void
cuda_mparticles_convert_to_cuda(struct cuda_mparticles *cmprts)
{
  unsigned int n_recv = 0;
  for (int p = 0; p < cmprts->n_patches; p++) {
    n_recv += psc_particle_cuda_buf_size(&cmprts->bnd.bpatch[p].buf);
  }

  cmprts->bnd.h_bnd_xi4  = (float4 *) malloc(n_recv * sizeof(*cmprts->bnd.h_bnd_xi4));
  cmprts->bnd.h_bnd_pxi4 = (float4 *) malloc(n_recv * sizeof(*cmprts->bnd.h_bnd_pxi4));
  cmprts->bnd.h_bnd_idx  = (unsigned int *) malloc(n_recv * sizeof(*cmprts->bnd.h_bnd_idx));
  cmprts->bnd.h_bnd_off  = (unsigned int *) malloc(n_recv * sizeof(*cmprts->bnd.h_bnd_off));

  memset(cmprts->bnd.h_bnd_cnt, 0,
	 cmprts->n_blocks * sizeof(*cmprts->bnd.h_bnd_cnt));

  unsigned int off = 0;
  for (int p = 0; p < cmprts->n_patches; p++) {
    int n_recv = psc_particle_cuda_buf_size(&cmprts->bnd.bpatch[p].buf);
    cmprts->bnd.bpatch[p].n_recv = n_recv;
    
    float4 *h_bnd_xi4 = cmprts->bnd.h_bnd_xi4 + off;
    float4 *h_bnd_pxi4 = cmprts->bnd.h_bnd_pxi4 + off;
    unsigned int *h_bnd_idx = cmprts->bnd.h_bnd_idx + off;
    unsigned int *h_bnd_off = cmprts->bnd.h_bnd_off + off;
    for (int n = 0; n < n_recv; n++) {
      particle_cuda_t *prt = &cmprts->bnd.bpatch[p].buf.m_data[n];
      h_bnd_xi4[n].x  = prt->xi;
      h_bnd_xi4[n].y  = prt->yi;
      h_bnd_xi4[n].z  = prt->zi;
      h_bnd_xi4[n].w  = cuda_int_as_float(prt->kind);
      h_bnd_pxi4[n].x = prt->pxi;
      h_bnd_pxi4[n].y = prt->pyi;
      h_bnd_pxi4[n].z = prt->pzi;
      h_bnd_pxi4[n].w = prt->qni_wni;

      int b_pos[3];
      for (int d = 0; d < 3; d++) {
	float *xi = &h_bnd_xi4[n].x;
	b_pos[d] = particle_cuda_real_fint(xi[d] * cmprts->b_dxi[d]);
	if (b_pos[d] < 0 || b_pos[d] >= cmprts->b_mx[d]) {
	  printf("!!! xi %g %g %g\n", xi[0], xi[1], xi[2]);
	  printf("!!! d %d xi4[n] %g biy %d // %d\n",
		 d, xi[d], b_pos[d], cmprts->b_mx[d]);
	  if (b_pos[d] < 0) {
	    xi[d] = 0.f;
	  } else {
	    xi[d] *= (1. - 1e-6);
	  }
	}
	b_pos[d] = particle_cuda_real_fint(xi[d] * cmprts->b_dxi[d]);
	assert(b_pos[d] >= 0 && b_pos[d] < cmprts->b_mx[d]);
      }
      unsigned int b = (b_pos[2] * cmprts->b_mx[1] + b_pos[1]) * cmprts->b_mx[0] + b_pos[0];
      assert(b < cmprts->n_blocks_per_patch);
      b += p * cmprts->n_blocks_per_patch;
      h_bnd_idx[n] = b;
      h_bnd_off[n] = cmprts->bnd.h_bnd_cnt[b]++;
    }
    off += n_recv;
  }
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

// ----------------------------------------------------------------------
// cuda_mparticles_find_block_indices_3

void
cuda_mparticles_find_block_indices_3(struct cuda_mparticles *cmprts)
{
  cudaError_t ierr;
  
  unsigned int nr_recv = cmprts->bnd.n_prts_recv;
  unsigned int nr_prts_prev = cmprts->n_prts - nr_recv;

  // for consistency, use same block indices that we counted earlier
  // OPT unneeded?
  ierr = cudaMemcpy(cmprts->d_bidx + nr_prts_prev, cmprts->bnd.h_bnd_idx,
		    nr_recv * sizeof(*cmprts->d_bidx), cudaMemcpyHostToDevice); cudaCheck(ierr);
  // slight abuse of the now unused last part of spine_cnts
  ierr = cudaMemcpy(cmprts->bnd.d_bnd_spine_cnts + 10 * cmprts->n_blocks,
		    cmprts->bnd.h_bnd_cnt, cmprts->n_blocks * sizeof(*cmprts->bnd.d_bnd_spine_cnts),
		    cudaMemcpyHostToDevice); cudaCheck(ierr);
  ierr = cudaMemcpy(cmprts->bnd.d_alt_bidx + nr_prts_prev, cmprts->bnd.h_bnd_off,
		    nr_recv * sizeof(*cmprts->bnd.d_alt_bidx), cudaMemcpyHostToDevice); cudaCheck(ierr);

  free(cmprts->bnd.h_bnd_idx);
  free(cmprts->bnd.h_bnd_off);
}

// ----------------------------------------------------------------------
// cuda_mparticles_count_received

__global__ static void
mprts_count_received(int nr_total_blocks, unsigned int *d_alt_bidx, unsigned int *d_spine_cnts)
{
  int bid = threadIdx.x + THREADS_PER_BLOCK * blockIdx.x;

  if (bid < nr_total_blocks) {
    d_spine_cnts[bid * 10 + CUDA_BND_S_NEW] = d_alt_bidx[bid];
  }
}

void
cuda_mparticles_count_received(struct cuda_mparticles *cmprts)
{
  unsigned int n_blocks = cmprts->n_blocks;
  
  mprts_count_received<<<n_blocks, THREADS_PER_BLOCK>>>
    (n_blocks, cmprts->bnd.d_bnd_spine_cnts + 10 * n_blocks, cmprts->bnd.d_bnd_spine_cnts);
}

void
cuda_mparticles_count_received_gold(struct cuda_mparticles *cmprts)
{
  int n_blocks = cmprts->n_blocks;

  thrust::device_ptr<unsigned int> d_spine_cnts(cmprts->bnd.d_bnd_spine_cnts);

  thrust::host_vector<unsigned int> h_spine_cnts(1 + n_blocks * (10 + 1));

  thrust::copy(d_spine_cnts, d_spine_cnts + 1 + n_blocks * (10 + 1), h_spine_cnts.begin());

  for (int bid = 0; bid < n_blocks; bid++) {
    h_spine_cnts[bid * 10 + CUDA_BND_S_NEW] = h_spine_cnts[10 * n_blocks + bid];
  }

  thrust::copy(h_spine_cnts.begin(), h_spine_cnts.end(), d_spine_cnts);
}

void
cuda_mparticles_count_received_v1(struct cuda_mparticles *cmprts)
{
  int n_blocks = cmprts->n_blocks;

  thrust::device_ptr<unsigned int> d_bidx(cmprts->d_bidx);
  thrust::device_ptr<unsigned int> d_spine_cnts(cmprts->bnd.d_bnd_spine_cnts);

  thrust::host_vector<unsigned int> h_bidx(cmprts->n_prts);
  thrust::host_vector<unsigned int> h_spine_cnts(1 + n_blocks * (10 + 1));

  thrust::copy(d_bidx, d_bidx + cmprts->n_prts, h_bidx.begin());
  thrust::copy(d_spine_cnts, d_spine_cnts + 1 + n_blocks * (10 + 1), h_spine_cnts.begin());
  for (int n = cmprts->n_prts - cmprts->bnd.n_prts_recv; n < cmprts->n_prts; n++) {
    assert(h_bidx[n] < n_blocks);
    h_spine_cnts[h_bidx[n] * 10 + CUDA_BND_S_NEW]++;
  }
  thrust::copy(h_spine_cnts.begin(), h_spine_cnts.end(), d_spine_cnts);
}

// ----------------------------------------------------------------------
// cuda_mparticles_scan_scatter_received

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
cuda_mparticles_scan_scatter_received(struct cuda_mparticles *cmprts)
{
  int nr_recv = cmprts->bnd.n_prts_recv;

  if (nr_recv == 0) {
    return;
  }
  
  int nr_prts_prev = cmprts->n_prts - nr_recv;

  int dimGrid = (nr_recv + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;

  mprts_scan_scatter_received<<<dimGrid, THREADS_PER_BLOCK>>>
    (nr_recv, nr_prts_prev, cmprts->bnd.d_bnd_spine_sums, cmprts->bnd.d_alt_bidx,
     cmprts->d_bidx, cmprts->d_id);
  cuda_sync_if_enabled();
}

void
cuda_mparticles_scan_scatter_received_gold(struct cuda_mparticles *cmprts)
{
  unsigned int n_blocks = cmprts->n_blocks;

  thrust::device_ptr<unsigned int> d_bidx(cmprts->d_bidx);
  thrust::device_ptr<unsigned int> d_alt_bidx(cmprts->bnd.d_alt_bidx);
  thrust::device_ptr<unsigned int> d_id(cmprts->d_id);
  thrust::device_ptr<unsigned int> d_spine_sums(cmprts->bnd.d_bnd_spine_sums);

  thrust::host_vector<unsigned int> h_bidx(cmprts->n_prts);
  thrust::host_vector<unsigned int> h_alt_bidx(cmprts->n_prts);
  thrust::host_vector<unsigned int> h_id(cmprts->n_prts);
  thrust::host_vector<unsigned int> h_spine_sums(1 + n_blocks * (10 + 1));

  thrust::copy(d_spine_sums, d_spine_sums + n_blocks * 11, h_spine_sums.begin());
  thrust::copy(d_bidx, d_bidx + cmprts->n_prts, h_bidx.begin());
  thrust::copy(d_alt_bidx, d_alt_bidx + cmprts->n_prts, h_alt_bidx.begin());
  for (int n = cmprts->n_prts - cmprts->bnd.n_prts_recv; n < cmprts->n_prts; n++) {
    int nn = h_spine_sums[h_bidx[n] * 10 + CUDA_BND_S_NEW] + h_alt_bidx[n];
    h_id[nn] = n;
  }
  thrust::copy(h_id.begin(), h_id.end(), d_id);
}

// ----------------------------------------------------------------------
// cuda_mparticles_sort_pairs_device

void
cuda_mparticles_sort_pairs_device(struct cuda_mparticles *cmprts)
{
  static int pr_A, pr_B, pr_C, pr_D;
  if (!pr_B) {
    pr_A = prof_register("xchg_cnt_recvd", 1., 0, 0);
    pr_B = prof_register("xchg_top_scan", 1., 0, 0);
    pr_C = prof_register("xchg_ss_recvd", 1., 0, 0);
    pr_D = prof_register("xchg_bottom_scan", 1., 0, 0);
  }

  unsigned int n_blocks = cmprts->n_blocks;

  prof_start(pr_A);
  cuda_mparticles_count_received(cmprts);
  prof_stop(pr_A);

  prof_start(pr_B);
  // FIXME why isn't 10 + 0 enough?
  thrust::device_ptr<unsigned int> d_spine_cnts(cmprts->bnd.d_bnd_spine_cnts);
  thrust::device_ptr<unsigned int> d_spine_sums(cmprts->bnd.d_bnd_spine_sums);
  thrust::exclusive_scan(d_spine_cnts, d_spine_cnts + 1 + n_blocks * (10 + 1), d_spine_sums);
  prof_stop(pr_B);

  prof_start(pr_C);
  cuda_mparticles_scan_scatter_received(cmprts);
  prof_stop(pr_C);

  prof_start(pr_D);
  int *b_mx = cmprts->b_mx;
  if (b_mx[0] == 1 && b_mx[1] == 8 && b_mx[2] == 8) {
    ScanScatterDigits3x<K, V, 0, RADIX_BITS, 0,
			NopFunctor<K>,
			NopFunctor<K>,
			8, 8> 
      <<<n_blocks, B40C_RADIXSORT_THREADS>>>
      (cmprts->bnd.d_bnd_spine_sums, cmprts->d_bidx,
       cmprts->d_id, cmprts->d_off, n_blocks);
  } else if (b_mx[0] == 1 && b_mx[1] == 16 && b_mx[2] == 16) {
    ScanScatterDigits3x<K, V, 0, RADIX_BITS, 0,
			NopFunctor<K>,
			NopFunctor<K>,
			16, 16> 
      <<<n_blocks, B40C_RADIXSORT_THREADS>>>
      (cmprts->bnd.d_bnd_spine_sums, cmprts->d_bidx,
       cmprts->d_id, cmprts->d_off, n_blocks);
  } else if (b_mx[0] == 1 && b_mx[1] == 32 && b_mx[2] == 32) {
    ScanScatterDigits3x<K, V, 0, RADIX_BITS, 0,
			NopFunctor<K>,
			NopFunctor<K>,
			32, 32> 
      <<<n_blocks, B40C_RADIXSORT_THREADS>>>
      (cmprts->bnd.d_bnd_spine_sums, cmprts->d_bidx,
       cmprts->d_id, cmprts->d_off, n_blocks);
  } else if (b_mx[0] == 1 && b_mx[1] == 64 && b_mx[2] == 64) {
    ScanScatterDigits3x<K, V, 0, RADIX_BITS, 0,
			NopFunctor<K>,
			NopFunctor<K>,
			64, 64> 
      <<<n_blocks, B40C_RADIXSORT_THREADS>>>
      (cmprts->bnd.d_bnd_spine_sums, cmprts->d_bidx,
       cmprts->d_id, cmprts->d_off, n_blocks);
  } else if (b_mx[0] == 1 && b_mx[1] == 128 && b_mx[2] == 128) {
    ScanScatterDigits3x<K, V, 0, RADIX_BITS, 0,
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
  prof_stop(pr_D);

  // d_ids now contains the indices to reorder by
}

void
cuda_mprts_sort_pairs_gold(struct cuda_mparticles *cmprts)
{
  unsigned int n_blocks_per_patch = cmprts->n_blocks_per_patch;
  unsigned int n_blocks = cmprts->n_blocks;
  int *b_mx = cmprts->b_mx;

  thrust::device_ptr<unsigned int> d_bidx(cmprts->d_bidx);
  thrust::device_ptr<unsigned int> d_id(cmprts->d_id);
  thrust::device_ptr<unsigned int> d_off(cmprts->d_off);
  thrust::device_ptr<unsigned int> d_spine_cnts(cmprts->bnd.d_bnd_spine_cnts);
  thrust::device_ptr<unsigned int> d_spine_sums(cmprts->bnd.d_bnd_spine_sums);

  thrust::host_vector<unsigned int> h_bidx(d_bidx, d_bidx + cmprts->n_prts);
  thrust::host_vector<unsigned int> h_id(cmprts->n_prts);
  thrust::host_vector<unsigned int> h_off(d_off, d_off + n_blocks + 1);
  thrust::host_vector<unsigned int> h_spine_cnts(d_spine_cnts, d_spine_cnts + 1 + n_blocks * (10 + 1));

  thrust::host_vector<unsigned int> h_spine_sums(1 + n_blocks * (10 + 1));

  for (int n = cmprts->n_prts - cmprts->bnd.n_prts_recv; n < cmprts->n_prts; n++) {
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
  for (int n = cmprts->n_prts - cmprts->bnd.n_prts_recv; n < cmprts->n_prts; n++) {
      int nn = h_spine_sums[h_bidx[n] * 10 + CUDA_BND_S_NEW]++;
      h_id[nn] = n;
  }

  thrust::copy(h_id.begin(), h_id.end(), d_id);
  // d_ids now contains the indices to reorder by
}

// ----------------------------------------------------------------------
// cuda_mparticles_sort

void
cuda_mparticles_sort(struct cuda_mparticles *cmprts, int *n_prts_by_patch)
{
  cuda_mparticles_sort_pairs_device(cmprts);

  for (int p = 0; p < cmprts->n_patches; p++) {
    n_prts_by_patch[p] += cmprts->bnd.bpatch[p].n_recv - cmprts->bnd.bpatch[p].n_send;
  }
  cmprts->n_prts -= cmprts->bnd.n_prts_send;
}

// ----------------------------------------------------------------------
// cuda_mparticles_update_offsets

__global__ static void
mprts_update_offsets(int nr_total_blocks, unsigned int *d_off, unsigned int *d_spine_sums)
{
  int bid = threadIdx.x + THREADS_PER_BLOCK * blockIdx.x;
  
  if (bid <= nr_total_blocks) {
    d_off[bid] = d_spine_sums[bid * CUDA_BND_STRIDE + 0];
  }
}

void
cuda_mparticles_update_offsets(struct cuda_mparticles *cmprts)
{
  unsigned int n_blocks = cmprts->n_blocks;
  int dimGrid = (n_blocks + 1 + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;

  mprts_update_offsets<<<dimGrid, THREADS_PER_BLOCK>>>
    (n_blocks, cmprts->d_off, cmprts->bnd.d_bnd_spine_sums);
  cuda_sync_if_enabled();
}

void
cuda_mparticles_update_offsets_gold(struct cuda_mparticles *cmprts)
{
  unsigned int n_blocks = cmprts->n_blocks;

  thrust::device_ptr<unsigned int> d_spine_sums(cmprts->bnd.d_bnd_spine_sums);
  thrust::device_ptr<unsigned int> d_off(cmprts->d_off);

  thrust::host_vector<unsigned int> h_spine_sums(d_spine_sums, d_spine_sums + 1 + n_blocks * (10 + 1));
  thrust::host_vector<unsigned int> h_off(n_blocks + 1);

  for (int bid = 0; bid <= n_blocks; bid++) {
    h_off[bid] = h_spine_sums[bid * 10];
  }

  thrust::copy(h_off.begin(), h_off.end(), d_off);
}

// ----------------------------------------------------------------------
// cuda_mparticles_bnd_prep

void
cuda_mparticles_bnd_prep(struct cuda_mparticles *cmprts)
{
  static int pr_A, pr_B, pr_D, pr_E, pr_B0, pr_B1;
  if (!pr_A) {
    pr_A = prof_register("xchg_bidx", 1., 0, 0);
    pr_B0= prof_register("xchg_reduce", 1., 0, 0);
    pr_B1= prof_register("xchg_n_send", 1., 0, 0);
    pr_B = prof_register("xchg_scan_send", 1., 0, 0);
    pr_D = prof_register("xchg_from_dev", 1., 0, 0);
    pr_E = prof_register("xchg_cvt_from", 1., 0, 0);
  }

  //prof_start(pr_A);
  //cuda_mprts_find_block_keys(mprts);
  //prof_stop(pr_A);
  
  prof_start(pr_B0);
  cuda_mparticles_spine_reduce(cmprts);
  prof_stop(pr_B0);

  prof_start(pr_B1);
  cuda_mparticles_find_n_send(cmprts);
  prof_stop(pr_B1);

  prof_start(pr_B);
  cuda_mparticles_scan_send_buf_total(cmprts);
  prof_stop(pr_B);

  prof_start(pr_D);
  cuda_mparticles_copy_from_dev(cmprts);
  prof_stop(pr_D);

  prof_start(pr_E);
  cuda_mparticles_convert_from_cuda(cmprts);
  prof_stop(pr_E);
}

// ----------------------------------------------------------------------
// cuda_mparticles_bnd_post

void
cuda_mparticles_bnd_post(struct cuda_mparticles *cmprts)
{
  static int pr_A, pr_B, pr_C, pr_D, pr_E, pr_D1;
  if (!pr_A) {
    pr_A = prof_register("xchg_cvt_to", 1., 0, 0);
    pr_B = prof_register("xchg_to_dev", 1., 0, 0);
    pr_C = prof_register("xchg_bidx", 1., 0, 0);
    pr_D = prof_register("xchg_sort", 1., 0, 0);
    pr_D1= prof_register("xchg_upd_off", 1., 0, 0);
    pr_E = prof_register("xchg_reorder", 1., 0, 0);
  }

  prof_start(pr_A);
  cuda_mparticles_convert_to_cuda(cmprts);
  prof_stop(pr_A);

  prof_start(pr_B);
  cuda_mparticles_copy_to_dev(cmprts);
  prof_stop(pr_B);

  prof_start(pr_C);
  cuda_mparticles_find_block_indices_3(cmprts);
  prof_stop(pr_C);
  
  prof_start(pr_D);
  unsigned int n_prts_by_patch[cmprts->n_patches];
  cuda_mparticles_get_size_all(cmprts, n_prts_by_patch);
  cuda_mparticles_sort(cmprts, (int *) n_prts_by_patch); // FIXME cast
  // FIXME, is this necessary, or doesn't update_offsets() do this, too?
  cuda_mparticles_resize_all(cmprts, n_prts_by_patch);
  prof_stop(pr_D);

  prof_start(pr_D1);
  cuda_mparticles_update_offsets(cmprts);
  prof_stop(pr_D1);
  
  prof_start(pr_E);
#if 0
  cuda_mparticles_reorder(cmprts);
  //  cuda_mprts_check_ordered_total(mprts);
#else
  cmprts->need_reorder = true;
#endif
  prof_stop(pr_E);

  for (int p = 0; p < cmprts->n_patches; p++) {
    psc_particle_cuda_buf_dtor(&cmprts->bnd.bpatch[p].buf);
  }
}
