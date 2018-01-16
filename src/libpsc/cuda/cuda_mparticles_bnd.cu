
#include "cuda_mparticles.h"
#include "cuda_bits.h"

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/scan.h>

#include "psc_bits.h"

#include <mrc_profile.h>

#define THREADS_PER_BLOCK 256

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

// ======================================================================
// cuda_mparticles_bnd

// ----------------------------------------------------------------------
// setup

void cuda_mparticles_bnd::setup(struct cuda_mparticles *cmprts)
{
  cudaError_t ierr;

  ierr = cudaMalloc((void **) &d_bnd_spine_cnts,
		    (1 + cmprts->n_blocks * (CUDA_BND_STRIDE + 1)) * sizeof(uint)); cudaCheck(ierr);
  ierr = cudaMalloc((void **) &d_bnd_spine_sums,
		    (1 + cmprts->n_blocks * (CUDA_BND_STRIDE + 1)) * sizeof(uint)); cudaCheck(ierr);

  bpatch = new cuda_bnd[cmprts->n_patches];
}  

// ----------------------------------------------------------------------
// dtor

cuda_mparticles_bnd::~cuda_mparticles_bnd()
{
  cudaError_t ierr;

  ierr = cudaFree(d_bnd_spine_cnts); cudaCheck(ierr);
  ierr = cudaFree(d_bnd_spine_sums); cudaCheck(ierr);

  delete[] bpatch;
}

// ----------------------------------------------------------------------
// free_particle_mem

void cuda_mparticles_bnd::free_particle_mem()
{
  cudaError_t ierr;

  ierr = cudaFree(d_sums); cudaCheck(ierr);
}

// ----------------------------------------------------------------------
// reserve_all

void cuda_mparticles_bnd::reserve_all(cuda_mparticles *cmprts)
{
  cudaError_t ierr;

  int n_alloced = cmprts->n_alloced;
  d_alt_bidx.resize(n_alloced);
  ierr = cudaMalloc((void **) &d_sums, n_alloced * sizeof(uint)); cudaCheck(ierr);
}

// ----------------------------------------------------------------------
// find_n_send

void cuda_mparticles_bnd::find_n_send(cuda_mparticles *cmprts)
{
  uint n_blocks = cmprts->n_blocks;

  thrust::device_ptr<uint> d_spine_sums(d_bnd_spine_sums);
  thrust::host_vector<uint> h_spine_sums(n_blocks + 1);

  thrust::copy(d_spine_sums + n_blocks * 10,
	       d_spine_sums + n_blocks * 11 + 1,
	       h_spine_sums.begin());

  uint off = 0;
  for (int p = 0; p < cmprts->n_patches; p++) {
    uint n_send = h_spine_sums[(p + 1) * cmprts->n_blocks_per_patch];
    bpatch[p].n_send = n_send - off;
    off = n_send;
  }
  n_prts_send = off;
}

// ----------------------------------------------------------------------
// copy_from_dev_and_convert

void cuda_mparticles_bnd::copy_from_dev_and_convert(cuda_mparticles *cmprts)
{
  thrust::device_ptr<float4> d_xi4(cmprts->d_xi4);
  thrust::device_ptr<float4> d_pxi4(cmprts->d_pxi4);

  thrust::host_vector<float4> h_bnd_xi4(n_prts_send);
  thrust::host_vector<float4> h_bnd_pxi4(n_prts_send);

  assert(cmprts->n_prts + n_prts_send <= cmprts->n_alloced);

  thrust::copy(d_xi4  + cmprts->n_prts, d_xi4  + cmprts->n_prts + n_prts_send, h_bnd_xi4.begin());
  thrust::copy(d_pxi4 + cmprts->n_prts, d_pxi4 + cmprts->n_prts + n_prts_send, h_bnd_pxi4.begin());

  uint off = 0;
  for (int p = 0; p < cmprts->n_patches; p++) {
    psc_particle_cuda_buf_t& buf = bpatch[p].buf;
    uint n_send = bpatch[p].n_send;
    buf.reserve(n_send);
    buf.resize(n_send);

    for (int n = 0; n < n_send; n++) {
      particle_cuda_t *prt = &buf[n];
      prt->xi      = h_bnd_xi4[n + off].x;
      prt->yi      = h_bnd_xi4[n + off].y;
      prt->zi      = h_bnd_xi4[n + off].z;
      prt->kind_   = cuda_float_as_int(h_bnd_xi4[n + off].w);
      prt->pxi     = h_bnd_pxi4[n + off].x;
      prt->pyi     = h_bnd_pxi4[n + off].y;
      prt->pzi     = h_bnd_pxi4[n + off].z;
      prt->qni_wni = h_bnd_pxi4[n + off].w;
    }
    off += n_send;
  }
}

// ----------------------------------------------------------------------
// convert_and_copy_to_dev

void cuda_mparticles_bnd::convert_and_copy_to_dev(cuda_mparticles *cmprts)
{
  thrust::device_ptr<float4> d_xi4(cmprts->d_xi4);
  thrust::device_ptr<float4> d_pxi4(cmprts->d_pxi4);
  thrust::device_ptr<uint> d_bidx(cmprts->d_bidx);
  thrust::device_ptr<uint> d_bnd_spine_cnts(this->d_bnd_spine_cnts);
  
  uint n_recv = 0;
  for (int p = 0; p < cmprts->n_patches; p++) {
    n_recv += bpatch[p].buf.size();
  }

  thrust::host_vector<float4> h_bnd_xi4(n_recv);
  thrust::host_vector<float4> h_bnd_pxi4(n_recv);
  thrust::host_vector<uint> h_bnd_idx(n_recv);
  thrust::host_vector<uint> h_bnd_off(n_recv);

  thrust::host_vector<uint> h_bnd_cnt(cmprts->n_blocks, 0);
  
  uint off = 0;
  for (int p = 0; p < cmprts->n_patches; p++) {
    int n_recv = bpatch[p].buf.size();
    bpatch[p].n_recv = n_recv;
    
    for (int n = 0; n < n_recv; n++) {
      particle_cuda_t *prt = &bpatch[p].buf[n];
      h_bnd_xi4[n + off].x  = prt->xi;
      h_bnd_xi4[n + off].y  = prt->yi;
      h_bnd_xi4[n + off].z  = prt->zi;
      h_bnd_xi4[n + off].w  = cuda_int_as_float(prt->kind_);
      h_bnd_pxi4[n + off].x = prt->pxi;
      h_bnd_pxi4[n + off].y = prt->pyi;
      h_bnd_pxi4[n + off].z = prt->pzi;
      h_bnd_pxi4[n + off].w = prt->qni_wni;

      int b_pos[3];
      for (int d = 0; d < 3; d++) {
	float *xi = &h_bnd_xi4[n + off].x;
	b_pos[d] = fint(xi[d] * cmprts->b_dxi[d]);
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
	b_pos[d] = fint(xi[d] * cmprts->b_dxi[d]);
	assert(b_pos[d] >= 0 && b_pos[d] < cmprts->b_mx[d]);
      }
      uint b = (b_pos[2] * cmprts->b_mx[1] + b_pos[1]) * cmprts->b_mx[0] + b_pos[0];
      assert(b < cmprts->n_blocks_per_patch);
      b += p * cmprts->n_blocks_per_patch;
      h_bnd_idx[n + off] = b;
      h_bnd_off[n + off] = h_bnd_cnt[b]++;
    }
    off += n_recv;
  }

  assert(cmprts->n_prts + n_recv <= cmprts->n_alloced);

  thrust::copy(h_bnd_xi4.begin(), h_bnd_xi4.end(), d_xi4 + cmprts->n_prts);
  thrust::copy(h_bnd_pxi4.begin(), h_bnd_pxi4.end(), d_pxi4 + cmprts->n_prts);

  // for consistency, use same block indices that we counted earlier
  // OPT unneeded?
  thrust::copy(h_bnd_idx.begin(), h_bnd_idx.end(), d_bidx + cmprts->n_prts);
  // slight abuse of the now unused last part of spine_cnts
  thrust::copy(h_bnd_cnt.begin(), h_bnd_cnt.end(),
	       d_bnd_spine_cnts + 10 * cmprts->n_blocks);
  thrust::copy(h_bnd_off.begin(), h_bnd_off.end(), d_alt_bidx.begin() + cmprts->n_prts);

  cmprts->n_prts += n_recv;
  n_prts_recv = n_recv;
}

// ----------------------------------------------------------------------
// cuda_mparticles_bnd::sort

void cuda_mparticles_bnd::sort(cuda_mparticles *cmprts, int *n_prts_by_patch)
{
  sort_pairs_device(cmprts);

  for (int p = 0; p < cmprts->n_patches; p++) {
    n_prts_by_patch[p] += bpatch[p].n_recv - bpatch[p].n_send;
  }
  cmprts->n_prts -= n_prts_send;
}

// ----------------------------------------------------------------------
// update_offsets

__global__ static void
mprts_update_offsets(int nr_total_blocks, uint *d_off, uint *d_spine_sums)
{
  int bid = threadIdx.x + THREADS_PER_BLOCK * blockIdx.x;
  
  if (bid <= nr_total_blocks) {
    d_off[bid] = d_spine_sums[bid * CUDA_BND_STRIDE + 0];
  }
}

void cuda_mparticles_bnd::update_offsets(cuda_mparticles *cmprts)
{
  uint n_blocks = cmprts->n_blocks;
  int dimGrid = (n_blocks + 1 + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;

  mprts_update_offsets<<<dimGrid, THREADS_PER_BLOCK>>>
    (n_blocks, cmprts->d_off, d_bnd_spine_sums);
  cuda_sync_if_enabled();
}

void cuda_mparticles_bnd::update_offsets_gold(cuda_mparticles *cmprts)
{
  uint n_blocks = cmprts->n_blocks;

  thrust::device_ptr<uint> d_spine_sums(d_bnd_spine_sums);
  thrust::device_ptr<uint> d_off(cmprts->d_off);

  thrust::host_vector<uint> h_spine_sums(d_spine_sums, d_spine_sums + 1 + n_blocks * (10 + 1));
  thrust::host_vector<uint> h_off(n_blocks + 1);

  for (int bid = 0; bid <= n_blocks; bid++) {
    h_off[bid] = h_spine_sums[bid * 10];
  }

  thrust::copy(h_off.begin(), h_off.end(), d_off);
}

// ----------------------------------------------------------------------
// bnd_prep

void cuda_mparticles::bnd_prep()
{
  static int pr_A, pr_B, pr_D, pr_B0, pr_B1;
  if (!pr_A) {
    pr_A = prof_register("xchg_bidx", 1., 0, 0);
    pr_B0= prof_register("xchg_reduce", 1., 0, 0);
    pr_B1= prof_register("xchg_n_send", 1., 0, 0);
    pr_B = prof_register("xchg_scan_send", 1., 0, 0);
    pr_D = prof_register("xchg_from_dev", 1., 0, 0);
  }

  //prof_start(pr_A);
  //cuda_mprts_find_block_keys(mprts);
  //prof_stop(pr_A);
  
  prof_start(pr_B0);
  spine_reduce(this);
  prof_stop(pr_B0);

  prof_start(pr_B1);
  find_n_send(this);
  prof_stop(pr_B1);

  prof_start(pr_B);
  scan_send_buf_total(this);
  prof_stop(pr_B);

  prof_start(pr_D);
  copy_from_dev_and_convert(this);
  prof_stop(pr_D);
}

// ----------------------------------------------------------------------
// cuda_mparticles_bnd_post

void cuda_mparticles::bnd_post()
{
  static int pr_A, pr_D, pr_E, pr_D1;
  if (!pr_A) {
    pr_A = prof_register("xchg_to_dev", 1., 0, 0);
    pr_D = prof_register("xchg_sort", 1., 0, 0);
    pr_D1= prof_register("xchg_upd_off", 1., 0, 0);
    pr_E = prof_register("xchg_reorder", 1., 0, 0);
  }

  prof_start(pr_A);
  convert_and_copy_to_dev(this);
  prof_stop(pr_A);

  prof_start(pr_D);
  uint n_prts_by_patch[n_patches];
  get_size_all(n_prts_by_patch);
  sort(this, (int *) n_prts_by_patch); // FIXME cast
  // FIXME, is this necessary, or doesn't update_offsets() do this, too?
  resize_all(n_prts_by_patch);
  prof_stop(pr_D);

  prof_start(pr_D1);
  update_offsets(this);
  prof_stop(pr_D1);
  
  prof_start(pr_E);
#if 0
  cuda_mparticles_reorder(this);
  //  check_ordered_total();
#else
  need_reorder = true;
#endif
  prof_stop(pr_E);
}

// ----------------------------------------------------------------------
// cuda_mparticles_bnd_get_buffer

psc_particle_cuda_buf_t *cuda_mparticles::bnd_get_buffer(int p)
{
  return &bpatch[p].buf;
}

