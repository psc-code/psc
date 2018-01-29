
#include "cuda_particles_bnd.h"
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
// reserve_all

void cuda_mparticles_bnd::reserve_all(cuda_mparticles *cmprts)
{
  int n_alloced = cmprts->n_alloced;
  d_sums.resize(n_alloced);
}

// ----------------------------------------------------------------------
// find_n_send

uint cuda_particles_bnd::find_n_send(cuda_mparticles *cmprts)
{
  uint n_blocks = cmprts->n_blocks;

  thrust::host_vector<uint> h_spine_sums(n_blocks + 1);

  thrust::copy(d_spine_sums.data() + n_blocks * 10,
	       d_spine_sums.data() + n_blocks * 11 + 1,
	       h_spine_sums.begin());

  uint off = 0;
  for (int p = 0; p < cmprts->n_patches; p++) {
    uint n_send = h_spine_sums[(p + 1) * cmprts->n_blocks_per_patch];
    bpatch[p].n_send = n_send - off;
    off = n_send;
  }
  return off;
}

// ----------------------------------------------------------------------
// copy_from_dev_and_convert

void cuda_particles_bnd::copy_from_dev_and_convert(cuda_mparticles *cmprts, uint n_prts_send)
{
  uint n_prts = cmprts->n_prts;
  thrust::host_vector<float4> h_bnd_xi4(n_prts_send);
  thrust::host_vector<float4> h_bnd_pxi4(n_prts_send);

  assert(n_prts + n_prts_send <= cmprts->n_alloced);

  thrust::copy(cmprts->d_xi4.data()  + n_prts, cmprts->d_xi4.data()  + n_prts + n_prts_send, h_bnd_xi4.begin());
  thrust::copy(cmprts->d_pxi4.data() + n_prts, cmprts->d_pxi4.data() + n_prts + n_prts_send, h_bnd_pxi4.begin());

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

uint cuda_particles_bnd::convert_and_copy_to_dev(cuda_mparticles *cmprts)
{
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
	b_pos[d] = fint(xi[d] * cmprts->indexer.b_dxi_[d]);
	if (b_pos[d] < 0 || b_pos[d] >= cmprts->indexer.b_mx_[d]) {
	  printf("!!! xi %g %g %g\n", xi[0], xi[1], xi[2]);
	  printf("!!! d %d xi4[n] %g biy %d // %d\n",
		 d, xi[d], b_pos[d], cmprts->indexer.b_mx_[d]);
	  if (b_pos[d] < 0) {
	    xi[d] = 0.f;
	  } else {
	    xi[d] *= (1. - 1e-6);
	  }
	}
	b_pos[d] = fint(xi[d] * cmprts->indexer.b_dxi_[d]);
	assert(b_pos[d] >= 0 && b_pos[d] < cmprts->indexer.b_mx_[d]);
      }
      uint b = (b_pos[2] * cmprts->indexer.b_mx_[1] + b_pos[1]) * cmprts->indexer.b_mx_[0] + b_pos[0];
      assert(b < cmprts->n_blocks_per_patch);
      b += p * cmprts->n_blocks_per_patch;
      h_bnd_idx[n + off] = b;
      h_bnd_off[n + off] = h_bnd_cnt[b]++;
    }
    off += n_recv;
  }

  assert(cmprts->n_prts + n_recv <= cmprts->n_alloced);

  thrust::copy(h_bnd_xi4.begin(), h_bnd_xi4.end(), cmprts->d_xi4.data() + cmprts->n_prts);
  thrust::copy(h_bnd_pxi4.begin(), h_bnd_pxi4.end(), cmprts->d_pxi4.data() + cmprts->n_prts);

  // for consistency, use same block indices that we counted earlier
  // OPT unneeded?
  thrust::copy(h_bnd_idx.begin(), h_bnd_idx.end(), cmprts->d_bidx.data() + cmprts->n_prts);
  // slight abuse of the now unused last part of spine_cnts
  thrust::copy(h_bnd_cnt.begin(), h_bnd_cnt.end(),
	       d_spine_cnts.data() + 10 * cmprts->n_blocks);

  d_bnd_off.resize(n_recv);
  thrust::copy(h_bnd_off.begin(), h_bnd_off.end(), d_bnd_off.begin());

  return n_recv;
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

void cuda_particles_bnd::update_offsets(cuda_mparticles *cmprts)
{
  uint n_blocks = cmprts->n_blocks;
  int dimGrid = (n_blocks + 1 + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;

  mprts_update_offsets<<<dimGrid, THREADS_PER_BLOCK>>>
    (n_blocks, cmprts->d_off.data().get(), d_spine_sums.data().get());
  cuda_sync_if_enabled();
}

void cuda_particles_bnd::update_offsets_gold(cuda_mparticles *cmprts)
{
  uint n_blocks = cmprts->n_blocks;

  thrust::host_vector<uint> h_spine_sums(d_spine_sums.data(), d_spine_sums.data() + 1 + n_blocks * (10 + 1));
  thrust::host_vector<uint> h_off(n_blocks + 1);

  for (int bid = 0; bid <= n_blocks; bid++) {
    h_off[bid] = h_spine_sums[bid * 10];
  }

  thrust::copy(h_off.begin(), h_off.end(), cmprts->d_off.begin());
}

