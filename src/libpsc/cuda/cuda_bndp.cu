
#include "cuda_bndp.h"
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

// ----------------------------------------------------------------------
// ctor

cuda_bndp::cuda_bndp(const Grid_t& grid)
  : cuda_mparticles_indexer(grid)
{
  d_spine_cnts.resize(1 + n_blocks * (CUDA_BND_STRIDE + 1));
  d_spine_sums.resize(1 + n_blocks * (CUDA_BND_STRIDE + 1));

  bpatch.resize(n_patches);
}

// ----------------------------------------------------------------------
// prep

void cuda_bndp::prep(ddcp_t* ddcp, cuda_mparticles* cmprts)
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
  spine_reduce(cmprts);
  prof_stop(pr_B0);

  prof_start(pr_B1);
  n_prts_send = find_n_send(cmprts);
  prof_stop(pr_B1);

  prof_start(pr_B);
  scan_send_buf_total(cmprts, n_prts_send);
  prof_stop(pr_B);

  prof_start(pr_D);
  copy_from_dev_and_convert(cmprts, n_prts_send);
  prof_stop(pr_D);

  if (!ddcp) return; // FIXME testing hack
  for (int p = 0; p < ddcp->nr_patches; p++) {
    ddcp_patch *dpatch = &ddcp->patches[p];
    dpatch->m_buf = &bpatch[p].buf;
    dpatch->m_begin = 0;
  }
}

// ----------------------------------------------------------------------
// post

void cuda_bndp::post(ddcp_t* ddcp, cuda_mparticles* cmprts)
{
  static int pr_A, pr_D, pr_E, pr_D1;
  if (!pr_A) {
    pr_A = prof_register("xchg_to_dev", 1., 0, 0);
    pr_D = prof_register("xchg_sort", 1., 0, 0);
    pr_D1= prof_register("xchg_upd_off", 1., 0, 0);
    pr_E = prof_register("xchg_reorder", 1., 0, 0);
  }

  prof_start(pr_A);
  uint n_prts_recv = convert_and_copy_to_dev(cmprts);
  cmprts->n_prts += n_prts_recv;
  prof_stop(pr_A);

  prof_start(pr_D);
  sort_pairs_device(cmprts, n_prts_recv);
  cmprts->n_prts -= n_prts_send;
  prof_stop(pr_D);

  prof_start(pr_D1);
  update_offsets(cmprts);
  prof_stop(pr_D1);
  
  prof_start(pr_E);
#if 0
  cmprts->reorder(cmprts);
  assert(cmprts->check_ordered());
#else
  cmprts->need_reorder = true;
#endif
  prof_stop(pr_E);

  if (!ddcp) return; // FIXME, testing hack
  for (int p = 0; p < ddcp->nr_patches; p++) {
    ddcp->patches[p].m_buf = NULL;
  }
}

// ----------------------------------------------------------------------
// find_n_send

uint cuda_bndp::find_n_send(cuda_mparticles *cmprts)
{
  thrust::host_vector<uint> h_spine_sums(n_blocks + 1);

  thrust::copy(d_spine_sums.data() + n_blocks * 10,
	       d_spine_sums.data() + n_blocks * 11 + 1,
	       h_spine_sums.begin());

  uint off = 0;
  for (int p = 0; p < n_patches; p++) {
    uint n_send = h_spine_sums[(p + 1) * n_blocks_per_patch];
    bpatch[p].n_send = n_send - off;
    off = n_send;
  }
  return off;
}

// ----------------------------------------------------------------------
// copy_from_dev_and_convert

void cuda_bndp::copy_from_dev_and_convert(cuda_mparticles *cmprts, uint n_prts_send)
{
  uint n_prts = cmprts->n_prts;
  thrust::host_vector<float4> h_bnd_xi4(n_prts_send);
  thrust::host_vector<float4> h_bnd_pxi4(n_prts_send);

  assert(cmprts->d_xi4.begin() + n_prts + n_prts_send == cmprts->d_xi4.end());

  thrust::copy(cmprts->d_xi4.begin()  + n_prts, cmprts->d_xi4.end(), h_bnd_xi4.begin());
  thrust::copy(cmprts->d_pxi4.begin() + n_prts, cmprts->d_pxi4.end(), h_bnd_pxi4.begin());

  uint off = 0;
  for (int p = 0; p < n_patches; p++) {
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
      prt->qni_wni_ = h_bnd_pxi4[n + off].w;
    }
    off += n_send;
  }

  cmprts->resize(n_prts);
}

// ----------------------------------------------------------------------
// convert_and_copy_to_dev

uint cuda_bndp::convert_and_copy_to_dev(cuda_mparticles *cmprts)
{
  uint n_recv = 0;
  for (int p = 0; p < n_patches; p++) {
    n_recv += bpatch[p].buf.size();
  }

  thrust::host_vector<float4> h_bnd_xi4(n_recv);
  thrust::host_vector<float4> h_bnd_pxi4(n_recv);
  thrust::host_vector<uint> h_bnd_idx(n_recv);
  thrust::host_vector<uint> h_bnd_off(n_recv);

  thrust::host_vector<uint> h_bnd_cnt(n_blocks, 0);
  
  uint off = 0;
  for (int p = 0; p < n_patches; p++) {
    int n_recv = bpatch[p].buf.size();
    bpatch[p].n_recv = n_recv;
    
    for (int n = 0; n < n_recv; n++) {
      const particle_cuda_t& prt = bpatch[p].buf[n];

      h_bnd_xi4[n + off].x  = prt.xi;
      h_bnd_xi4[n + off].y  = prt.yi;
      h_bnd_xi4[n + off].z  = prt.zi;
      h_bnd_xi4[n + off].w  = cuda_int_as_float(prt.kind_);
      h_bnd_pxi4[n + off].x = prt.pxi;
      h_bnd_pxi4[n + off].y = prt.pyi;
      h_bnd_pxi4[n + off].z = prt.pzi;
      h_bnd_pxi4[n + off].w = prt.qni_wni_;

      pi_.checkInPatchMod(&h_bnd_xi4[n + off].x);
      uint b = blockIndex(h_bnd_xi4[n + off], p);
      assert(b < n_blocks);
      h_bnd_idx[n + off] = b;
      h_bnd_off[n + off] = h_bnd_cnt[b]++;
    }
    off += n_recv;
  }

  cmprts->resize(cmprts->n_prts + n_recv);

  thrust::copy(h_bnd_xi4.begin(), h_bnd_xi4.end(), cmprts->d_xi4.begin() + cmprts->n_prts);
  thrust::copy(h_bnd_pxi4.begin(), h_bnd_pxi4.end(), cmprts->d_pxi4.begin() + cmprts->n_prts);

  // for consistency, use same block indices that we counted earlier
  // OPT unneeded?
  thrust::copy(h_bnd_idx.begin(), h_bnd_idx.end(), cmprts->d_bidx.begin() + cmprts->n_prts);
  // slight abuse of the now unused last part of spine_cnts
  thrust::copy(h_bnd_cnt.begin(), h_bnd_cnt.end(), d_spine_cnts.begin() + 10 * n_blocks);

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

void cuda_bndp::update_offsets(cuda_mparticles *cmprts)
{
  int dimGrid = (n_blocks + 1 + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;

  mprts_update_offsets<<<dimGrid, THREADS_PER_BLOCK>>>
    (n_blocks, cmprts->d_off.data().get(), d_spine_sums.data().get());
  cuda_sync_if_enabled();
}

void cuda_bndp::update_offsets_gold(cuda_mparticles *cmprts)
{
  thrust::host_vector<uint> h_spine_sums(d_spine_sums.data(), d_spine_sums.data() + 1 + n_blocks * (10 + 1));
  thrust::host_vector<uint> h_off(n_blocks + 1);

  for (int bid = 0; bid <= n_blocks; bid++) {
    h_off[bid] = h_spine_sums[bid * 10];
  }

  thrust::copy(h_off.begin(), h_off.end(), cmprts->d_off.begin());
}

