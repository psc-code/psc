
#include "cuda_bits.h"
#include "cuda_bndp.h"
#include "cuda_mparticles.hxx"

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/scan.h>
#include <thrust/sequence.h>

#include "psc_bits.h"

#include <mrc_profile.h>

#define THREADS_PER_BLOCK 256

// clang-format off
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
// clang-format on

#include <cstdio>
#include <cassert>

#ifdef CUDA_BNDP_DIM_YZ_SPECIAL

// ----------------------------------------------------------------------
// ctor

template <typename CudaMparticles>
cuda_bndp<CudaMparticles, dim_yz>::cuda_bndp(const Grid_t& grid)
  : cuda_mparticles_indexer<BS>(grid)
{
  d_spine_cnts.resize(1 + n_blocks * (CUDA_BND_STRIDE + 1));
  d_spine_sums.resize(1 + n_blocks * (CUDA_BND_STRIDE + 1));

  bufs.resize(n_patches());
  n_sends.resize(n_patches());
  n_recvs.resize(n_patches());
}

// ----------------------------------------------------------------------
// prep

template <typename CudaMparticles>
auto cuda_bndp<CudaMparticles, dim_yz>::prep(CudaMparticles* cmprts)
  -> BndBuffers&
{
  static int pr_A, pr_B, pr_D, pr_B0, pr_B1;
  if (!pr_A) {
    pr_A = prof_register("xchg_bidx", 1., 0, 0);
    pr_B0 = prof_register("xchg_reduce", 1., 0, 0);
    pr_B1 = prof_register("xchg_n_send", 1., 0, 0);
    pr_B = prof_register("xchg_scan_send", 1., 0, 0);
    pr_D = prof_register("xchg_from_dev", 1., 0, 0);
  }

  // prof_start(pr_A);
  // cuda_mprts_find_block_keys(mprts);
  // prof_stop(pr_A);

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

  return bufs;
}

// ----------------------------------------------------------------------
// post

template <typename CudaMparticles>
void cuda_bndp<CudaMparticles, dim_yz>::post(CudaMparticles* cmprts)
{
  static int pr_A, pr_D, pr_E, pr_D1;
  if (!pr_A) {
    pr_A = prof_register("xchg_to_dev", 1., 0, 0);
    pr_D = prof_register("xchg_sort", 1., 0, 0);
    pr_D1 = prof_register("xchg_upd_off", 1., 0, 0);
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
}

// ----------------------------------------------------------------------
// find_n_send

template <typename CudaMparticles>
uint cuda_bndp<CudaMparticles, dim_yz>::find_n_send(CudaMparticles* cmprts)
{
  thrust::host_vector<uint> h_spine_sums(n_blocks + 1);

  thrust::copy(d_spine_sums.data() + n_blocks * 10,
               d_spine_sums.data() + n_blocks * 11 + 1, h_spine_sums.begin());

  uint off = 0;
  for (int p = 0; p < n_patches(); p++) {
    uint n_send = h_spine_sums[(p + 1) * n_blocks_per_patch];
    n_sends[p] = n_send - off;
    off = n_send;
  }
  return off;
}

// ----------------------------------------------------------------------
// copy_from_dev_and_convert

template <typename CudaMparticles>
void cuda_bndp<CudaMparticles, dim_yz>::copy_from_dev_and_convert(
  CudaMparticles* cmprts, uint n_prts_send)
{
  uint n_prts = cmprts->n_prts;
  HMparticlesCudaStorage h_bnd_storage(n_prts_send);

  assert(cmprts->storage.size() == n_prts + n_prts_send);

  thrust::copy(cmprts->storage.begin() + n_prts, cmprts->storage.end(),
               h_bnd_storage.begin());

  uint off = 0;
  for (int p = 0; p < n_patches(); p++) {
    auto& buf = bufs[p];
    uint n_send = n_sends[p];
    buf.resize(n_send);

    for (int n = 0; n < n_send; n++) {
      buf[n] = h_bnd_storage.load(n + off);
    }
    off += n_send;
  }

  cmprts->resize(n_prts);
}

// ----------------------------------------------------------------------
// convert_and_copy_to_dev

template <typename CudaMparticles>
uint cuda_bndp<CudaMparticles, dim_yz>::convert_and_copy_to_dev(
  CudaMparticles* cmprts)
{
  uint n_recv = 0;
  for (int p = 0; p < n_patches(); p++) {
    n_recv += bufs[p].size();
  }

  HMparticlesCudaStorage h_bnd_storage{n_recv};
  thrust::host_vector<uint> h_bnd_idx(n_recv);
  thrust::host_vector<uint> h_bnd_off(n_recv);

  thrust::host_vector<uint> h_bnd_cnt(n_blocks, 0);

  uint off = 0;
  for (int p = 0; p < n_patches(); p++) {
    int n_recv = bufs[p].size();
    n_recvs[p] = n_recv;

    for (int n = 0; n < n_recv; n++) {
      h_bnd_storage.store(bufs[p][n], n + off);
      checkInPatchMod(h_bnd_storage[n + off].x);
      uint b = blockIndex(h_bnd_storage[n + off], p);
      assert(b < n_blocks);
      h_bnd_idx[n + off] = b;
      h_bnd_off[n + off] = h_bnd_cnt[b]++;
    }
    off += n_recv;
  }

  cmprts->resize(cmprts->n_prts + n_recv);

#ifdef PSC_HAVE_RMM
  thrust::copy(rmm::exec_policy(), h_bnd_storage.begin(), h_bnd_storage.end(),
               cmprts->storage.begin() + cmprts->n_prts);
#else
  thrust::copy(h_bnd_storage.begin(), h_bnd_storage.end(),
               cmprts->storage.begin() + cmprts->n_prts);
#endif

  // for consistency, use same block indices that we counted earlier
  // OPT unneeded?
  thrust::copy(h_bnd_idx.begin(), h_bnd_idx.end(),
               cmprts->by_block_.d_idx.begin() + cmprts->n_prts);
  // slight abuse of the now unused last part of spine_cnts
  thrust::copy(h_bnd_cnt.begin(), h_bnd_cnt.end(),
               d_spine_cnts.begin() + 10 * n_blocks);

  d_bnd_off.resize(n_recv);
  thrust::copy(h_bnd_off.begin(), h_bnd_off.end(), d_bnd_off.begin());

  return n_recv;
}

// ----------------------------------------------------------------------
// update_offsets

__global__ static void mprts_update_offsets(int nr_total_blocks, uint* d_off,
                                            uint* d_spine_sums)
{
  int bid = threadIdx.x + THREADS_PER_BLOCK * blockIdx.x;

  if (bid <= nr_total_blocks) {
    d_off[bid] = d_spine_sums[bid * CUDA_BND_STRIDE + 0];
  }
}

template <typename CudaMparticles>
void cuda_bndp<CudaMparticles, dim_yz>::update_offsets(CudaMparticles* cmprts)
{
  int dimGrid = (n_blocks + 1 + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;

  mprts_update_offsets<<<dimGrid, THREADS_PER_BLOCK>>>(
    n_blocks, cmprts->by_block_.d_off.data().get(), d_spine_sums.data().get());
  cuda_sync_if_enabled();
}

template <typename CudaMparticles>
void cuda_bndp<CudaMparticles, dim_yz>::update_offsets_gold(
  CudaMparticles* cmprts)
{
  thrust::host_vector<uint> h_spine_sums(
    d_spine_sums.data(), d_spine_sums.data() + 1 + n_blocks * (10 + 1));
  thrust::host_vector<uint> h_off(n_blocks + 1);

  for (int bid = 0; bid <= n_blocks; bid++) {
    h_off[bid] = h_spine_sums[bid * 10];
  }

  thrust::copy(h_off.begin(), h_off.end(), cmprts->by_block_.d_off.begin());
}

#endif

// ----------------------------------------------------------------------
// convert_and_copy_to_dev

template <typename CudaMparticles, typename DIM>
uint cuda_bndp<CudaMparticles, DIM>::convert_and_copy_to_dev(
  CudaMparticles* cmprts)
{
  uint n_recv = 0;
  for (int p = 0; p < n_patches(); p++) {
    n_recv += bufs[p].size();
  }

  HMparticlesCudaStorage h_bnd(n_recv);
  thrust::host_vector<uint> h_bnd_idx(n_recv);
  // thrust::host_vector<uint> h_bnd_off(n_recv);

  thrust::host_vector<uint> h_bnd_cnt(n_blocks, 0);

  uint off = 0;
  for (int p = 0; p < n_patches(); p++) {
    int n_recv = bufs[p].size();
    n_recvs[p] = n_recv;

    for (int n = 0; n < n_recv; n++) {
      const auto& prt = bufs[p][n];

      h_bnd.store(prt, n + off);
      checkInPatchMod(h_bnd[n + off].x);
      uint b = blockIndex(h_bnd[n + off], p);
      assert(b < n_blocks);
      h_bnd_idx[n + off] = b;
      // h_bnd_off[n + off] = h_bnd_cnt[b]++;
    }
    off += n_recv;
  }

  cmprts->resize(cmprts->n_prts + n_recv);

#ifdef PSC_HAVE_RMM
  thrust::copy(rmm::exec_policy(), h_bnd.begin(), h_bnd.end(),
               cmprts->storage.begin() + cmprts->n_prts);
#else
  thrust::copy(h_bnd.begin(), h_bnd.end(),
               cmprts->storage.begin() + cmprts->n_prts);
#endif
  thrust::copy(h_bnd_idx.begin(), h_bnd_idx.end(),
               cmprts->by_block_.d_idx.begin() + cmprts->n_prts);
  // // slight abuse of the now unused last part of spine_cnts
  // thrust::copy(h_bnd_cnt.begin(), h_bnd_cnt.end(), d_spine_cnts.begin() + 10
  // * n_blocks);

  // d_bnd_off.resize(n_recv);
  // thrust::copy(h_bnd_off.begin(), h_bnd_off.end(), d_bnd_off.begin());

  return n_recv;
}

template <typename CudaMparticles, typename DIM>
void cuda_bndp<CudaMparticles, DIM>::post(CudaMparticles* _cmprts)
{
  auto& cmprts = *_cmprts;
  auto& d_bidx = cmprts.by_block_.d_idx;

  auto n_prts_send = d_bidx.size() - cmprts.n_prts;
  auto n_prts_recv = convert_and_copy_to_dev(&cmprts);
  cmprts.n_prts += n_prts_recv;

  thrust::sequence(cmprts.by_block_.d_id.begin(), cmprts.by_block_.d_id.end());
#ifdef PSC_HAVE_RMM
  thrust::stable_sort_by_key(rmm::exec_policy(), d_bidx.begin(), d_bidx.end(),
                             cmprts.by_block_.d_id.begin());
#else
  thrust::stable_sort_by_key(d_bidx.begin(), d_bidx.end(),
                             cmprts.by_block_.d_id.begin());
#endif

  // drop the previously sent particles, which have been sorted to the end of
  // the array, now
  cmprts.n_prts -= n_prts_send;
  // FIXME, this is evil, but done in yz case, too: even though the arrays are
  // still have 6 elements because they contain previously sent particles
  // (a.k.a., gaps), the actual # particles is only 4, but they are pointed to
  // correctly by d_id
  // cmprts.resize(cmprts.n_prts);

  // find offsets
  thrust::counting_iterator<uint> search_begin(0);
  thrust::upper_bound(d_bidx.begin(), d_bidx.end(), search_begin,
                      search_begin + cmprts.n_blocks,
                      cmprts.by_block_.d_off.begin() + 1);
  // d_off[0] was set to zero during d_off initialization

  cmprts.need_reorder = true;
}

template struct cuda_bndp<cuda_mparticles<BS144>, dim_yz>;
template struct cuda_bndp<cuda_mparticles<BS444>, dim_xyz>;
