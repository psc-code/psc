
#pragma once

#include "cuda_bits.h"
#include "rng_state.hxx"
#include "bs.hxx"
#include "cuda_base.hxx"

#include <thrust/binary_search.h>
#include <thrust/device_vector.h>
#include <thrust/random.h>
#include <thrust/sort.h>

#include <curand_kernel.h>

extern std::size_t mem_randomize_sort;
extern std::size_t mem_sort_by_block;

template <typename BS>
struct cuda_mparticles;

template <typename BS>
struct DMparticlesCuda;

#define THREADS_PER_BLOCK 128

namespace detail
{
template <typename BS>
struct bs_to_dim;

template <>
struct bs_to_dim<BS144>
{
  using type = dim_yz;
};

template <>
struct bs_to_dim<BS444>
{
  using type = dim_xyz;
};
} // namespace detail

template <typename BS>
using bs_to_dim_t = typename detail::bs_to_dim<BS>::type;

// ----------------------------------------------------------------------
// find_cell_indices_ids

template <typename BS, typename Block>
__global__ static void k_find_cell_indices_ids(DMparticlesCuda<BS> dmprts,
                                               uint* d_cidx, uint* d_id,
                                               int n_blocks)
{
  int bid = blockIdx.x;

  Block current_block;
  for (; bid < n_blocks; bid += gridDim.x) {
    current_block.init(dmprts, bid);

    int block_begin = dmprts.off_[current_block.bid];
    int block_end = dmprts.off_[current_block.bid + 1];
    for (int n : in_block_loop(block_begin, block_end)) {
      if (n < block_begin) {
        continue;
      }
      auto prt = dmprts.storage[n];
      d_cidx[n] = dmprts.validCellIndex(prt, current_block.p);
      d_id[n] = n;
    }
  }
}

template <typename CMPRTS>
inline void find_cell_indices_ids(CMPRTS& cmprts,
                                  psc::device_vector<uint>& d_cidx,
                                  psc::device_vector<uint>& d_id)
{
  if (cmprts.n_patches() == 0) {
    return;
  }

  using BS = typename CMPRTS::BS;
  using dim = bs_to_dim_t<BS>;
  using Block = BlockSimple2<BS, dim>;
  dim3 dimGrid = Block::dimGrid(cmprts);

  ::k_find_cell_indices_ids<BS, Block><<<dimGrid, THREADS_PER_BLOCK>>>(
    cmprts, d_cidx.data().get(), d_id.data().get(), cmprts.n_blocks);
}

// ----------------------------------------------------------------------
// find_random_cell_indices_ids

template <typename BS, typename Block>
__global__ static void k_find_random_cell_indices_ids(
  DMparticlesCuda<BS> dmprts, double* d_random_idx, uint* d_id,
  RngStateCuda::Device rng_state, int n_blocks)
{
  Block current_block;
  int bid = blockIdx.x;
  int id = blockIdx.x * blockDim.x + threadIdx.x;
  auto rng = rng_state[id];

  for (; bid < n_blocks; bid += gridDim.x) {
    current_block.init(dmprts, bid);

    int block_begin = dmprts.off_[current_block.bid];
    int block_end = dmprts.off_[current_block.bid + 1];
    for (int n : in_block_loop(block_begin, block_end)) {
      if (n < block_begin) {
        continue;
      }

      auto prt = dmprts.storage[n];
      d_random_idx[n] =
        dmprts.validCellIndex(prt, current_block.p) + .5 * rng.uniform();
      d_id[n] = n;
    }
  }
  rng_state[id] = rng;
}

// ----------------------------------------------------------------------
// find_block_indices_ids

template <typename BS>
__global__ static void k_find_block_indices_ids(DMparticlesCuda<BS> dmprts,
                                                uint* d_bidx, uint* d_id,
                                                int n_patches,
                                                int n_blocks_per_patch)
{
  for (int p = 0; p < n_patches; p++) {
    uint off = dmprts.off_[p * n_blocks_per_patch];
    uint n_prts = dmprts.off_[(p + 1) * n_blocks_per_patch] - off;

    int n = threadIdx.x + blockDim.x * blockIdx.x;
    for (; n < n_prts; n += gridDim.x * blockDim.x) {
      auto prt = dmprts.storage[n + off];
      d_bidx[n + off] = dmprts.blockIndex(prt, p);
      d_id[n + off] = n + off;
    }
  }
}

template <typename BS>
inline void find_block_indices_ids(cuda_mparticles<BS>& cmprts,
                                   psc::device_vector<uint>& d_idx,
                                   psc::device_vector<uint>& d_id)
{
  if (cmprts.n_patches() == 0) {
    return;
  }

  // OPT: if we didn't need max_n_prts, we wouldn't have to get the
  // sizes / offsets at all, and it seems likely we could do a better
  // job here in general
  auto n_prts_by_patch = cmprts.sizeByPatch();

  int max_n_prts = 0;
  for (int p = 0; p < cmprts.n_patches(); p++) {
    if (n_prts_by_patch[p] > max_n_prts) {
      max_n_prts = n_prts_by_patch[p];
    }
  }

  if (max_n_prts == 0) {
    return;
  }

  int n_blocks = (max_n_prts + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
  // if (n_blocks > 32768) n_blocks = 32768;
  dim3 dimGrid(n_blocks);
  dim3 dimBlock(THREADS_PER_BLOCK);

  k_find_block_indices_ids<BS>
    <<<dimGrid, dimBlock>>>(cmprts, d_idx.data().get(), d_id.data().get(),
                            cmprts.n_patches(), cmprts.n_blocks_per_patch);
  cuda_sync_if_enabled();
}

// ======================================================================
// cuda_mparticles_sort
//
// by cell

struct cuda_mparticles_sort
{
  cuda_mparticles_sort(uint n_cells) : d_off(n_cells + 1) {}

  template <typename CMPRTS>
  void find_indices_ids(CMPRTS& cmprts)
  {
    d_idx.resize(cmprts.n_prts);
    d_id.resize(cmprts.n_prts);
    find_cell_indices_ids(cmprts, d_idx, d_id);
  }

  void stable_sort_cidx()
  {
#ifdef PSC_HAVE_RMM
    thrust::stable_sort_by_key(rmm::exec_policy(), d_idx.begin(), d_idx.end(),
                               d_id.begin());
#else
    thrust::stable_sort_by_key(d_idx.begin(), d_idx.end(), d_id.begin());
#endif
  }

  void find_offsets()
  {
    int n_cells = d_off.size() - 1;
    thrust::counting_iterator<uint> search_begin(0);
    thrust::upper_bound(d_idx.begin(), d_idx.end(), search_begin,
                        search_begin + n_cells, d_off.begin() + 1);
    // d_off[0] was set to zero during d_off initialization
  }

  template <typename BS>
  void reorder(cuda_mparticles<BS>& cmprts)
  {
    cmprts.reorder(d_id);
  }

public:
  psc::device_vector<uint> d_idx; // cell index (incl patch) per particle
  psc::device_vector<uint> d_id;  // particle id used for reordering
  psc::device_vector<uint>
    d_off; // particles per cell
           // are at indices [offsets[cell] .. offsets[cell+1][
};

// ======================================================================
// cuda_mparticles_randomize_sort
//
// by cell

struct cuda_mparticles_randomize_sort
{
  cuda_mparticles_randomize_sort() : rng_state_{get_rng_state()} {}

  ~cuda_mparticles_randomize_sort() {}

  template <typename CMPRTS>
  void operator()(CMPRTS& cmprts, psc::device_vector<uint>& d_off,
                  psc::device_vector<uint>& d_id)
  {
    if (cmprts.n_patches() == 0) {
      return;
    }

    assert(d_off.size() == cmprts.n_cells() + 1);
    assert(d_id.size() == cmprts.n_prts);

    psc::device_vector<double> d_random_idx(cmprts.n_prts);

    find_indices_ids(cmprts, d_random_idx, d_id);
    sort(d_random_idx, d_id);
    find_offsets(d_random_idx, d_off);
  }

  template <typename CMPRTS>
  void find_indices_ids(CMPRTS& cmprts,
                        psc::device_vector<double>& d_random_idx,
                        psc::device_vector<uint>& d_id)
  {
    using BS = typename CMPRTS::BS;
    using dim = bs_to_dim_t<BS>;
    using Block = BlockSimple2<BS, dim>;
    dim3 dimGrid = Block::dimGrid(cmprts);

    assert(d_random_idx.size() == cmprts.n_prts);
    assert(d_id.size() == cmprts.n_prts);

    if (dimGrid.x * THREADS_PER_BLOCK > rng_state_.size()) {
      rng_state_.resize(dimGrid.x * THREADS_PER_BLOCK);
    }

    ::k_find_random_cell_indices_ids<BS, Block><<<dimGrid, THREADS_PER_BLOCK>>>(
      cmprts, d_random_idx.data().get(), d_id.data().get(), rng_state_,
      cmprts.n_blocks);
    cuda_sync_if_enabled();
  }

  void sort(psc::device_vector<double>& d_random_idx,
            psc::device_vector<uint>& d_id)
  {
#ifdef PSC_HAVE_RMM
    thrust::sort_by_key(rmm::exec_policy(), d_random_idx.begin(),
                        d_random_idx.end(), d_id.begin());
#else
    thrust::sort_by_key(d_random_idx.begin(), d_random_idx.end(), d_id.begin());
#endif
  }

  void find_offsets(psc::device_vector<double>& d_random_idx,
                    psc::device_vector<uint>& d_off)
  {
    int n_cells = d_off.size() - 1;
    thrust::counting_iterator<uint> search_begin(0);
    thrust::lower_bound(d_random_idx.begin(), d_random_idx.end(), search_begin,
                        search_begin + n_cells, d_off.begin());
    d_off[n_cells] = d_random_idx.size();
  }

private:
  RngStateCuda& rng_state_;
};

// ======================================================================
// cuda_mparticles_sort_by_block
//
// by block

struct cuda_mparticles_sort_by_block
{
  cuda_mparticles_sort_by_block(uint n_blocks) : d_off(n_blocks + 1)
  {
    mem_sort_by_block += allocated_bytes(d_off);
  }

  ~cuda_mparticles_sort_by_block()
  {
    mem_sort_by_block -= allocated_bytes(d_off);
  }

  template <typename BS>
  void find_indices_ids(cuda_mparticles<BS>& cmprts)
  {
    mem_sort_by_block -= allocated_bytes(d_idx);
    d_idx.resize(cmprts.n_prts);
    mem_sort_by_block += allocated_bytes(d_idx);

    mem_sort_by_block -= allocated_bytes(d_id);
    d_id.resize(cmprts.n_prts);
    mem_sort_by_block += allocated_bytes(d_id);

    find_block_indices_ids(cmprts, d_idx, d_id);
  }

  void stable_sort()
  {
#ifdef PSC_HAVE_RMM
    thrust::stable_sort_by_key(rmm::exec_policy(), d_idx.begin(), d_idx.end(),
                               d_id.begin());
#else
    thrust::stable_sort_by_key(d_idx.begin(), d_idx.end(), d_id.begin());
#endif
  }

  void find_offsets()
  {
    int n_blocks = d_off.size() - 1;
    thrust::counting_iterator<uint> search_begin(0);
    thrust::upper_bound(d_idx.begin(), d_idx.end(), search_begin,
                        search_begin + n_blocks, d_off.begin() + 1);
    // d_off[0] was set to zero during d_off initialization
  }

  template <typename BS>
  void reorder(cuda_mparticles<BS>& cmprts)
  {
    cmprts.reorder(d_id);
  }

  template <typename BS>
  void reorder_and_offsets(cuda_mparticles<BS>& cmprts)
  {
    cmprts.reorder_and_offsets(d_idx, d_id, d_off);
  }

  void clear() { thrust::fill(d_off.begin(), d_off.end(), 0); }

public:
  psc::device_vector<uint> d_idx; // block index (incl patch) per particle
  psc::device_vector<uint> d_id;  // particle id used for reordering
  psc::device_vector<uint>
    d_off; // particles per cell
           // are at indices [offsets[block] .. offsets[block+1][
};

#undef THREADS_PER_BLOCK
