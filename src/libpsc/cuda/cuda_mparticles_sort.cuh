
#pragma once

#include "cuda_bits.h"
#include "rng_state.cuh"

#include <thrust/binary_search.h>
#include <thrust/device_vector.h>
#include <thrust/random.h>
#include <thrust/sort.h>

#include <curand_kernel.h>

template <typename BS>
struct cuda_mparticles;

template <typename BS>
struct DMparticlesCuda;

#define THREADS_PER_BLOCK 256

// ----------------------------------------------------------------------
// find_cell_indices_ids

template <typename BS>
__global__ static void k_find_cell_indices_ids(DMparticlesCuda<BS> dmprts,
                                               uint* d_cidx, uint* d_id,
                                               int n_patches,
                                               int n_blocks_per_patch)
{
  int n = threadIdx.x + THREADS_PER_BLOCK * blockIdx.x;

  for (int p = 0; p < n_patches; p++) {
    uint off = dmprts.off_[p * n_blocks_per_patch];
    uint n_prts = dmprts.off_[(p + 1) * n_blocks_per_patch] - off;
    if (n < n_prts) {
      float4 xi4 = dmprts.storage.xi4[n + off];
      d_cidx[n + off] = dmprts.validCellIndex(xi4, p);
      d_id[n + off] = n + off;
    }
  }
}

template <typename BS>
inline void find_cell_indices_ids(cuda_mparticles<BS>& cmprts,
                                  thrust::device_vector<uint>& d_cidx,
                                  thrust::device_vector<uint>& d_id)
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

  dim3 dimGrid((max_n_prts + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK);
  dim3 dimBlock(THREADS_PER_BLOCK);

  k_find_cell_indices_ids<BS>
    <<<dimGrid, dimBlock>>>(cmprts, d_cidx.data().get(), d_id.data().get(),
                            cmprts.n_patches(), cmprts.n_blocks_per_patch);
  cuda_sync_if_enabled();
}

// ----------------------------------------------------------------------
// find_random_cell_indices_ids

template <typename BS>
__global__ static void k_find_random_cell_indices_ids(
  DMparticlesCuda<BS> dmprts, float* d_random_idx, uint* d_id, int n_patches,
  int n_blocks_per_patch, RngStateCuda::Device rng_state)
{
  int n = threadIdx.x + THREADS_PER_BLOCK * blockIdx.x;

  auto rng = rng_state[n];

  for (int p = 0; p < n_patches; p++) {
    uint off = dmprts.off_[p * n_blocks_per_patch];
    uint n_prts = dmprts.off_[(p + 1) * n_blocks_per_patch] - off;
    if (n < n_prts) {
      float4 xi4 = dmprts.storage.xi4[n + off];
      d_random_idx[n + off] = dmprts.validCellIndex(xi4, p) + rng.uniform();
      d_id[n + off] = n + off;
    }
  }

  rng_state[n] = rng;
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
      float4 xi4 = dmprts.storage.xi4[n + off];
      d_bidx[n + off] = dmprts.blockIndex(xi4, p);
      d_id[n + off] = n + off;
    }
  }
}

template <typename BS>
inline void find_block_indices_ids(cuda_mparticles<BS>& cmprts,
                                   thrust::device_vector<uint>& d_idx,
                                   thrust::device_vector<uint>& d_id)
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

  template <typename BS>
  void find_indices_ids(cuda_mparticles<BS>& cmprts)
  {
    d_idx.resize(cmprts.n_prts);
    d_id.resize(cmprts.n_prts);
    find_cell_indices_ids(cmprts, d_idx, d_id);
  }

  void stable_sort_cidx()
  {
    thrust::stable_sort_by_key(d_idx.begin(), d_idx.end(), d_id.begin());
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
  thrust::device_vector<uint> d_idx; // cell index (incl patch) per particle
  thrust::device_vector<uint> d_id;  // particle id used for reordering
  thrust::device_vector<uint>
    d_off; // particles per cell
           // are at indices [offsets[cell] .. offsets[cell+1][
};

// ======================================================================
// cuda_mparticles_randomize_sort
//
// by cell

struct cuda_mparticles_randomize_sort
{
  cuda_mparticles_randomize_sort(uint n_cells) : d_off(n_cells + 1) {}

  template <typename BS>
  void find_indices_ids(cuda_mparticles<BS>& cmprts)
  {
    d_id.resize(cmprts.n_prts);
    d_random_idx.resize(cmprts.n_prts);

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

    if (max_n_prts > rng_state_.size()) {
      rng_state_.resize(max_n_prts);
    }
    
    dim3 dimGrid((max_n_prts + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK);
    dim3 dimBlock(THREADS_PER_BLOCK);

    k_find_random_cell_indices_ids<BS><<<dimGrid, dimBlock>>>(
      cmprts, d_random_idx.data().get(), d_id.data().get(), cmprts.n_patches(),
      cmprts.n_blocks_per_patch, rng_state_);
    cuda_sync_if_enabled();
  }

  void sort()
  {
    thrust::sort_by_key(d_random_idx.begin(), d_random_idx.end(), d_id.begin());
  }

  void find_offsets()
  {
    int n_cells = d_off.size() - 1;
    thrust::counting_iterator<uint> search_begin(0);
    thrust::lower_bound(d_random_idx.begin(), d_random_idx.end(), search_begin,
                        search_begin + n_cells, d_off.begin());
    d_off[n_cells] = d_random_idx.size();
  }

public:
  thrust::device_vector<float> d_random_idx; // randomized cell index
  thrust::device_vector<uint> d_id;          // particle id used for reordering
  thrust::device_vector<uint>
    d_off; // particles per cell
           // are at indices [offsets[cell] .. offsets[cell+1][
  RngStateCuda rng_state_;
};

// ======================================================================
// cuda_mparticles_sort_by_block
//
// by block

struct cuda_mparticles_sort_by_block
{
  cuda_mparticles_sort_by_block(uint n_blocks) : d_off(n_blocks + 1) {}

  template <typename BS>
  void find_indices_ids(cuda_mparticles<BS>& cmprts)
  {
    d_idx.resize(cmprts.n_prts);
    d_id.resize(cmprts.n_prts);
    find_block_indices_ids(cmprts, d_idx, d_id);
  }

  void stable_sort()
  {
    thrust::stable_sort_by_key(d_idx.begin(), d_idx.end(), d_id.begin());
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

public:
  thrust::device_vector<uint> d_idx; // block index (incl patch) per particle
  thrust::device_vector<uint> d_id;  // particle id used for reordering
  thrust::device_vector<uint>
    d_off; // particles per cell
           // are at indices [offsets[block] .. offsets[block+1][
};

#undef THREADS_PER_BLOCK

