
#ifndef CUDA_BNDP_H
#define CUDA_BNDP_H

#include "cuda_bits.h"
#include "cuda_mparticles.hxx"
#include "cuda_mparticles_indexer.h"
#include "ddc_particles.hxx"
#include "mparticles_cuda.hxx"

#include <thrust/device_vector.h>
#include <thrust/partition.h>

// ----------------------------------------------------------------------
// cuda_bndp

template <typename CudaMparticles, typename DIM>
struct cuda_bndp : cuda_mparticles_indexer<typename CudaMparticles::BS>
{
  using BS = typename CudaMparticles::BS;
  using BndBuffer = std::vector<typename CudaMparticles::Particle>;
  using BndBuffers = typename MparticlesCuda<BS>::BndBuffers;

  using cuda_mparticles_indexer<BS>::n_blocks;
  using cuda_mparticles_indexer<BS>::n_blocks_per_patch;
  using cuda_mparticles_indexer<BS>::n_patches;
  using cuda_mparticles_indexer<BS>::checkInPatchMod;
  using cuda_mparticles_indexer<BS>::blockIndex;
  using cuda_mparticles_indexer<BS>::b_mx;

  cuda_bndp(const Grid_t& grid) : cuda_mparticles_indexer<BS>{grid}
  {
    bufs.resize(n_patches());
    n_sends.resize(n_patches());
    n_recvs.resize(n_patches());
  }

  // ----------------------------------------------------------------------
  // prep

  BndBuffers& prep(CudaMparticles* _cmprts)
  {
    auto& cmprts = *_cmprts;
    auto& d_bidx = cmprts.by_block_.d_idx;

#ifdef PSC_HAVE_RMM
    auto oob = thrust::count_if(rmm::exec_policy(), d_bidx.begin(),
                                d_bidx.end(), is_outside(cmprts.n_blocks));
#else
    auto oob = thrust::count_if(d_bidx.begin(), d_bidx.end(),
                                is_outside(cmprts.n_blocks));
#endif
    auto sz = d_bidx.size();
    assert(cmprts.storage.size() == sz);
    assert(cmprts.n_prts == sz);
    mem_bndp -= allocated_bytes(d_bidx);
    d_bidx.resize(sz + oob);
    mem_bndp += allocated_bytes(d_bidx);
    cmprts.storage.resize(sz + oob);

    auto begin = thrust::make_zip_iterator(
      thrust::make_tuple(d_bidx.begin(), cmprts.storage.begin()));
    auto end = begin + sz;

#ifdef PSC_HAVE_RMM
    auto oob_end = thrust::copy_if(rmm::exec_policy(), begin, end, begin + sz,
                                   is_outside(cmprts.n_blocks));
#else
    auto oob_end =
      thrust::copy_if(begin, end, begin + sz, is_outside(cmprts.n_blocks));
#endif
    assert(oob_end == begin + sz + oob);

    n_prts_send = oob;

    copy_from_dev_and_convert(&cmprts, n_prts_send);

    return bufs;
  }

  // ----------------------------------------------------------------------
  // post

  void post(CudaMparticles* _cmprts);

  // ----------------------------------------------------------------------
  // copy_from_dev_and_convert

  void copy_from_dev_and_convert(CudaMparticles* cmprts, uint n_prts_send)
  {
    uint n_prts = cmprts->n_prts;
    HMparticlesCudaStorage h_bnd_storage(n_prts_send);
    thrust::host_vector<uint> h_bidx(n_prts_send);

    assert(cmprts->storage.size() == n_prts + n_prts_send);

#ifdef PSC_HAVE_RMM
    thrust::copy(rmm::exec_policy(), cmprts->storage.begin() + n_prts,
                 cmprts->storage.end(), h_bnd_storage.begin());
#else
    thrust::copy(cmprts->storage.begin() + n_prts, cmprts->storage.end(),
                 h_bnd_storage.begin());
#endif
    thrust::copy(cmprts->by_block_.d_idx.begin() + n_prts,
                 cmprts->by_block_.d_idx.end(), h_bidx.begin());

    for (int p = 0; p < n_patches(); p++) {
      bufs[p].clear();
      n_sends[p] = 0;
    }
    for (int n = 0; n < n_prts_send; n++) {
      auto prt = h_bnd_storage[n];
      int p = h_bidx[n] - cmprts->n_blocks;
      bufs[p].push_back(prt);
      n_sends[p]++;
    }
  }

  uint convert_and_copy_to_dev(CudaMparticles* cmprts);

  struct is_outside
  {
    is_outside(int n_blocks) : n_blocks_(n_blocks) {}

    __host__ __device__ bool operator()(uint bidx) { return bidx >= n_blocks_; }

    __host__ __device__ bool operator()(
      thrust::tuple<uint, thrust::tuple<float4, float4>> tup)
    {
      uint bidx = thrust::get<0>(tup);
      return (*this)(bidx);
    }

  private:
    int n_blocks_;
  };

  std::vector<BndBuffer> bufs;
  std::vector<int> n_sends;
  std::vector<int> n_recvs;
  BndBuffers bufs_;
  uint n_prts_send;
};

#ifdef CUDA_BNDP_DIM_YZ_SPECIAL

// ----------------------------------------------------------------------
// specialized for dim_yz

template <typename CudaMparticles>
struct cuda_bndp<CudaMparticles, dim_yz>
  : cuda_mparticles_indexer<typename CudaMparticles::BS>
{
  using BS = typename CudaMparticles::BS;
  using BndBuffer = std::vector<typename CudaMparticles::Particle>;
  using BndBuffers = typename CudaMparticles::BndBuffers;

  using cuda_mparticles_indexer<BS>::n_blocks;
  using cuda_mparticles_indexer<BS>::n_blocks_per_patch;
  using cuda_mparticles_indexer<BS>::n_patches;
  using cuda_mparticles_indexer<BS>::checkInPatchMod;
  using cuda_mparticles_indexer<BS>::blockIndex;
  using cuda_mparticles_indexer<BS>::b_mx;

  cuda_bndp(const Grid_t& grid);

  BndBuffers& prep(CudaMparticles* cmprts);
  void post(CudaMparticles* cmprts);

  // pieces for prep
  void spine_reduce(CudaMparticles* cmprts);
  uint find_n_send(CudaMparticles* cmprts);
  void scan_send_buf_total(CudaMparticles* cmprts, uint n_prts_send);
  void reorder_send_by_id(CudaMparticles* cmprts, uint n_prts_send);
  void reorder_send_buf_total(CudaMparticles* cmprts, uint n_prts_send);
  void copy_from_dev_and_convert(CudaMparticles* cmprts, uint n_prts_send);

  // pieces for post
  uint convert_and_copy_to_dev(CudaMparticles* cmprts);
  void sort_pairs_device(CudaMparticles* cmprts, uint n_prts_recv);
  void count_received(CudaMparticles* cmprts);
  void scan_scatter_received(CudaMparticles* cmprts, uint n_prts_recv);
  void update_offsets(CudaMparticles* cmprts);

  // gold
  void spine_reduce_gold(CudaMparticles* cmprts);
  void scan_send_buf_total_gold(CudaMparticles* cmprts, uint n_prts_send);
  void reorder_send_by_id_gold(CudaMparticles* cmprts, uint n_prts_send);
  void sort_pairs_gold(CudaMparticles* cmprts, uint n_prts_recv);
  void count_received_gold(CudaMparticles* cmprts);
  void scan_scatter_received_gold(CudaMparticles* cmprts, uint n_prts_recv);
  void update_offsets_gold(CudaMparticles* cmprts);

  psc::device_vector<uint> d_spine_cnts;
  psc::device_vector<uint> d_spine_sums;
  uint n_prts_send;
  psc::device_vector<uint> d_bnd_off;

  psc::device_vector<uint>
    d_sums; // FIXME, should go away (only used in some gold stuff)

  BndBuffers bufs;
  std::vector<int> n_sends;
  std::vector<int> n_recvs;
};

#endif

#endif
