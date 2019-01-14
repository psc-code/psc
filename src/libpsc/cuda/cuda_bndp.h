
#ifndef CUDA_BNDP_H
#define CUDA_BNDP_H

#include "cuda_mparticles_indexer.h"
#include "cuda_mparticles.cuh"
#include "mparticles_cuda.hxx"
#include "ddc_particles.hxx"
#include "cuda_bits.h"

#include <thrust/device_vector.h>
#include <thrust/partition.h>

// ======================================================================
// bnd

// ----------------------------------------------------------------------
// cuda_bnd

template<typename buf_t>
struct cuda_bnd {
  buf_t buf;
  int n_recv;
  int n_send;
};

// ----------------------------------------------------------------------
// cuda_bndp

template<typename CudaMparticles, typename DIM>
struct cuda_bndp : cuda_mparticles_indexer<typename CudaMparticles::BS>
{
  using BS = typename CudaMparticles::BS;
  using buf_t = std::vector<typename CudaMparticles::Particle>;

  using cuda_mparticles_indexer<BS>::n_blocks;
  using cuda_mparticles_indexer<BS>::n_blocks_per_patch;
  using cuda_mparticles_indexer<BS>::n_patches;
  using cuda_mparticles_indexer<BS>::checkInPatchMod;
  using cuda_mparticles_indexer<BS>::blockIndex;
  using cuda_mparticles_indexer<BS>::b_mx;

  cuda_bndp(const Grid_t& grid);
  
  void prep(CudaMparticles* cmprts);
  void post(CudaMparticles* cmprts);

  // pieces for prep
  void spine_reduce(CudaMparticles *cmprts);
  uint find_n_send(CudaMparticles *cmprts);
  void scan_send_buf_total(CudaMparticles *cmprts, uint n_prts_send);
  void reorder_send_by_id(CudaMparticles *cmprts, uint n_prts_send);
  void reorder_send_buf_total(CudaMparticles *cmprts, uint n_prts_send);
  void copy_from_dev_and_convert(CudaMparticles *cmprts, uint n_prts_send);

  // pieces for post
  uint convert_and_copy_to_dev(CudaMparticles *cmprts);
  void sort_pairs_device(CudaMparticles *cmprts, uint n_prts_recv);
  void count_received(CudaMparticles *cmprts);
  void scan_scatter_received(CudaMparticles *cmprts, uint n_prts_recv);
  void update_offsets(CudaMparticles *cmprts);

  // gold
  void spine_reduce_gold(CudaMparticles *cmprts);
  void scan_send_buf_total_gold(CudaMparticles *cmprts, uint n_prts_send);
  void reorder_send_by_id_gold(CudaMparticles *cmprts, uint n_prts_send);
  void sort_pairs_gold(CudaMparticles *cmprts, uint n_prts_recv);
  void count_received_gold(CudaMparticles *cmprts);
  void scan_scatter_received_gold(CudaMparticles *cmprts, uint n_prts_recv);
  void update_offsets_gold(CudaMparticles *cmprts);

  thrust::device_vector<uint> d_spine_cnts;
  thrust::device_vector<uint> d_spine_sums;
  uint n_prts_send;
  thrust::device_vector<uint> d_bnd_off;

  thrust::device_vector<uint> d_sums; // FIXME, should go away (only used in some gold stuff)

  std::vector<cuda_bnd<buf_t>> bpatch;
  std::vector<buf_t*> bufs_;
};

template<typename CudaMparticles>
struct cuda_bndp<CudaMparticles, dim_xyz> : cuda_mparticles_indexer<typename CudaMparticles::BS>
{
  using BS = typename CudaMparticles::BS;
  using buf_t = typename MparticlesCuda<BS>::buf_t;

  using cuda_mparticles_indexer<BS>::n_blocks;
  using cuda_mparticles_indexer<BS>::n_blocks_per_patch;
  using cuda_mparticles_indexer<BS>::n_patches;
  using cuda_mparticles_indexer<BS>::checkInPatchMod;
  using cuda_mparticles_indexer<BS>::blockIndex;
  using cuda_mparticles_indexer<BS>::b_mx;

  cuda_bndp(const Grid_t& grid)
    : cuda_mparticles_indexer<BS>{grid}
  {
    bpatch.resize(n_patches());
    bufs_.reserve(n_patches());
    for (int p = 0; p < n_patches(); p++) {
      bufs_.push_back(&bpatch[p].buf);
    }
  }

  // ----------------------------------------------------------------------
  // prep

  void prep(CudaMparticles* _cmprts)
  {
    auto& cmprts = *_cmprts;
    auto& d_bidx = cmprts.by_block_.d_idx;
    
    auto oob = thrust::count_if(d_bidx.begin(), d_bidx.end(), is_outside(cmprts.n_blocks));
    auto sz = d_bidx.size();
    assert(cmprts.storage.xi4.size() == sz);
    assert(cmprts.storage.pxi4.size() == sz);
    assert(cmprts.n_prts == sz);
    d_bidx.resize(sz + oob);
    cmprts.storage.xi4.resize(sz + oob);
    cmprts.storage.pxi4.resize(sz + oob);

    auto begin = thrust::make_zip_iterator(thrust::make_tuple(d_bidx.begin(), cmprts.storage.xi4.begin(), cmprts.storage.pxi4.begin()));
    auto end = begin + sz;
    
    auto oob_end = thrust::copy_if(begin, end, begin + sz, is_outside(cmprts.n_blocks));
    assert(oob_end == begin + sz + oob);

    n_prts_send = oob;

    copy_from_dev_and_convert(&cmprts, n_prts_send);
  }

  // ----------------------------------------------------------------------
  // post
  
  void post(CudaMparticles* _cmprts);

  // ----------------------------------------------------------------------
  // copy_from_dev_and_convert

  void copy_from_dev_and_convert(CudaMparticles *cmprts, uint n_prts_send)
  {
    uint n_prts = cmprts->n_prts;
    HMparticlesCudaStorage h_bnd_storage(n_prts_send);
    thrust::host_vector<uint> h_bidx(n_prts_send);
    
    assert(cmprts->storage.xi4.begin() + n_prts + n_prts_send == cmprts->storage.xi4.end());
    
    thrust::copy(cmprts->storage.xi4.begin()  + n_prts, cmprts->storage.xi4.end(), h_bnd_storage.xi4.begin());
    thrust::copy(cmprts->storage.pxi4.begin() + n_prts, cmprts->storage.pxi4.end(), h_bnd_storage.pxi4.begin());
    thrust::copy(cmprts->by_block_.d_idx.begin() + n_prts, cmprts->by_block_.d_idx.end(),
		 h_bidx.begin());

    for (int p = 0; p < n_patches(); p++) {
      bpatch[p].buf.clear();
      bpatch[p].n_send = 0;
    }
    for (int n = 0; n < n_prts_send; n++) {
      auto prt = h_bnd_storage.load(n);
      int p = h_bidx[n] - cmprts->n_blocks;
      auto& buf = bpatch[p].buf;
      bpatch[p].buf.push_back(prt);
      bpatch[p].n_send++;
    }
  }

  uint convert_and_copy_to_dev(CudaMparticles *cmprts);

  struct is_outside
  {
    is_outside(int n_blocks) : n_blocks_(n_blocks) {}
    
    __host__ __device__
    bool operator()(uint bidx) { return bidx >= n_blocks_; }
    
    __host__ __device__
    bool operator()(thrust::tuple<uint, float4, float4> tup)
    {
      uint bidx = thrust::get<0>(tup);
      return (*this)(bidx);
    }
    
  private:
    int n_blocks_;
  };

  std::vector<cuda_bnd<buf_t>> bpatch;
  std::vector<buf_t*> bufs_;
  uint n_prts_send;
};
  
#endif

