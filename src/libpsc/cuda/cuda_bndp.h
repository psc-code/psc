
#ifndef CUDA_BNDP_H
#define CUDA_BNDP_H

#include "cuda_mparticles_indexer.h"
#include "cuda_mparticles.h"
#include "psc_particles_cuda.h"
#include "ddc_particles.hxx"
#include "cuda_bits.h"

#include <thrust/device_vector.h>
#include <thrust/partition.h>

// ======================================================================
// bnd

// ----------------------------------------------------------------------
// cuda_bnd

struct cuda_bnd {
  std::vector<particle_cuda_t> buf;
  int n_recv;
  int n_send;
};

// ----------------------------------------------------------------------
// cuda_bndp

template<typename CudaMparticles, typename DIM>
struct cuda_bndp : cuda_mparticles_indexer<typename CudaMparticles::BS>
{
  using BS = typename CudaMparticles::BS;
  using buf_t = std::vector<particle_cuda_t>;

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

  std::vector<cuda_bnd> bpatch;
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
    bpatch.resize(n_patches);
    bufs_.reserve(n_patches);
    for (int p = 0; p < n_patches; p++) {
      bufs_.push_back(&bpatch[p].buf);
    }
  }

  // ----------------------------------------------------------------------
  // prep

  void prep(CudaMparticles* _cmprts)
  {
    auto& cmprts = *_cmprts;
    auto& d_bidx = cmprts.by_block_.d_idx;
    
    auto begin = thrust::make_zip_iterator(thrust::make_tuple(d_bidx.begin(), cmprts.d_xi4.begin(), cmprts.d_pxi4.begin()));
    auto end = thrust::make_zip_iterator(thrust::make_tuple(d_bidx.end(), cmprts.d_xi4.end(), cmprts.d_pxi4.end()));
    auto oob = thrust::stable_partition(begin, end, is_inside(cmprts.n_blocks));

    n_prts_send = end - oob;
    cmprts.n_prts -= n_prts_send;

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
    thrust::host_vector<float4> h_bnd_xi4(n_prts_send);
    thrust::host_vector<float4> h_bnd_pxi4(n_prts_send);
    thrust::host_vector<uint> h_bidx(n_prts_send);
    
    assert(cmprts->d_xi4.begin() + n_prts + n_prts_send == cmprts->d_xi4.end());
    
    thrust::copy(cmprts->d_xi4.begin()  + n_prts, cmprts->d_xi4.end(), h_bnd_xi4.begin());
    thrust::copy(cmprts->d_pxi4.begin() + n_prts, cmprts->d_pxi4.end(), h_bnd_pxi4.begin());
    thrust::copy(cmprts->by_block_.d_idx.begin() + n_prts, cmprts->by_block_.d_idx.end(),
		 h_bidx.begin());

    for (int p = 0; p < n_patches; p++) {
      bpatch[p].buf.clear();
      bpatch[p].n_send = 0;
    }
    for (int n = 0; n < n_prts_send; n++) {
      int kind = cuda_float_as_int(h_bnd_xi4[n].w);
      auto prt = particle_cuda_t{{h_bnd_xi4[n].x, h_bnd_xi4[n].y, h_bnd_xi4[n].z},
				 {h_bnd_pxi4[n].x, h_bnd_pxi4[n].y, h_bnd_pxi4[n].z},
				 h_bnd_pxi4[n].w / float(cmprts->grid_.kinds[kind].q),
				 kind};

      int p = h_bidx[n] - cmprts->n_blocks;
      auto& buf = bpatch[p].buf;
      bpatch[p].buf.push_back(prt);
      bpatch[p].n_send++;
    }
  }

  uint convert_and_copy_to_dev(CudaMparticles *cmprts);

  struct is_inside
  {
    is_inside(int n_blocks) : n_blocks_(n_blocks) {}
    
    __host__ __device__
    bool operator()(thrust::tuple<uint, float4, float4> tup)
    {
      uint bidx = thrust::get<0>(tup);
      return bidx < n_blocks_;
    }
    
    int n_blocks_;
  };

  std::vector<cuda_bnd> bpatch;
  std::vector<buf_t*> bufs_;
  uint n_prts_send;
};
  
#endif

