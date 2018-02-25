
#ifndef CUDA_BNDP_H
#define CUDA_BNDP_H

#include "cuda_mparticles_indexer.h"
#include "psc_particles_cuda.h"
#include "ddc_particles.hxx"

#include <thrust/device_vector.h>

// ======================================================================
// bnd

#define CUDA_BND_S_NEW (9)
#define CUDA_BND_S_OOB (10)
#define CUDA_BND_STRIDE (10)

// ----------------------------------------------------------------------
// cuda_bnd

struct cuda_bnd {
  psc_particle_cuda_buf_t buf;
  int n_recv;
  int n_send;
};

// ----------------------------------------------------------------------
// cuda_bndp

struct cuda_bndp : cuda_mparticles_indexer
{
  using mparticles_t = PscMparticlesCuda;
  using ddcp_t = ddc_particles<mparticles_t>;
  using ddcp_patch = typename ddcp_t::patch;

  cuda_bndp(const Grid_t& grid);
  
  void prep(ddcp_t* ddcp, cuda_mparticles* cmprts);
  void post(ddcp_t* ddcp, cuda_mparticles* cmprts);

  // pieces for prep
  void spine_reduce(cuda_mparticles *cmprts);
  uint find_n_send(cuda_mparticles *cmprts);
  void scan_send_buf_total(cuda_mparticles *cmprts, uint n_prts_send);
  void reorder_send_by_id(cuda_mparticles *cmprts, uint n_prts_send);
  void reorder_send_buf_total(cuda_mparticles *cmprts, uint n_prts_send);
  void copy_from_dev_and_convert(cuda_mparticles *cmprts, uint n_prts_send);

  // pieces for post
  uint convert_and_copy_to_dev(cuda_mparticles *cmprts);
  void sort_pairs_device(cuda_mparticles *cmprts, uint n_prts_recv);
  void count_received(cuda_mparticles *cmprts);
  void scan_scatter_received(cuda_mparticles *cmprts, uint n_prts_recv);
  void update_offsets(cuda_mparticles *cmprts);

  // gold
  void spine_reduce_gold(cuda_mparticles *cmprts);
  void scan_send_buf_total_gold(cuda_mparticles *cmprts, uint n_prts_send);
  void reorder_send_by_id_gold(cuda_mparticles *cmprts, uint n_prts_send);
  void sort_pairs_gold(cuda_mparticles *cmprts, uint n_prts_recv);
  void count_received_gold(cuda_mparticles *cmprts);
  void scan_scatter_received_gold(cuda_mparticles *cmprts, uint n_prts_recv);
  void update_offsets_gold(cuda_mparticles *cmprts);

  thrust::device_vector<uint> d_spine_cnts;
  thrust::device_vector<uint> d_spine_sums;
  uint n_prts_send;
  thrust::device_vector<uint> d_bnd_off;

  thrust::device_vector<uint> d_sums; // FIXME, should go away (only used in some gold stuff)

  std::vector<cuda_bnd> bpatch;
};

#endif

