
#ifndef CUDA_MPARTICLES_H
#define CUDA_MPARTICLES_H

#include "cuda_iface.h"

#include "grid.hxx"
#include "particles.hxx"

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
// cuda_mparticles_bnd

struct cuda_mparticles;

struct cuda_mparticles_bnd
{
  ~cuda_mparticles_bnd();

  void setup(cuda_mparticles *cmprts);
  void free_particle_mem();
  void reserve_all(cuda_mparticles *cmprts);

  void scan_send_buf_total(cuda_mparticles *cmprts);
  void scan_send_buf_total_gold(cuda_mparticles *cmprts);
  void spine_reduce(cuda_mparticles *cmprts);
  void sort_pairs_device(cuda_mparticles *cmprts);

  void spine_reduce_gold(cuda_mparticles *cmprts);
  void sort_pairs_gold(cuda_mparticles *cmprts);

  void reorder_send_by_id(struct cuda_mparticles *cmprts);
  void find_n_send(cuda_mparticles *cmprts);
  void copy_from_dev_and_convert(cuda_mparticles *cmprts);
  void reorder_send_by_id_gold(cuda_mparticles *cmprts);
  void convert_and_copy_to_dev(cuda_mparticles *cmprts);
  void sort(cuda_mparticles *cmprts, int *n_prts_by_patch);
  void update_offsets(cuda_mparticles *cmprts);
  void update_offsets_gold(cuda_mparticles *cmprts);
  void count_received(cuda_mparticles *cmprts);
  void count_received_gold(cuda_mparticles *cmprts);
  void scan_scatter_received(cuda_mparticles *cmprts);
  void scan_scatter_received_gold(cuda_mparticles *cmprts);
  void reorder_send_buf_total(cuda_mparticles *cmprts);
  
public:
  thrust::device_vector<uint> d_alt_bidx;
  thrust::device_vector<uint> d_sums; // FIXME, too many arrays, consolidation would be good

  uint n_prts_send;
  uint n_prts_recv;

  thrust::device_vector<uint> d_spine_cnts;
  thrust::device_vector<uint> d_spine_sums;

  struct cuda_bnd *bpatch;
};

// ----------------------------------------------------------------------
// cuda_mparticles

struct cuda_mparticles : cuda_mparticles_bnd
{
public:
  using particle_t = particle_cuda_t;
  using real_t = particle_t::real_t;
  using Real3 = Vec3<real_t>;

  cuda_mparticles(const Grid_t& grid, const Int3& bs);
  cuda_mparticles(const cuda_mparticles&) = delete;
  ~cuda_mparticles();

  void free_particle_mem();

  void reserve_all(const uint *n_prts_by_patch);
  void get_size_all(uint *n_prts_by_patch);
  void resize_all(const uint *n_prts_by_patch);
  uint get_n_prts();
  void set_particles(uint n_prts, uint off,
		     void (*get_particle)(cuda_mparticles_prt *prt, int n, void *ctx),
		     void *ctx);
  void get_particles(uint n_prts, uint off,
		     void (*put_particle)(cuda_mparticles_prt *, int, void *),
		     void *ctx);
  void setup_internals();
  void inject(cuda_mparticles_prt *buf, uint *buf_n_by_patch);
  const particle_cuda_real_t *patch_get_b_dxi(int p);
  const int *patch_get_b_mx(int p);

  psc_particle_cuda_buf_t *bnd_get_buffer(int p);
  void bnd_prep();
  void bnd_post();
  
  void dump();
  void dump_by_patch(uint *n_prts_by_patch);

public:
  void to_device(float_4 *xi4, float_4 *pxi4,
		 uint n_prts, uint off);
  void from_device(float_4 *xi4, float_4 *pxi4,
		   uint n_prts, uint off);
  
  void find_block_indices_ids();
  void reorder_and_offsets();
  void check_in_patch_unordered_slow(uint *nr_prts_by_patch);
  void check_bidx_id_unordered_slow(uint *n_prts_by_patch);
  void check_ordered();
  void check_ordered_slow();
  
public:
  // per particle
  float4 *d_xi4 = {};             // current particle data
  float4 *d_pxi4 = {};
  float4 *d_alt_xi4  = {};        // storage for out-of-place reordering of particle data
  float4 *d_alt_pxi4 = {};
  uint *d_bidx = {};              // block index (incl patch) per particle
  uint *d_id = {};                // particle id for sorting

  // per block
  thrust::device_vector<uint> d_off;      // particles per block
                                  // are at indices [offsets[block] .. offsets[block+1]-1[

  uint n_prts = {};               // total # of particles across all patches
  uint n_alloced = {};            // size of particle-related arrays as allocated
  uint n_patches;                 // # of patches
  uint n_blocks_per_patch;        // number of blocks per patch
  uint n_blocks;                  // number of blocks in all patches in mprts

  Int3 b_mx;                      // number of blocks per direction in each patch
  Int3 bs;
  Real3 b_dxi;                    // inverse of block size (in actual length units)
  std::vector<Real3> xb_by_patch; // lower left corner for each patch

  bool need_reorder;              // particles haven't yet been put into their sorted order

public:
  const Grid_t& grid_;
};

void cuda_mparticles_swap_alt(struct cuda_mparticles *cmprts);
void cuda_mparticles_reorder(struct cuda_mparticles *cmprts);

#endif
