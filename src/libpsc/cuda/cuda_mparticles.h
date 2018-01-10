
#ifndef CUDA_MPARTICLES_H
#define CUDA_MPARTICLES_H

#include "cuda_iface.h"

#include "psc_particle_buf_cuda.h"

#ifndef __CUDACC__
struct float4 { float x; float y; float z; float w; };
#endif

// ----------------------------------------------------------------------
// float_3

typedef float float_3[3];

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

struct cuda_mparticles_bnd {
  unsigned int *d_alt_bidx;
  unsigned int *d_sums; // FIXME, too many arrays, consolidation would be good

  unsigned int n_prts_send;
  unsigned int n_prts_recv;

  unsigned int *d_bnd_spine_cnts;
  unsigned int *d_bnd_spine_sums;

  struct cuda_bnd *bpatch;
};

void cuda_mparticles_bnd_setup(struct cuda_mparticles *cmprts);
void cuda_mparticles_bnd_destroy(struct cuda_mparticles *cmprts);
void cuda_mparticles_bnd_reserve_all(struct cuda_mparticles *cmprts);
void cuda_mparticles_bnd_free_particle_mem(struct cuda_mparticles *cmprts);

void cuda_mparticles_scan_send_buf_total(struct cuda_mparticles *cmprts);
void cuda_mparticles_spine_reduce(struct cuda_mparticles *cmprts);
void cuda_mparticles_sort_pairs_device(struct cuda_mparticles *cmprts);

// ----------------------------------------------------------------------
// cuda_mparticles

struct cuda_mparticles
{
public:
  cuda_mparticles(mrc_json_t json);
  cuda_mparticles(const cuda_mparticles&) = delete;
  ~cuda_mparticles();
    
public:
  // per particle
  float4 *d_xi4, *d_pxi4;         // current particle data
  float4 *d_alt_xi4, *d_alt_pxi4; // storage for out-of-place reordering of particle data
  unsigned int *d_bidx;           // block index (incl patch) per particle
  unsigned int *d_id;             // particle id for sorting

  // per block
  unsigned int *d_off;            // particles per block
                                  // are at indices [offsets[block] .. offsets[block+1]-1[

  unsigned int n_prts;            // total # of particles across all patches
  unsigned int n_alloced;         // size of particle-related arrays as allocated
  unsigned int n_patches;         // # of patches
  unsigned int n_blocks_per_patch;// number of blocks per patch
  unsigned int n_blocks;          // number of blocks in all patches in mprts

  int ldims[3];                   // number of cells per direction in each patch
  int b_mx[3];                    // number of blocks per direction in each patch
  int bs[3];
  float dx[3];                    // cell size (in actual length units)
  float b_dxi[3];                 // inverse of block size (in actual length units)
  float_3 *xb_by_patch;           // lower left corner for each patch

  bool need_reorder;              // particles haven't yet been put into their sorted order

  float fnqs;                     // normalization constant
  float eta;                      // another normalization constant
  float dt;                       // timestep (FIXME, shouldn't really be here, but it's constant...)

  int n_kinds;
  float *kind_q;
  float *kind_m;

  struct cuda_mparticles_bnd bnd;
};

void cuda_mparticles_swap_alt(struct cuda_mparticles *cmprts);
void cuda_mparticles_find_block_indices_ids(struct cuda_mparticles *cmprts);
void cuda_mparticles_reorder_and_offsets(struct cuda_mparticles *cmprts);
void cuda_mparticles_check_ordered(struct cuda_mparticles *cmprts);
void cuda_mparticles_reorder(struct cuda_mparticles *cmprts);

#endif
