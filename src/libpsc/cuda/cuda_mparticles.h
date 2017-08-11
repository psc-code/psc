
#ifndef CUDA_MPARTICLES_H
#define CUDA_MPARTICLES_H

#include "cuda_iface.h"

#ifdef __cplusplus
#define EXTERN_C extern "C"
#else
#define EXTERN_C
#endif

// ----------------------------------------------------------------------
// float_3

typedef float float_3[3];

// ----------------------------------------------------------------------
// cuda_mparticles

struct cuda_mparticles {
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
  float_3 *xb_by_patch;            // lower left corner for each patch

  bool need_reorder;              // particles haven't yet been put into their sorted order
};

void cuda_mparticles_swap_alt(struct cuda_mparticles *cmprts);
void cuda_mparticles_find_block_indices_ids(struct cuda_mparticles *cmprts);
void cuda_mparticles_reorder_and_offsets(struct cuda_mparticles *cmprts);
void cuda_mparticles_check_ordered(struct cuda_mparticles *cmprts);
void cuda_mparticles_reorder(struct cuda_mparticles *cmprts);

#endif
