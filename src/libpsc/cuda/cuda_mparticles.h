
#ifndef CUDA_MPARTICLES_H
#define CUDA_MPARTICLES_H

// ----------------------------------------------------------------------
// cuda_domain_info

struct cuda_domain_info {
  int n_patches;
  int mx[3]; // number of cells per patch
  int bs[3]; // size of each block / super-cell
  double dx[3]; // size of a single cell
};

// ----------------------------------------------------------------------
// cuda_mparticles

struct cuda_mparticles {
  // per particle
  float4 *d_xi4, *d_pxi4;         // current particle data
  float4 *d_alt_xi4, *d_alt_pxi4; // storage for out-of-place reordering of particle data
  unsigned int *d_bidx;           // block index (incl patch) per particle
  unsigned int *d_id;             // particle id for sorting

  // per patch
  int *d_n_prts_by_patch;         // # of particles per batch

  // per block
  unsigned int *d_off;            // particles per block
                                  // are at indices [offsets[block] .. offsets[block+1]-1[

  unsigned int n_prts;            // total # of particles across all patches
  unsigned int n_alloced;         // size of particle-related arrays as allocated
  unsigned int n_patches;         // # of patches
  unsigned int n_blocks_per_patch;// number of blocks per patch
  unsigned int n_blocks;          // number of blocks in all patches in mprts

  int mx[3];                      // number of cells per direction in each patch
  int b_mx[3];                    // number of blocks per direction in each patch
  float dx[3];                    // cell size (in actual length units)
  float b_dxi[3];                 // inverse of block size (in actual length units)
};

struct cuda_mparticles *cuda_mparticles_create(void);
void cuda_mparticles_destroy(struct cuda_mparticles *cmprts);
void cuda_mparticles_set_domain_info(struct cuda_mparticles *cuda_mprts,
				     const struct cuda_domain_info *info);
void cuda_mparticles_alloc(struct cuda_mparticles *cmprts, unsigned int *n_prts_by_patch);
void cuda_mparticles_dealloc(struct cuda_mparticles *cmprts);
void cuda_mparticles_swap_alt(struct cuda_mparticles *cmprts);
void cuda_mparticles_find_block_indices_ids(struct cuda_mparticles *cmprts,
					    unsigned int *n_prts_by_patch);
void cuda_mparticles_reorder_and_offsets(struct cuda_mparticles *cmprts);

#endif
