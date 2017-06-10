
#ifndef CUDA_MPARTICLES_H
#define CUDA_MPARTICLES_H

#define cudaCheck(ierr) do {						\
    if (ierr != cudaSuccess)						\
      fprintf(stderr, "IERR = %d (%s)\n", ierr, cudaGetErrorName(ierr)); \
    assert(ierr == cudaSuccess);					\
  } while(0)

#define cuda_sync_if_enabled() do {					\
    if (1) {								\
      cudaError_t ierr = cudaThreadSynchronize(); cudaCheck(ierr);	\
    }									\
  } while(0)

// ----------------------------------------------------------------------
// cuda_domain_info

struct cuda_domain_info {
  int nr_patches;
  int mx[3]; // number of cells per patch
  int bs[3]; // size of each block / super-cell
  double dx[3]; // size of a single cell
};

// ----------------------------------------------------------------------
// cuda_mparticles

struct cuda_mparticles {
  unsigned int nr_prts;     // total # of particles in all patches
  unsigned int nr_alloced;  // arrays are alloced for this # of particles
  unsigned int nr_patches;
  unsigned int nr_blocks_per_patch;
  unsigned int nr_blocks;
  int mx[3];      // number of cells per direction in each patch
  int b_mx[3];    // number of blocks per direction in each patch
  float dx[3];    // cell size (in actual length units)
  float b_dxi[3]; // inverse of block size (in actual length units)

  // per particle
  float4 *d_xi4, *d_pxi4;
  float4 *d_alt_xi4, *d_alt_pxi4;
  unsigned int *d_bidx;
  unsigned int *d_id;

  // per block
  unsigned int *d_off;

  // temporary: particles in each patch, used before initial sort, ie., particles
  // are arranged by patch order only, but not sorted by block yet, so d_off isn't
  // set yet, either
  unsigned int *d_nr_prts_by_patch;

};

void cuda_mparticles_set_domain_info(struct cuda_mparticles *cuda_mprts,
				     const struct cuda_domain_info *info);
void cuda_mparticles_alloc(struct cuda_mparticles *cuda_mprts, unsigned int *nr_prts_by_patch);
void cuda_mparticles_free(struct cuda_mparticles *cuda_mprts);
void cuda_mparticles_dump(struct cuda_mparticles *cuda_mprts);
void cuda_mparticles_find_block_indices_ids_total(struct cuda_mparticles *cuda_mprts,
						  unsigned int *nr_parts_by_patch);
void cuda_mparticles_reorder_and_offsets(struct cuda_mparticles *cuda_mprts);

#endif

