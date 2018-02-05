
#ifndef CUDA_MPARTICLES_INDEXER_H
#define CUDA_MPARTICLES_INDEXER_H

#include "grid.hxx"
#include "psc_particles_cuda.h"

// ======================================================================
// cuda_mparticles_indexer

struct cuda_mparticles_indexer
{
  using particle_t = particle_cuda_t;
  using real_t = particle_t::real_t;
  using Real3 = Vec3<real_t>;

  cuda_mparticles_indexer(const Grid_t& grid)
    : pi_(grid)
  {
    n_patches = grid.patches.size();
    for (int d = 0; d < 3; d++) {
      assert(grid.ldims[d] % grid.bs[d] == 0);
      b_mx_[d] = grid.ldims[d] / grid.bs[d];
    }

    n_blocks_per_patch = b_mx_[0] * b_mx_[1] * b_mx_[2];
    n_blocks = n_patches * n_blocks_per_patch;
  }

  int get_bidx(Int3 bpos, int p)
  {
    if (uint(bpos[0]) >= b_mx_[0] ||
	uint(bpos[1]) >= b_mx_[1] ||
	uint(bpos[2]) >= b_mx_[2]) {
      return -1;
    } else {
      return ((p * b_mx_[2] + bpos[2]) * b_mx_[1] + bpos[1]) * b_mx_[0] + bpos[0];
    }
  }

  int get_bidx(Int3 bpos)
  {
    return get_bidx(bpos, 0);
  }

  uint n_patches;                // number of patches
  uint n_blocks_per_patch;       // number of blocks per patch
  uint n_blocks;                 // number of blocks in all patches in mprts
  Int3 b_mx_;                    // number of blocks per direction in each patch
  ParticleIndexer<real_t> pi_;
};


#endif


