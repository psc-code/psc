
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
    b_mx_ = pi_.b_mx_;
    n_patches = grid.n_patches();
    n_blocks_per_patch = b_mx_[0] * b_mx_[1] * b_mx_[2];
    n_blocks = n_patches * n_blocks_per_patch;
  }

  int blockIndex(const float4& xi4, int p)
  {
    Int3 bpos = pi_.blockPosition(&xi4.x);
    return blockIndex(bpos, p);
  }

private:
  int blockIndex(Int3 bpos, int p)
  {
    int bidx = pi_.blockIndex(bpos);
    if (bidx < 0) {
      return bidx;
    }

    return p * n_blocks_per_patch + bidx;
  }

public:
  uint n_patches;                // number of patches
  uint n_blocks_per_patch;       // number of blocks per patch
  uint n_blocks;                 // number of blocks in all patches in mprts
  Int3 b_mx_;                    // number of blocks per direction in each patch
  ParticleIndexer<real_t> pi_;
};


#endif


