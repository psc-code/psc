
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
    n_patches = grid.n_patches();
    n_blocks_per_patch = pi_.b_mx_[0] * pi_.b_mx_[1] * pi_.b_mx_[2];
    n_blocks = n_patches * n_blocks_per_patch;
  }

  Int3 blockPosition(const real_t xi[3]) const { return pi_.blockPosition(xi); }

  int blockIndex(const float4& xi4, int p) const
  {
    Int3 bpos = blockPosition(&xi4.x);
    return blockIndex(bpos, p);
  }

  void checkInPatchMod(real_t xi[3]) const { pi_.checkInPatchMod(xi); }
  const Real3& b_dxi() const { return pi_.b_dxi_; }
  const Int3& b_mx() const { return pi_.b_mx_; }

private:
  int blockIndex(Int3 bpos, int p) const
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
private:
  ParticleIndexer<real_t> pi_;
};


#endif


