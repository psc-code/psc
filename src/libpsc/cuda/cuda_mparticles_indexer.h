
#ifndef CUDA_MPARTICLES_INDEXER_H
#define CUDA_MPARTICLES_INDEXER_H

#include "grid.hxx"

// ======================================================================
// cuda_mparticles_indexer

struct cuda_mparticles_indexer
{
  using particle_t = particle_cuda_t;
  using real_t = particle_t::real_t;
  using Real3 = Vec3<real_t>;

  cuda_mparticles_indexer() = default; // FIXME, should go away

  cuda_mparticles_indexer(const Grid_t& grid)
  {
    for (int d = 0; d < 3; d++) {
      assert(grid.ldims[d] % grid.bs[d] == 0);
      b_mx_[d] = grid.ldims[d] / grid.bs[d];
      b_dxi_[d] = 1.f / (grid.bs[d] * grid.dx[d]);
    }
  }
  
  Int3 b_mx_;          // number of blocks per direction in each patch
  Real3 b_dxi_;        // inverse of block size (in actual length units)
};


#endif


