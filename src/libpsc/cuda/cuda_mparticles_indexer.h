
#ifndef CUDA_MPARTICLES_INDEXER_H
#define CUDA_MPARTICLES_INDEXER_H

#include "grid.hxx"
#include "psc_particles_cuda.h"

// ======================================================================
// DParticleIndexer

struct cuda_mparticles_indexer;

struct DParticleIndexer
{
  using real_t = float;

  DParticleIndexer(const cuda_mparticles_indexer& cpi);

  __device__ int blockIndex(float4 xi4, int p) const
  {
    uint block_pos_y = __float2int_rd(xi4.y * b_dxi[1]);
    uint block_pos_z = __float2int_rd(xi4.z * b_dxi[2]);
    
    //assert(block_pos_y < b_mx[1] && block_pos_z < b_mx[2]); FIXME, assert doesn't work (on macbook)
    return (p * b_mx[2] + block_pos_z) * b_mx[1] + block_pos_y;
  }
  
  uint b_mx[3];
  real_t b_dxi[3];
};

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

  operator DParticleIndexer()
  {
    return DParticleIndexer(*this);
  }
  
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


inline DParticleIndexer::DParticleIndexer(const cuda_mparticles_indexer& cpi)
{
  for (int d = 0; d < 3; d++) {
    b_mx[d]  = cpi.b_mx()[d];
    b_dxi[d] = cpi.b_dxi()[d];
  }
}

#endif


