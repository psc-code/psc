
#ifndef CUDA_MPARTICLES_INDEXER_H
#define CUDA_MPARTICLES_INDEXER_H

#include "grid.hxx"
#include "psc_particles_cuda.h"


#define CUDA_BND_S_NEW (9)
#define CUDA_BND_S_OOB (10)
#define CUDA_BND_STRIDE (10)

// ======================================================================
// DParticleIndexer

struct cuda_mparticles_indexer;

struct DParticleIndexer
{
  using real_t = float;

  DParticleIndexer() = default; // FIXME, delete
  DParticleIndexer(const int b_mx[3], const real_t b_dxi[3], const real_t dxi[3])
  {
    for (int d = 0; d < 3; d++) {
      b_mx_[d]  = b_mx[d];
      b_dxi_[d] = b_dxi[d];
      dxi_[d] = dxi[d];
    }
  }

  __device__ int blockIndex(float4 xi4, int p) const
  {
    uint block_pos_y = __float2int_rd(xi4.y * b_dxi_[1]);
    uint block_pos_z = __float2int_rd(xi4.z * b_dxi_[2]);
    
    //assert(block_pos_y < b_mx_[1] && block_pos_z < b_mx_[2]); FIXME, assert doesn't work (on macbook)
    return (p * b_mx_[2] + block_pos_z) * b_mx_[1] + block_pos_y;
  }

  __device__ int blockShift(float xi[3], int p_nr, int bid) const
  {
    uint block_pos_y = __float2int_rd(xi[1] * b_dxi_[1]);
    uint block_pos_z = __float2int_rd(xi[2] * b_dxi_[2]);
    
    if (block_pos_y >= b_mx_[1] || block_pos_z >= b_mx_[2]) {
      return CUDA_BND_S_OOB;
    } else {
      int bidx = (p_nr * b_mx_[2] + block_pos_z) * b_mx_[1] + block_pos_y;
      int b_diff = bid - bidx + b_mx_[1] + 1;
      int d1 = b_diff % b_mx_[1];
      int d2 = b_diff / b_mx_[1];
      return d2 * 3 + d1;
    }
  }
  
  __device__ uint block_pos_to_block_idx(int block_pos[3])
  {
#define NO_CHECKERBOARD
#ifdef NO_CHECKERBOARD
    return block_pos[2] * b_mx_[1] + block_pos[1];
#else
    int dimy = b_mx_[1] >> 1;
    return
      ((block_pos[1] & 1) << 0) |
      ((block_pos[2] & 1) << 1) | 
      (((block_pos[2] >> 1) * dimy + (block_pos[1] >> 1)) << 2);
#endif
  }

  __device__ int find_bid()
  {
    return blockIdx.y * b_mx_[1] + blockIdx.x;
  }

  __device__ int find_bid_q(int p, int *block_pos)
  {
    // FIXME won't work if b_mx_[1,2] not even (?)
    return block_pos_to_block_idx(block_pos) + p * b_mx_[1] * b_mx_[2];
  }

  template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
  __device__ int find_block_pos_patch(int *block_pos, int *ci0)
  {
    block_pos[1] = blockIdx.x;
    block_pos[2] = blockIdx.y % b_mx_[2];
    
    ci0[0] = 0;
    ci0[1] = block_pos[1] * BLOCKSIZE_Y;
    ci0[2] = block_pos[2] * BLOCKSIZE_Z;
    
    return blockIdx.y / b_mx_[2];
  }

  template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
  __device__ int find_block_pos_patch(int *block_pos)
  {
    block_pos[1] = blockIdx.x;
    block_pos[2] = blockIdx.y % b_mx_[2];
    
    return blockIdx.y / b_mx_[2];
  }

  template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
  __device__ int find_block_pos_patch_q(int *block_pos, int *ci0, int block_start)
  {
    int grid_dim_y = (b_mx_[2] + 1) / 2;
    block_pos[1] = blockIdx.x * 2;
    block_pos[2] = (blockIdx.y % grid_dim_y) * 2;
    block_pos[1] += block_start & 1;
    block_pos[2] += block_start >> 1;
    if (block_pos[1] >= b_mx_[1] || block_pos[2] >= b_mx_[2])
      return -1;
    
    ci0[0] = 0;
    ci0[1] = block_pos[1] * BLOCKSIZE_Y;
    ci0[2] = block_pos[2] * BLOCKSIZE_Z;
    
    return blockIdx.y / grid_dim_y;
  }

  // ======================================================================
  // cell related
  
  // ----------------------------------------------------------------------
  // find_idx_off_1st

  __device__ void find_idx_off_1st(const float xi[3], int j[3], float h[3], float shift)
  {
    for (int d = 0; d < 3; d++) {
      real_t pos = xi[d] * dxi_[d] + shift;
      j[d] = __float2int_rd(pos);
      h[d] = pos - j[d];
    }
  }

private:
  uint b_mx_[3];
  real_t b_dxi_[3];
  real_t dxi_[3];
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
  const Real3& dxi() const { return pi_.dxi_; }

  operator DParticleIndexer()
  {
    return DParticleIndexer(pi_.b_mx_, pi_.b_dxi_, pi_.dxi_);
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


#endif


