
#ifndef CUDA_MPARTICLES_INDEXER_H
#define CUDA_MPARTICLES_INDEXER_H

#include "grid.hxx"
#include "psc_particles_cuda.h"
#include "range.hxx"

#define CUDA_BND_S_NEW (9)
#define CUDA_BND_S_OOB (10)
#define CUDA_BND_STRIDE (10)

// ======================================================================
// cuda_mparticles_indexer

template<typename BS>
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

protected:
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

// ======================================================================
// DParticleIndexer

template<typename BS>
struct DParticleIndexer
{
  using real_t = float;

  DParticleIndexer() = default; // FIXME, delete

  DParticleIndexer(const cuda_mparticles_indexer<BS>& cpi)
  {
    for (int d = 0; d < 3; d++) {
      b_mx_[d]  = cpi.b_mx()[d];
      b_dxi_[d] = cpi.b_dxi()[d];
      dxi_[d]   = cpi.dxi()[d];
    }
  }

  __device__ int blockIndex(float4 xi4, int p) const
  {
    uint block_pos_y = __float2int_rd(xi4.y * b_dxi_[1]);
    uint block_pos_z = __float2int_rd(xi4.z * b_dxi_[2]);
    
    //assert(block_pos_y < b_mx_[1] && block_pos_z < b_mx_[2]); FIXME, assert doesn't work (on macbook)
    return (p * b_mx_[2] + block_pos_z) * b_mx_[1] + block_pos_y;
  }

  __device__ int blockShift(float xi[3], int p, int bid) const
  {
    uint block_pos_y = __float2int_rd(xi[1] * b_dxi_[1]);
    uint block_pos_z = __float2int_rd(xi[2] * b_dxi_[2]);
    
    if (block_pos_y >= b_mx_[1] || block_pos_z >= b_mx_[2]) {
      return CUDA_BND_S_OOB;
    } else {
      int bidx = (p * b_mx_[2] + block_pos_z) * b_mx_[1] + block_pos_y;
      int b_diff = bid - bidx + b_mx_[1] + 1;
      int d1 = b_diff % b_mx_[1];
      int d2 = b_diff / b_mx_[1];
      return d2 * 3 + d1;
    }
  }
  
  __device__ uint block_pos_to_block_idx(int block_pos[3]) const
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

  __device__ int find_bid() const
  {
    return blockIdx.y * b_mx_[1] + blockIdx.x;
  }

  __device__ int find_bid_q(int p, int *block_pos) const
  {
    // FIXME won't work if b_mx_[1,2] not even (?)
    return block_pos_to_block_idx(block_pos) + p * b_mx_[1] * b_mx_[2];
  }

  __device__ int find_block_pos_patch(int *block_pos, int *ci0) const
  {
    block_pos[1] = blockIdx.x;
    block_pos[2] = blockIdx.y % b_mx_[2];
    
    ci0[0] = 0;
    ci0[1] = block_pos[1] * BS::y::value;
    ci0[2] = block_pos[2] * BS::z::value;
    
    return blockIdx.y / b_mx_[2];
  }

  __device__ int find_block_pos_patch(int *block_pos) const
  {
    block_pos[1] = blockIdx.x;
    block_pos[2] = blockIdx.y % b_mx_[2];
    
    return blockIdx.y / b_mx_[2];
  }

  __device__ int find_block_pos_patch_q(int *block_pos, int *ci0, int block_start) const
  {
    int grid_dim_y = (b_mx_[2] + 1) / 2;
    block_pos[1] = blockIdx.x * 2;
    block_pos[2] = (blockIdx.y % grid_dim_y) * 2;
    block_pos[1] += block_start & 1;
    block_pos[2] += block_start >> 1;
    if (block_pos[1] >= b_mx_[1] || block_pos[2] >= b_mx_[2])
      return -1;
    
    ci0[0] = 0;
    ci0[1] = block_pos[1] * BS::y::value;
    ci0[2] = block_pos[2] * BS::z::value;
    
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

  // ----------------------------------------------------------------------
  // find_idx_off_pos_1st
  
  __device__ void find_idx_off_pos_1st(const float xi[3], int j[3], float h[3], float pos[3], float shift)
  {
    for (int d = 0; d < 3; d++) {
      pos[d] = xi[d] * dxi_[d] + shift;
      j[d] = __float2int_rd(pos[d]);
      h[d] = pos[d] - j[d];
    }
  }

  __device__ void scalePos(real_t xs[3], real_t xi[3])
  {
    xs[0] = 0.;
    xs[1] = xi[1] * dxi_[1];
    xs[2] = xi[2] * dxi_[2];
  }

  __device__ real_t dxi(int d) const { return dxi_[d]; }

  __device__ uint b_mx(int d) const { return b_mx_[d]; }

private:
  uint b_mx_[3];
  real_t b_dxi_[3];
  real_t dxi_[3];
};

struct BlockBase
{
  int bid;
  int p;
  int ci0[3];
};

template<typename BS>
struct BlockSimple : BlockBase
{
  __device__
  bool init(const DParticleIndexer<BS>& dpi, int block_start = 0)
  {
    int block_pos[3];
    p = dpi.find_block_pos_patch(block_pos, ci0);
    if (p < 0) {
      return false;
    }
    bid = dpi.find_bid();
    return true;
  }
};

template<typename BS>
struct BlockQ : BlockBase
{
  __device__
  int init(const DParticleIndexer<BS>& dpi, int block_start)
  {
    int block_pos[3];
    p = dpi.find_block_pos_patch_q(block_pos, ci0, block_start);
    if (p < 0) {
      return false;
    }
    bid = dpi.find_bid_q(p, block_pos);
    return true;
  }
};

__device__
inline RangeStrided<uint> in_block_loop(uint block_begin, uint block_end)
{
  return range((block_begin & ~31) + threadIdx.x, block_end, blockDim.x);
}

#endif


