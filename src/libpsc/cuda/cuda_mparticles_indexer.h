
#ifndef CUDA_MPARTICLES_INDEXER_H
#define CUDA_MPARTICLES_INDEXER_H

#include "grid.hxx"
#include "psc_particles_cuda.h"
#include "range.hxx"
#include "dim.hxx"

#define CUDA_BND_S_NEW (9)
#define CUDA_BND_S_OOB (10)
#define CUDA_BND_STRIDE (10)

template<typename BS>
struct DParticleIndexer;

template<typename BS, typename DIM>
struct BlockQ;

template<typename BS, typename DIM>
struct BlockSimple;

// ======================================================================
// cuda_mparticles_indexer

template<typename BS>
struct cuda_mparticles_indexer
{
  using real_t = float;
  using Real3 = Vec3<real_t>;

  cuda_mparticles_indexer(const Grid_t& grid)
    : pi_(grid),
      b_mx_{int(grid.ldims[0] / BS::x::value),
            int(grid.ldims[1] / BS::y::value),
            int(grid.ldims[2] / BS::z::value)}
  {
    n_patches = grid.n_patches();
    n_blocks_per_patch = b_mx_[0] * b_mx_[1] * b_mx_[2];
    n_blocks = n_patches * n_blocks_per_patch;
  }

  Int3 cellPosition(const real_t xi[3]) const
  {
    return pi_.cellPosition(xi);
  }

  int validCellIndex(const float4& xi4, int p) const
  {
    Int3 cpos = cellPosition(&xi4.x);
    return validCellIndex(cpos, p);
  }

  Int3 blockPosition(const real_t xi[3]) const
  {
    Int3 pos = pi_.cellPosition(xi);
    pos[0] /= BS::x::value;
    pos[1] /= BS::y::value;
    pos[2] /= BS::z::value;
    return pos;
  }

  int blockIndex(const float4& xi4, int p) const
  {
    Int3 bpos = blockPosition(&xi4.x);
    return blockIndex(bpos, p);
  }

  void checkInPatchMod(real_t xi[3]) const { pi_.checkInPatchMod(xi); }
  const Int3& b_mx() const { return b_mx_; }

protected:
  int blockIndex(Int3 bpos, int p) const
  {
    if (uint(bpos[0]) >= b_mx_[0] ||
	uint(bpos[1]) >= b_mx_[1] ||
	uint(bpos[2]) >= b_mx_[2]) {
      return -1;
    }
    
    return ((p * b_mx_[2] + bpos[2]) * b_mx_[1] + bpos[1]) * b_mx_[0] + bpos[0];
  }

  int validCellIndex(Int3 cpos, int p) const
  {
    assert(uint(cpos[0]) < pi_.ldims_[0]);
    assert(uint(cpos[1]) < pi_.ldims_[1]);
    assert(uint(cpos[2]) < pi_.ldims_[2]);
    
    return ((p * pi_.ldims_[2] + cpos[2]) * pi_.ldims_[1] + cpos[1]) * pi_.ldims_[0] + cpos[0];
  }

public:
  uint n_cells() const { return pi_.n_cells_ * n_patches; }
  
  uint n_patches;                // number of patches
  uint n_blocks_per_patch;       // number of blocks per patch
  uint n_blocks;                 // number of blocks in all patches in mprts
private:
  ParticleIndexer<real_t> pi_;
  Int3 b_mx_;

  friend struct DParticleIndexer<BS>;
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
      ldims_[d] = cpi.pi_.ldims_[d];
      b_mx_[d]  = cpi.b_mx_[d];
      dxi_[d]   = cpi.pi_.dxi_[d];
      n_blocks_ = cpi.n_blocks;
    }
  }

  __device__ int validCellIndex(float4 xi4, int p) const
  {
    uint pos_x = __float2int_rd(xi4.x * dxi_[0]);
    uint pos_y = __float2int_rd(xi4.y * dxi_[1]);
    uint pos_z = __float2int_rd(xi4.z * dxi_[2]);
    
    //assert(pos_y < ldims_[1] && pos_z < ldims_[2]); FIXME, assert doesn't work (on macbook)
    return ((p * ldims_[2] + pos_z) * ldims_[1] + pos_y) * ldims_[0] + pos_x;
  }

  __device__ int blockIndex(float4 xi4, int p) const
  {
    int block_pos[3] = { int(__float2int_rd(xi4.x * dxi_[0]) / BS::x::value),
			 int(__float2int_rd(xi4.y * dxi_[1]) / BS::y::value),
			 int(__float2int_rd(xi4.z * dxi_[2]) / BS::z::value) };

    return validBlockIndex(block_pos, p);
  }

  __device__ int validBlockIndex(const int* block_pos, int p) const
  {
    /* assert(block_pos[0] >= 0 && block_pos[0] < b_mx_[0]); */
    /* assert(block_pos[1] >= 0 && block_pos[1] < b_mx_[1]); */
    /* assert(block_pos[2] >= 0 && block_pos[2] < b_mx_[2]); */

    return ((p * b_mx_[2] + block_pos[2]) * b_mx_[1] + block_pos[1]) * b_mx_[0] + block_pos[0];
  }

  __device__ int blockIndexFromCellPosition(const int* cpos, int p) const
  {
    if (uint(cpos[0]) >= ldims_[0] ||
	uint(cpos[1]) >= ldims_[1] ||
	uint(cpos[2]) >= ldims_[2]) {
      return n_blocks_;
    }

    int bpos[3] = { int(cpos[0] / BS::x::value),
		    int(cpos[1] / BS::y::value),
		    int(cpos[2] / BS::z::value) };
    return validBlockIndex(bpos, p);
  }
  
  __device__ int blockShift(float xi[3], int p, int bid) const
  {
    static_assert(BS::x::value == 1, "blockShift needs work for dim_xyz");
    uint block_pos_y = __float2int_rd(xi[1] * dxi_[1]) / BS::y::value;
    uint block_pos_z = __float2int_rd(xi[2] * dxi_[2]) / BS::z::value;
    
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
  
  // ======================================================================
  // cell related
  
  // ----------------------------------------------------------------------
  // find_idx_off_1st

  __device__ void find_idx_off_1st(const float xi[3], int j[3], float h[3], float shift)
  {
    for (int d = 0; d < 3; d++) {
      real_t pos = scalePos(xi[d], d) + shift;
      j[d] = __float2int_rd(pos);
      h[d] = pos - j[d];
    }
  }

  // ----------------------------------------------------------------------
  // find_idx_off_pos_1st
  
  __device__ void find_idx_off_pos_1st(const float xi[3], int j[3], float h[3], float pos[3], float shift)
  {
    for (int d = 0; d < 3; d++) {
      pos[d] = scalePos(xi[d], d) + shift;
      j[d] = __float2int_rd(pos[d]);
      h[d] = pos[d] - j[d];
    }
  }

  __device__ real_t scalePos(real_t xi, int d)
  {
    return xi * dxi_[d];
  }

  template<typename DIM>
  __device__ void scalePos(real_t xs[3], real_t xi[3])
  {
    if (DIM::InvarX::value) xs[0] = 0.; else xs[0] = scalePos(xi[0], 0);
    xs[1] = scalePos(xi[1], 1);
    xs[2] = scalePos(xi[2], 2);
  }

private:
  uint ldims_[3];
  uint b_mx_[3];
  uint n_blocks_;
  real_t dxi_[3];

  friend class BlockQ<BS, dim_yz>;
  friend class BlockQ<BS, dim_xyz>;
  friend class BlockSimple<BS, dim_yz>;
  friend class BlockSimple<BS, dim_xyz>;
};

struct BlockBase
{
  int bid;
  int p;
  int ci0[3];
};

// ======================================================================
// BlockSimple

template<typename BS, typename DIM>
struct BlockSimple : BlockBase
{
  static Range<int> block_starts() { return range(1);  }

  template<typename CudaMparticles>
  static dim3 dimGrid(CudaMparticles& cmprts)
  {
    int gx = cmprts.b_mx()[0];
    int gy = cmprts.b_mx()[1];
    int gz = cmprts.b_mx()[2] * cmprts.n_patches;
    return dim3(gx, gy, gz);
  }

  __device__
  bool init(const DParticleIndexer<BS>& dpi, int block_start = 0)
  {
    int block_pos[3];
    block_pos[0] = blockIdx.x;
    block_pos[1] = blockIdx.y;
    block_pos[2] = blockIdx.z % dpi.b_mx_[2];
    
    ci0[0] = block_pos[0] * BS::x::value;
    ci0[1] = block_pos[1] * BS::y::value;
    ci0[2] = block_pos[2] * BS::z::value;
    
    p = blockIdx.z / dpi.b_mx_[2];
    bid = (blockIdx.z * dpi.b_mx_[1] + blockIdx.y) * dpi.b_mx_[0] + blockIdx.x;
    return true;
  }
};

// ======================================================================
// BlockSimple specialized for dim_yz

template<typename BS>
struct BlockSimple<BS, dim_yz> : BlockBase
{
  static Range<int> block_starts() { return range(1);  }

  template<typename CudaMparticles>
  static dim3 dimGrid(CudaMparticles& cmprts)
  {
    int gx = cmprts.b_mx()[1];
    int gy = cmprts.b_mx()[2] * cmprts.n_patches;
    return dim3(gx, gy);
  }

  __device__
  bool init(const DParticleIndexer<BS>& dpi, int block_start = 0)
  {
    int block_pos[3];
    block_pos[1] = blockIdx.x;
    block_pos[2] = blockIdx.y % dpi.b_mx_[2];
    
    ci0[0] = 0;
    ci0[1] = block_pos[1] * BS::y::value;
    ci0[2] = block_pos[2] * BS::z::value;
    
    p = blockIdx.y / dpi.b_mx_[2];
    bid = blockIdx.y * dpi.b_mx_[2] + blockIdx.x;
    return true;
  }
};

template<typename BS, typename DIM>
struct BlockQ : BlockBase
{
  static Range<int> block_starts() { return range(8); }

  template<typename CudaMparticles>
  static dim3 dimGrid(CudaMparticles& cmprts)
  {
    int gx =  (cmprts.b_mx()[0] + 1) / 2;
    int gy =  (cmprts.b_mx()[1] + 1) / 2;
    int gz = ((cmprts.b_mx()[2] + 1) / 2) * cmprts.n_patches;
    return dim3(gx, gy, gz);
  }

  __device__
  int init(const DParticleIndexer<BS>& dpi, int block_start)
  {
    int block_pos[3];
    int grid_dim_z = (dpi.b_mx_[2] + 1) / 2;
    block_pos[0] = blockIdx.x * 2;
    block_pos[1] = blockIdx.y * 2;
    block_pos[2] = (blockIdx.z % grid_dim_z) * 2;
    block_pos[0] += block_start & 1; block_start >>= 1;
    block_pos[1] += block_start & 1; block_start >>= 1;
    block_pos[2] += block_start;
    if (block_pos[0] >= dpi.b_mx_[0] ||
	block_pos[1] >= dpi.b_mx_[1] ||
	block_pos[2] >= dpi.b_mx_[2])
      return false;
    
    ci0[0] = block_pos[0] * BS::x::value;
    ci0[1] = block_pos[1] * BS::y::value;
    ci0[2] = block_pos[2] * BS::z::value;
    
    p = blockIdx.z / grid_dim_z;

    // FIXME won't work if b_mx_[0,1,2] not even (?)
    bid = (((p
	     * dpi.b_mx_[2] + block_pos[2])
	    * dpi.b_mx_[1] + block_pos[1])
	   * dpi.b_mx_[0] + block_pos[0]);
    return true;
  }
};

// ----------------------------------------------------------------------
// BlockQ specialized for dim_yz

template<typename BS>
struct BlockQ<BS, dim_yz> : BlockBase
{
  static Range<int> block_starts() { return range(4); }

  template<typename CudaMparticles>
  static dim3 dimGrid(CudaMparticles& cmprts)
  {
    int gx =  (cmprts.b_mx()[1] + 1) / 2;
    int gy = ((cmprts.b_mx()[2] + 1) / 2) * cmprts.n_patches;
    return dim3(gx, gy);
  }

  __device__
  int init(const DParticleIndexer<BS>& dpi, int block_start)
  {
    int block_pos[3];
    int grid_dim_y = (dpi.b_mx_[2] + 1) / 2;
    block_pos[1] = blockIdx.x * 2;
    block_pos[2] = (blockIdx.y % grid_dim_y) * 2;
    block_pos[1] += block_start & 1;
    block_pos[2] += block_start >> 1;
    if (block_pos[1] >= dpi.b_mx_[1] || block_pos[2] >= dpi.b_mx_[2])
      return false;
    
    ci0[0] = 0;
    ci0[1] = block_pos[1] * BS::y::value;
    ci0[2] = block_pos[2] * BS::z::value;
    
    p = blockIdx.y / grid_dim_y;

    // FIXME won't work if b_mx_[1,2] not even (?)
    bid = (p * dpi.b_mx_[2] + block_pos[2]) * dpi.b_mx_[1] + block_pos[1];
    return true;
  }
};

__device__
inline RangeStrided<uint> in_block_loop(uint block_begin, uint block_end)
{
  return range((block_begin & ~31) + threadIdx.x, block_end, blockDim.x);
}

#endif


