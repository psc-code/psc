
#ifndef CUDA_MPARTICLES_CONST_H
#define CUDA_MPARTICLES_CONST_H

#include <cuda_bits.h>
#include <cstdio>

// ----------------------------------------------------------------------
// cuda_mparticles_const
//
// cuda_mparticles parameters in CUDA constant memory
// (to avoid having to pass them in)

#define MAX_KINDS (4)

struct cuda_mparticles_const {
  float dxi[3];
  float b_dxi[3];
  int b_mx[3];
  float dt;
  float dqs;
  float dq[MAX_KINDS];
  float q_inv[MAX_KINDS];
  float fnqs;
  float fnqxs, fnqys, fnqzs;
};

__constant__ __device__ struct cuda_mparticles_const d_cmprts_const;

static void
cuda_mparticles_const_set(struct cuda_mparticles *cmprts)
{
  const Grid_t& grid = cmprts->grid_;
  
  struct cuda_mparticles_const c;
  for (int d = 0; d < 3; d++) {
    c.dxi[d] = 1.f / grid.dx[d];
    c.b_mx[d] = cmprts->b_mx_[d];
    c.b_dxi[d] = cmprts->b_dxi_[d];
  }

  c.dt = grid.dt;
  c.fnqs  = grid.fnqs;
  c.fnqxs = grid.dx[0] * c.fnqs / grid.dt;
  c.fnqys = grid.dx[1] * c.fnqs / grid.dt;
  c.fnqzs = grid.dx[2] * c.fnqs / grid.dt;
  c.dqs = .5f * grid.eta * grid.dt;

  int n_kinds = grid.kinds.size();
  assert(n_kinds <= MAX_KINDS);
  for (int k = 0; k < n_kinds; k++) {
    c.dq[k] = c.dqs * grid.kinds[k].q / grid.kinds[k].m;
    c.q_inv[k] = 1.f / grid.kinds[k].q;
  }

  cudaError_t ierr = cudaMemcpyToSymbol(d_cmprts_const, &c, sizeof(c)); cudaCheck(ierr);
}

// ----------------------------------------------------------------------
// struct d_particle

struct d_particle {
  float xi[3];
  float kind_as_float;
  float pxi[3];
  float qni_wni;
};

#define LOAD_PARTICLE_POS(pp, d_xi4, n) do {				\
    float4 _xi4 = d_xi4[n];						\
    (pp).xi[0]         = _xi4.x;					\
    (pp).xi[1]         = _xi4.y;					\
    (pp).xi[2]         = _xi4.z;					\
    (pp).kind_as_float = _xi4.w;					\
} while (0)

#define LOAD_PARTICLE_MOM(pp, d_pxi4, n) do {				\
    float4 _pxi4 = d_pxi4[n];						\
    (pp).pxi[0]        = _pxi4.x;					\
    (pp).pxi[1]        = _pxi4.y;					\
    (pp).pxi[2]        = _pxi4.z;					\
    (pp).qni_wni       = _pxi4.w;					\
} while (0)

#define STORE_PARTICLE_POS(pp, d_xi4, n) do {				\
    float4 xi4 = { (pp).xi[0], (pp).xi[1], (pp).xi[2], (pp).kind_as_float }; \
    d_xi4[n] = xi4;							\
} while (0)

#define STORE_PARTICLE_MOM(pp, d_pxi4, n) do {				\
    float4 pxi4 = { (pp).pxi[0], (pp).pxi[1], (pp).pxi[2], (pp).qni_wni }; \
    d_pxi4[n] = pxi4;							\
} while (0)

// ----------------------------------------------------------------------
// find_idx_off_1st

__device__ static void
find_idx_off_1st(const float xi[3], int j[3], float h[3], float shift)
{
  for (int d = 0; d < 3; d++) {
    float pos = xi[d] * d_cmprts_const.dxi[d] + shift;
    j[d] = __float2int_rd(pos);
    h[d] = pos - j[d];
  }
}

// ----------------------------------------------------------------------
// find_idx_off_pos_1st

__device__ static void
find_idx_off_pos_1st(const float xi[3], int j[3], float h[3], float pos[3], float shift)
{
  for (int d = 0; d < 3; d++) {
    pos[d] = xi[d] * d_cmprts_const.dxi[d] + shift;
    j[d] = __float2int_rd(pos[d]);
    h[d] = pos[d] - j[d];
  }
}

// ----------------------------------------------------------------------
// block_pos_to_block_idx

// OPT, we're not using the log / shift versions anymore, they might be faster
// (but the compiler should figure it out on its own)

template<int NBLOCKS_X, int NBLOCKS_Y, int NBLOCKS_Z>
static __device__ __host__ __forceinline__ uint
block_pos_to_block_idx(int block_pos[3])
{
  return (block_pos[2] *  NBLOCKS_Y) | block_pos[1];
}

#define NO_CHECKERBOARD
static __device__ inline uint
block_pos_to_block_idx(int block_pos[3], int b_mx[3])
{
#if 1 // DIM == DIM_YZ FIXME
#ifdef NO_CHECKERBOARD
  return block_pos[2] * b_mx[1] + block_pos[1];
#else
  int dimy = b_mx[1] >> 1;
  return
    ((block_pos[1] & 1) << 0) |
    ((block_pos[2] & 1) << 1) | 
    (((block_pos[2] >> 1) * dimy + (block_pos[1] >> 1)) << 2);
#endif
#else
#error TBD
#endif
}

// ----------------------------------------------------------------------
// find_bid
//
// FIXME, this is here to consolidate between moments / particle pusher
// but it's really an implementation detail

__device__ static int
find_bid()
{
  return blockIdx.y * d_cmprts_const.b_mx[1] + blockIdx.x;
}

__device__ static int
find_bid_q(int p, int *block_pos)
{
  // FIXME won't work if b_mx[1,2] not even (?)
  return block_pos_to_block_idx(block_pos, d_cmprts_const.b_mx) + p * d_cmprts_const.b_mx[1] * d_cmprts_const.b_mx[2];
}

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
__device__ static int
find_block_pos_patch(int *block_pos, int *ci0)
{
  block_pos[1] = blockIdx.x;
  block_pos[2] = blockIdx.y % d_cmprts_const.b_mx[2];

  ci0[0] = 0;
  ci0[1] = block_pos[1] * BLOCKSIZE_Y;
  ci0[2] = block_pos[2] * BLOCKSIZE_Z;

  return blockIdx.y / d_cmprts_const.b_mx[2];
}

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
__device__ static int
find_block_pos_patch_q(int *block_pos, int *ci0, int block_start)
{
  int grid_dim_y = (d_cmprts_const.b_mx[2] + 1) / 2;
  block_pos[1] = blockIdx.x * 2;
  block_pos[2] = (blockIdx.y % grid_dim_y) * 2;
  block_pos[1] += block_start & 1;
  block_pos[2] += block_start >> 1;
  if (block_pos[1] >= d_cmprts_const.b_mx[1] ||
      block_pos[2] >= d_cmprts_const.b_mx[2])
    return -1;

  ci0[0] = 0;
  ci0[1] = block_pos[1] * BLOCKSIZE_Y;
  ci0[2] = block_pos[2] * BLOCKSIZE_Z;

  return blockIdx.y / grid_dim_y;
}

#endif

