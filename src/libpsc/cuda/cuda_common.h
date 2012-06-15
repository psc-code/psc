
#ifndef CUDA_COMMON_H
#define CUDA_COMMON_H

#define DIM_Z 1
#define DIM_YZ 2

// OPT, we're not using the log / shift versions anymore, they might be faster
// (but the compiler should figure it out on its own)

template<int NBLOCKS_X, int NBLOCKS_Y, int NBLOCKS_Z>
static __device__ __host__ __forceinline__ unsigned int
block_pos_to_block_idx(int block_pos[3])
{
  return (block_pos[2] *  NBLOCKS_Y) | block_pos[1];
}

static __device__ inline unsigned int
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

template<int NBLOCKS_X, int NBLOCKS_Y, int NBLOCKS_Z>
static __device__ __host__ __forceinline__ void
block_idx_to_block_pos(int block_idx, int block_pos[3])
{
  block_pos[0] = 0;
  block_pos[1] = block_idx % NBLOCKS_Y;
  block_pos[2] = block_idx / NBLOCKS_Y;
}

__device__ static void
find_idx_off_1st(const real xi[3], int j[3], real h[3], real shift, real dxi[3])
{
  int d;
  for (d = 0; d < 3; d++) {
    real pos = xi[d] * dxi[d] + shift;
    j[d] = __float2int_rd(pos);
    h[d] = pos - j[d];
  }
}

#endif
