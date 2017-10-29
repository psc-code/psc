
#ifndef CUDA_COMMON_H
#define CUDA_COMMON_H

#define DIM_Z 1
#define DIM_YZ 2

template<int NBLOCKS_X, int NBLOCKS_Y, int NBLOCKS_Z>
static __device__ __host__ __forceinline__ void
block_idx_to_block_pos(int block_idx, int block_pos[3])
{
  block_pos[0] = 0;
  block_pos[1] = block_idx % NBLOCKS_Y;
  block_pos[2] = block_idx / NBLOCKS_Y;
}

#endif
