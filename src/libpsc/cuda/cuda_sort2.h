
#ifndef CUDA_SORT2_H
#define CUDA_SORT2_H

#include "psc_cuda.h"
#include "cuda_common.h"

#include <mrc_profile.h>

#include <thrust/device_vector.h>
#include <b40c/kernel_utils.h>

#undef FAST_COMPILE
#ifdef FAST_COMPILE
#define _NBLOCKS_X 1
#define _NBLOCKS_Y 16
#define _NBLOCKS_Z 16
#endif

template<int NBLOCKS_X, int NBLOCKS_Y, int NBLOCKS_Z>
static __device__ __host__ __forceinline__ unsigned int
block_idx_to_dir1(unsigned int b_base, unsigned int b2)
{
  int bp_base[3], bp[3];
  block_idx_to_block_pos<NBLOCKS_X, NBLOCKS_Y, NBLOCKS_Z>(b_base, bp_base);
  block_idx_to_block_pos<NBLOCKS_X, NBLOCKS_Y, NBLOCKS_Z>(b2, bp);
  unsigned int dir1 = bp[1] - bp_base[1] + 1;
  unsigned int dir2 = bp[2] - bp_base[2] + 1;
  dir1 &= 3;
  dir2 &= 3;
  return (dir2 << 2) | dir1;
}

template<int NBLOCKS_X, int NBLOCKS_Y, int NBLOCKS_Z>
static __device__ __host__ __forceinline__ unsigned int
dir1_to_block_idx(unsigned int b_base, unsigned int d)
{
  int bp_base[3], bp[3];
  block_idx_to_block_pos<NBLOCKS_X, NBLOCKS_Y, NBLOCKS_Z>(b_base, bp_base);
  int dir1 = (d & 3) - 1;
  int dir2 = (d >> 2) - 1;
  bp[1] = (bp_base[1] + dir1) & (NBLOCKS_Y - 1);
  bp[2] = (bp_base[2] + dir2) & (NBLOCKS_Z - 1);
  return block_pos_to_block_idx<NBLOCKS_X, NBLOCKS_Y, NBLOCKS_Z>(bp);
}

#endif
