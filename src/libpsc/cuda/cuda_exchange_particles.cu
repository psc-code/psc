
#undef _GLIBCXX_USE_INT128

#include "cuda_mparticles.h"
#include "cuda_mparticles_const.h"

#include "psc_bnd_cuda.h"

#define THREADS_PER_BLOCK 256

// ----------------------------------------------------------------------
// cuda_mprts_find_block_indices_2_total
//
// like cuda_find_block_indices, but handles out-of-bound
// particles

__global__ static void
mprts_find_block_indices_2_total(float4 *d_xi4, unsigned int *d_off,
				 unsigned int *d_bidx, int nr_patches)
{
  int tid = threadIdx.x;

  int block_pos[3];
  block_pos[1] = blockIdx.x;
  block_pos[2] = blockIdx.y % d_cmprts_const.b_mx[2];
  int bid = block_pos_to_block_idx(block_pos, d_cmprts_const.b_mx);
  int p = blockIdx.y / d_cmprts_const.b_mx[2];

  int nr_blocks = d_cmprts_const.b_mx[1] * d_cmprts_const.b_mx[2];

  // FIXME/OPT, could be done better like reorder_send_buf
  int block_begin = d_off[bid + p * nr_blocks];
  int block_end   = d_off[bid + p * nr_blocks + 1];

  for (int n = block_begin + tid; n < block_end; n += THREADS_PER_BLOCK) {
    float4 xi4 = d_xi4[n];
    unsigned int block_pos_y = __float2int_rd(xi4.y * d_cmprts_const.b_dxi[1]);
    unsigned int block_pos_z = __float2int_rd(xi4.z * d_cmprts_const.b_dxi[2]);

    int block_idx;
    if (block_pos_y >= d_cmprts_const.b_mx[1] || block_pos_z >= d_cmprts_const.b_mx[2]) {
      block_idx = nr_blocks * nr_patches;
    } else {
      block_idx = block_pos_z * d_cmprts_const.b_mx[1] + block_pos_y + p * nr_blocks;
    }
    d_bidx[n] = block_idx;
  }
}

EXTERN_C void
cuda_mprts_find_block_indices_2_total(struct psc_mparticles *mprts)
{
  struct cuda_mparticles *cmprts = psc_mparticles_cuda(mprts)->cmprts;

  if (mprts->nr_patches == 0) {
    return;
  }

  cuda_mparticles_const_set(cmprts);
    
  dim3 dimBlock(THREADS_PER_BLOCK);
  dim3 dimGrid(cmprts->b_mx[1], cmprts->b_mx[2] * cmprts->n_patches);
  
  mprts_find_block_indices_2_total<<<dimGrid, dimBlock>>>(cmprts->d_xi4, cmprts->d_off,
							  cmprts->d_bidx, mprts->nr_patches);
  cuda_sync_if_enabled();
}

// ----------------------------------------------------------------------
// cuda_mprts_find_block_keys

__global__ static void
mprts_find_block_keys(float4 *d_xi4, unsigned int *d_off,
		      unsigned int *d_bidx, int nr_total_blocks)
{
  int tid = threadIdx.x;
  int bid = blockIdx.x;

  int nr_blocks = d_cmprts_const.b_mx[1] * d_cmprts_const.b_mx[2];
  int p = bid / nr_blocks;

  int block_begin = d_off[bid];
  int block_end   = d_off[bid + 1];

  for (int n = block_begin + tid; n < block_end; n += THREADS_PER_BLOCK) {
    float4 xi4 = d_xi4[n];
    unsigned int block_pos_y = __float2int_rd(xi4.y * d_cmprts_const.b_dxi[1]);
    unsigned int block_pos_z = __float2int_rd(xi4.z * d_cmprts_const.b_dxi[2]);

    int block_idx;
    if (block_pos_y >= d_cmprts_const.b_mx[1] || block_pos_z >= d_cmprts_const.b_mx[2]) {
      block_idx = CUDA_BND_S_OOB;
    } else {
      int bidx = block_pos_z * d_cmprts_const.b_mx[1] + block_pos_y + p * nr_blocks;
      int b_diff = bid - bidx + d_cmprts_const.b_mx[1] + 1;
      int d1 = b_diff % d_cmprts_const.b_mx[1];
      int d2 = b_diff / d_cmprts_const.b_mx[1];
      block_idx = d2 * 3 + d1;
    }
    d_bidx[n] = block_idx;
  }
}

EXTERN_C void
cuda_mprts_find_block_keys(struct psc_mparticles *mprts)
{
  struct cuda_mparticles *cmprts = psc_mparticles_cuda(mprts)->cmprts;

  if (mprts->nr_patches == 0) {
    return;
  }

  dim3 dimBlock(THREADS_PER_BLOCK);
  dim3 dimGrid(cmprts->n_blocks);
  
  mprts_find_block_keys<<<dimGrid, dimBlock>>>(cmprts->d_xi4, cmprts->d_off,
					       cmprts->d_bidx, cmprts->n_blocks);
}

