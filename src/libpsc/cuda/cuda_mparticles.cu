
#include "cuda_mparticles.h"

#include <cstdio>
#include <cassert>

#define cudaCheck(ierr) do {						\
    if (ierr != cudaSuccess)						\
      fprintf(stderr, "IERR = %d (%s)\n", ierr, cudaGetErrorName(ierr)); \
    assert(ierr == cudaSuccess);					\
  } while(0)

#define cuda_sync_if_enabled() do {					\
    if (1) {								\
      cudaError_t ierr = cudaThreadSynchronize(); cudaCheck(ierr);	\
    }									\
  } while(0)

// ----------------------------------------------------------------------
// cuda_mparticles_create

struct cuda_mparticles *
cuda_mparticles_create()
{
  struct cuda_mparticles *cmprts = 
    (struct cuda_mparticles *) calloc(1, sizeof(*cmprts));

  return cmprts;
}

// ----------------------------------------------------------------------
// cuda_mparticles_destroy

void
cuda_mparticles_destroy(struct cuda_mparticles *cmprts)
{
  free(cmprts);
}

// ----------------------------------------------------------------------
// cuda_mparticles_set_domain_info

void
cuda_mparticles_set_domain_info(struct cuda_mparticles *cmprts,
				const struct cuda_domain_info *info)
{
  cmprts->n_patches = info->n_patches;
  for (int d = 0; d < 3; d++) {
    cmprts->mx[d] = info->mx[d];
    cmprts->b_mx[d] = info->mx[d] / info->bs[d];
    cmprts->dx[d] = info->dx[d];
    cmprts->b_dxi[d] = 1.f / (info->bs[d] * info->dx[d]);
  }
  cmprts->n_blocks_per_patch = cmprts->b_mx[0] * cmprts->b_mx[1] * cmprts->b_mx[2];
  cmprts->n_blocks = cmprts->n_patches * cmprts->n_blocks_per_patch;
}

// ----------------------------------------------------------------------
// cuda_mparticles_alloc

void
cuda_mparticles_alloc(struct cuda_mparticles *cmprts, unsigned int *n_prts_by_patch)
{
  cudaError_t ierr;

  cmprts->n_prts = 0;
  for (int p = 0; p < cmprts->n_patches; p++) {
    cmprts->n_prts += n_prts_by_patch[p];
  }

  cmprts->n_alloced = cmprts->n_prts * 1.4;
  unsigned int n_alloced = cmprts->n_alloced;

  ierr = cudaMalloc((void **) &cmprts->d_xi4, n_alloced * sizeof(float4)); cudaCheck(ierr);
  ierr = cudaMalloc((void **) &cmprts->d_pxi4, n_alloced * sizeof(float4)); cudaCheck(ierr);
  ierr = cudaMalloc((void **) &cmprts->d_alt_xi4, n_alloced * sizeof(float4)); cudaCheck(ierr);
  ierr = cudaMalloc((void **) &cmprts->d_alt_pxi4, n_alloced * sizeof(float4)); cudaCheck(ierr);
  ierr = cudaMalloc((void **) &cmprts->d_bidx, n_alloced * sizeof(unsigned int)); cudaCheck(ierr);
  ierr = cudaMalloc((void **) &cmprts->d_id, n_alloced * sizeof(unsigned int)); cudaCheck(ierr);

  ierr = cudaMalloc(&cmprts->d_n_prts_by_patch, cmprts->n_patches * sizeof(unsigned int)); cudaCheck(ierr);

  ierr = cudaMalloc(&cmprts->d_off, (cmprts->n_blocks + 1) * sizeof(*cmprts->d_off)); cudaCheck(ierr);
}

// ----------------------------------------------------------------------
// cuda_mparticles_dealloc

void
cuda_mparticles_dealloc(struct cuda_mparticles *cmprts)
{
  cudaError_t ierr;

  ierr = cudaFree(cmprts->d_xi4); cudaCheck(ierr);
  ierr = cudaFree(cmprts->d_pxi4); cudaCheck(ierr);
  ierr = cudaFree(cmprts->d_alt_xi4); cudaCheck(ierr);
  ierr = cudaFree(cmprts->d_alt_pxi4); cudaCheck(ierr);
  ierr = cudaFree(cmprts->d_bidx); cudaCheck(ierr);
  ierr = cudaFree(cmprts->d_id); cudaCheck(ierr);

  ierr = cudaFree(cmprts->d_n_prts_by_patch); cudaCheck(ierr);

  ierr = cudaFree(cmprts->d_off); cudaCheck(ierr);
}

// ----------------------------------------------------------------------
// cuda_mparticles_swap_alt

void
cuda_mparticles_swap_alt(struct cuda_mparticles *cmprts)
{
  float4 *tmp_xi4 = cmprts->d_alt_xi4;
  float4 *tmp_pxi4 = cmprts->d_alt_pxi4;
  cmprts->d_alt_xi4 = cmprts->d_xi4;
  cmprts->d_alt_pxi4 = cmprts->d_pxi4;
  cmprts->d_xi4 = tmp_xi4;
  cmprts->d_pxi4 = tmp_pxi4;
}

// ----------------------------------------------------------------------
// cuda_params2

struct cuda_params2 {
  unsigned int b_mx[3];
  float b_dxi[3];
};

static void
cuda_params2_set(struct cuda_params2 *prm, const struct cuda_mparticles *cuda_mprts)
{
  for (int d = 0; d < 3; d++) {
    prm->b_mx[d]  = cuda_mprts->b_mx[d];
    prm->b_dxi[d] = cuda_mprts->b_dxi[d];
  }
}

static void
cuda_params2_free(struct cuda_params2 *prm)
{
}

#define THREADS_PER_BLOCK 256

// ----------------------------------------------------------------------
// cuda_mprts_find_block_indices_ids

__global__ static void
mprts_find_block_indices_ids(struct cuda_params2 prm, float4 *d_xi4,
			     int *d_n_prts_by_patch, unsigned int *d_bidx,
			     unsigned int *d_ids, int nr_patches)
{
  int n = threadIdx.x + THREADS_PER_BLOCK * blockIdx.x;
  int nr_blocks = prm.b_mx[1] * prm.b_mx[2];

  unsigned int off = 0;
  for (int p = 0; p < nr_patches; p++) {
    if (n < d_n_prts_by_patch[p]) {
      float4 xi4 = d_xi4[n + off];
      unsigned int block_pos_y = __float2int_rd(xi4.y * prm.b_dxi[1]);
      unsigned int block_pos_z = __float2int_rd(xi4.z * prm.b_dxi[2]);
      
      int block_idx;
      if (block_pos_y >= prm.b_mx[1] || block_pos_z >= prm.b_mx[2]) {
	block_idx = -1; // not supposed to happen here!
      } else {
	block_idx = block_pos_z * prm.b_mx[1] + block_pos_y + p * nr_blocks;
      }
      d_bidx[n + off] = block_idx;
      d_ids[n + off] = n + off;
    }
    off += d_n_prts_by_patch[p];
  }
}

void
cuda_mparticles_find_block_indices_ids(struct cuda_mparticles *cmprts,
				       unsigned int *n_prts_by_patch)
{
  if (cmprts->n_patches == 0) {
    return;
  }

  int max_n_prts = 0;
  for (int p = 0; p < cmprts->n_patches; p++) {
    if (n_prts_by_patch[p] > max_n_prts) {
      max_n_prts = n_prts_by_patch[p];
    }
  }

  cudaError_t ierr;
  ierr = cudaMemcpy(cmprts->d_n_prts_by_patch, n_prts_by_patch,
		    cmprts->n_patches * sizeof(unsigned int),
		    cudaMemcpyHostToDevice); cudaCheck(ierr);

  struct cuda_params2 prm;
  cuda_params2_set(&prm, cmprts);
    
  dim3 dimGrid((max_n_prts + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK);
  dim3 dimBlock(THREADS_PER_BLOCK);

  mprts_find_block_indices_ids<<<dimGrid, dimBlock>>>(prm,
						      cmprts->d_xi4, 
						      cmprts->d_n_prts_by_patch,
						      cmprts->d_bidx,
						      cmprts->d_id,
						      cmprts->n_patches);
  cuda_sync_if_enabled();
  cuda_params2_free(&prm);
}


