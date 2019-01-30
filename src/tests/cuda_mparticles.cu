
#include "cuda_mparticles.h"

#include <cstdio>
#include <cassert>

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>

// ----------------------------------------------------------------------
// cuda_mparticles_set_domain_info

void
cuda_mparticles_set_domain_info(struct cuda_mparticles *cuda_mprts,
				const struct cuda_domain_info *info)
{
  cuda_mprts->nr_patches = info->nr_patches;
  for (int d = 0; d < 3; d++) {
    cuda_mprts->mx[d] = info->mx[d];
    cuda_mprts->b_mx[d] = info->mx[d] / info->bs[d];
    cuda_mprts->dx[d] = info->dx[d];
    cuda_mprts->b_dxi[d] = 1.f / (info->bs[d] * info->dx[d]);
  }
  cuda_mprts->nr_blocks_per_patch =
    cuda_mprts->b_mx[0] * cuda_mprts->b_mx[1] * cuda_mprts->b_mx[2];
  cuda_mprts->nr_blocks = info->nr_patches * cuda_mprts->nr_blocks_per_patch;
}

// ----------------------------------------------------------------------
// cuda_mparticles_alloc

void
cuda_mparticles_alloc(struct cuda_mparticles *cuda_mprts, unsigned int *nr_prts_by_patch)
{
  cudaError_t ierr;

  cuda_mprts->nr_prts = 0;
  for (int p = 0; p < cuda_mprts->nr_patches; p++) {
    cuda_mprts->nr_prts += nr_prts_by_patch[p];
  }

  unsigned int nr_alloced = cuda_mprts->nr_prts * 1.4;
  cuda_mprts->nr_alloced = nr_alloced;

  ierr = cudaMalloc((void **) &cuda_mprts->d_xi4, nr_alloced * sizeof(float4)); cudaCheck(ierr);
  ierr = cudaMalloc((void **) &cuda_mprts->d_pxi4, nr_alloced * sizeof(float4)); cudaCheck(ierr);
  ierr = cudaMalloc((void **) &cuda_mprts->d_alt_xi4, nr_alloced * sizeof(float4)); cudaCheck(ierr);
  ierr = cudaMalloc((void **) &cuda_mprts->d_alt_pxi4, nr_alloced * sizeof(float4)); cudaCheck(ierr);
  ierr = cudaMalloc((void **) &cuda_mprts->d_bidx, nr_alloced * sizeof(unsigned int)); cudaCheck(ierr);
  ierr = cudaMalloc((void **) &cuda_mprts->d_id, nr_alloced * sizeof(unsigned int)); cudaCheck(ierr);

  ierr = cudaMalloc((void **) &cuda_mprts->d_nr_prts_by_patch, cuda_mprts->nr_patches * sizeof(unsigned int)); cudaCheck(ierr);

  ierr = cudaMalloc((void **) &cuda_mprts->d_off, (cuda_mprts->nr_blocks + 1) * sizeof(unsigned int)); cudaCheck(ierr);
}

// ----------------------------------------------------------------------
// cuda_mparticles_free

void
cuda_mparticles_free(struct cuda_mparticles *cuda_mprts)
{
  cudaError_t ierr;

  ierr = cudaFree(cuda_mprts->d_xi4); cudaCheck(ierr);
  ierr = cudaFree(cuda_mprts->d_pxi4); cudaCheck(ierr);
  ierr = cudaFree(cuda_mprts->d_alt_xi4); cudaCheck(ierr);
  ierr = cudaFree(cuda_mprts->d_alt_pxi4); cudaCheck(ierr);
  ierr = cudaFree(cuda_mprts->d_bidx); cudaCheck(ierr);
  ierr = cudaFree(cuda_mprts->d_id); cudaCheck(ierr);

  ierr = cudaFree(cuda_mprts->d_nr_prts_by_patch); cudaCheck(ierr);

  ierr = cudaFree(cuda_mprts->d_off); cudaCheck(ierr);
}

// ----------------------------------------------------------------------
// cuda_mparticles_dump

void
cuda_mparticles_dump(struct cuda_mparticles *cuda_mprts)
{
  int nr_prts = cuda_mprts->nr_prts;
  
  thrust::device_ptr<float4> d_xi4(cuda_mprts->d_xi4);
  thrust::device_ptr<float4> d_pxi4(cuda_mprts->d_pxi4);
  thrust::device_ptr<unsigned int> d_bidx(cuda_mprts->d_bidx);
  thrust::device_ptr<unsigned int> d_id(cuda_mprts->d_id);
  thrust::device_ptr<unsigned int> d_off(cuda_mprts->d_off);

  printf("cuda_mparticles_dump: nr_prts = %d\n", nr_prts); 
  for (int n = 0; n < nr_prts; n++) {
    float4 xi4 = d_xi4[n], pxi4 = d_pxi4[n];
    unsigned int bidx = d_bidx[n], id = d_id[n];
    printf("cuda_mparticles_dump: [%d] %g %g %g // %g // %g %g %g // %g || bidx %d id %d\n",
	   n, xi4.x, xi4.y, xi4.z, xi4.w, pxi4.x, pxi4.y, pxi4.z, pxi4.w,
	   bidx, id);
  }

  for (int b = 0; b <= cuda_mprts->nr_blocks; b++) {
    unsigned int off = d_off[b];
    printf("cuda_mparticles_dump: off[%d] = %d\n", b, off);
  }
}

// ----------------------------------------------------------------------
// cuda_mparticles_swap_alt

void
cuda_mparticles_swap_alt(struct cuda_mparticles *cuda_mprts)
{
  float4 *tmp_xi4 = cuda_mprts->d_alt_xi4;
  float4 *tmp_pxi4 = cuda_mprts->d_alt_pxi4;
  cuda_mprts->d_alt_xi4 = cuda_mprts->d_xi4;
  cuda_mprts->d_alt_pxi4 = cuda_mprts->d_pxi4;
  cuda_mprts->d_xi4 = tmp_xi4;
  cuda_mprts->d_pxi4 = tmp_pxi4;
}

// ----------------------------------------------------------------------
// cuda_params

struct cuda_params {
  unsigned int b_mx[3];
  float b_dxi[3];
};

static void
cuda_params_set(struct cuda_params *prm, const struct cuda_mparticles *cuda_mprts)
{
  for (int d = 0; d < 3; d++) {
    prm->b_mx[d]  = cuda_mprts->b_mx[d];
    prm->b_dxi[d] = cuda_mprts->b_dxi[d];
  }
}

static void
cuda_params_free(struct cuda_params *prm)
{
}

// ----------------------------------------------------------------------
// cuda_mparticles_find_block_indices_ids_total

#define THREADS_PER_BLOCK 512

__global__ static void
mprts_find_block_indices_ids_total(struct cuda_params prm, float4 *d_xi4,
				   unsigned int *d_nr_prts_by_patch,
				   unsigned int *d_bidx, unsigned int *d_id,
				   int nr_patches)
{
  int n = threadIdx.x + THREADS_PER_BLOCK * blockIdx.x;
  int nr_blocks = prm.b_mx[1] * prm.b_mx[2];

  unsigned int off = 0;
  for (int p = 0; p < nr_patches; p++) {
    if (n < d_nr_prts_by_patch[p]) {
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
      d_id[n + off] = n + off;
    }
    off += d_nr_prts_by_patch[p];
  }
}

void
cuda_mparticles_find_block_indices_ids_total(struct cuda_mparticles *cuda_mprts,
					     unsigned int *nr_prts_by_patch)
{
  cudaError_t ierr;
  
  if (cuda_mprts->nr_patches == 0) {
    return;
  }

  int max_nr_prts = 0;
  for (int p = 0; p < cuda_mprts->nr_patches; p++) {
    if (nr_prts_by_patch[p] > max_nr_prts) {
      max_nr_prts = nr_prts_by_patch[p];
    }
  }

  ierr = cudaMemcpy(cuda_mprts->d_nr_prts_by_patch, nr_prts_by_patch,
		    cuda_mprts->nr_patches * sizeof(unsigned int),
		    cudaMemcpyHostToDevice); cudaCheck(ierr);

  struct cuda_params prm;
  cuda_params_set(&prm, cuda_mprts);
    
  dim3 dimGrid((max_nr_prts + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK);
  dim3 dimBlock(THREADS_PER_BLOCK);

  mprts_find_block_indices_ids_total<<<dimGrid, dimBlock>>>(prm, cuda_mprts->d_xi4, 
							    cuda_mprts->d_nr_prts_by_patch,
							    cuda_mprts->d_bidx,
							    cuda_mprts->d_id,
							    cuda_mprts->nr_patches);
  cuda_sync_if_enabled();
  cuda_params_free(&prm);
}

// ----------------------------------------------------------------------
// cuda_mparticles_reorder_and_offsets

__global__ static void
mprts_reorder_and_offsets(int nr_prts, float4 *xi4, float4 *pxi4, float4 *alt_xi4, float4 *alt_pxi4,
			  unsigned int *d_bidx, unsigned int *d_ids, unsigned int *d_off, int last_block)
{
  int i = threadIdx.x + THREADS_PER_BLOCK * blockIdx.x;

  if (i > nr_prts)
    return;

  int block, prev_block;
  if (i < nr_prts) {
    alt_xi4[i] = xi4[d_ids[i]];
    alt_pxi4[i] = pxi4[d_ids[i]];
    
    block = d_bidx[i];
  } else { // needed if there is no particle in the last block
    block = last_block;
  }

  // OPT: d_bidx[i-1] could use shmem
  // create offsets per block into particle array
  prev_block = -1;
  if (i > 0) {
    prev_block = d_bidx[i-1];
  }
  for (int b = prev_block + 1; b <= block; b++) {
    d_off[b] = i;
  }
}

void
cuda_mparticles_reorder_and_offsets(struct cuda_mparticles *cuda_mprts)
{
  if (cuda_mprts->nr_patches == 0) {
    return;
  }

  dim3 dimGrid((cuda_mprts->nr_prts + 1 + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK);
  dim3 dimBlock(THREADS_PER_BLOCK);

  mprts_reorder_and_offsets<<<dimGrid, dimBlock>>>(cuda_mprts->nr_prts, cuda_mprts->d_xi4, cuda_mprts->d_pxi4,
						   cuda_mprts->d_alt_xi4, cuda_mprts->d_alt_pxi4,
						   cuda_mprts->d_bidx, cuda_mprts->d_id,
						   cuda_mprts->d_off, cuda_mprts->nr_blocks);
  cuda_sync_if_enabled();

  cuda_mparticles_swap_alt(cuda_mprts);
}



