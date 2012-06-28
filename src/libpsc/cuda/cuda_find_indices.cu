
#include <psc_cuda.h>

// FIXME, do this always?
#define NO_CHECKERBOARD

#define BLOCKSIZE_X 1
#define BLOCKSIZE_Y 4
#define BLOCKSIZE_Z 4

#define PFX(x) sort_##x
#include "constants.c"

// FIXME, use const mem for some params

__global__ static void find_cell_indices_by_cell(int n_part, particles_cuda_dev_t d_part,
						 int *d_cnis, int *d_ids, int ldims_y)
{
  int i = threadIdx.x + THREADS_PER_BLOCK * blockIdx.x;
  if (i < n_part) {
    particle_cuda_real_t xi[3] = {
      d_part.xi4[i].x * d_consts.dxi[0],
      d_part.xi4[i].y * d_consts.dxi[1],
      d_part.xi4[i].z * d_consts.dxi[2] };
    int pos[3];
    for (int d = 0; d < 3; d++) {
      pos[d] = cuda_fint(xi[d]);
    }
    
    int idx = (((pos[2] / 8) * (ldims_y / 8) + (pos[1] / 8)) << 6);
    idx |=
      ((pos[2] & 4) << 3) |
      ((pos[1] & 4) << 2);
    idx |=
      ((pos[2] & 2) << 2) |
      ((pos[1] & 2) << 1) |
      ((pos[2] & 1) << 1) |
      ((pos[1] & 1) << 0);
    d_cnis[i] = idx;
    d_ids[i] = i;
  }
}

static void
sort_find_cell_indices_by_cell_device(particles_cuda_t *pp, struct psc_patch *patch,
				      int *d_cnis, int *d_ids)
{
  int dimBlock[2] = { THREADS_PER_BLOCK, 1 };
  int dimGrid[2]  = { (pp->n_part + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK, 1 };
  RUN_KERNEL(dimGrid, dimBlock,
	     find_cell_indices_by_cell, (pp->n_part, pp->d_part, d_cnis, d_ids,
					 patch->ldims[1]));
}

// FIXME, specific to 1x8x8, should be in ! .cu, so that cell_map works

static void __unused
sort_find_cell_indices_host(particles_cuda_t *pp, struct psc_patch *patch,
			    int *d_cnis, int *d_ids)
{
  int n_part = pp->n_part;
  particles_cuda_dev_t *d_part = &pp->d_part;
  int *h_cnis = (int *) malloc(n_part * sizeof(*h_cnis));
  int *h_ids = (int *) malloc(n_part * sizeof(*h_ids));

  float4 *h_xi4 = (float4 *) malloc(n_part * sizeof(*h_xi4));
  check(cudaMemcpy(h_xi4, d_part->xi4, n_part * sizeof(float4),
		   cudaMemcpyDeviceToHost));

  for (int i = 0; i < n_part; i++) {
    particle_cuda_real_t dxi = 1.f / ppsc->dx[0];
    particle_cuda_real_t dyi = 1.f / ppsc->dx[1];
    particle_cuda_real_t dzi = 1.f / ppsc->dx[2];
    particle_cuda_real_t xi[3] = {
      (h_xi4[i].x - patch->xb[0]) * dxi,
      (h_xi4[i].y - patch->xb[1]) * dyi,
      (h_xi4[i].z - patch->xb[2]) * dzi };
    int pos[3];
    for (int d = 0; d < 3; d++) {
      pos[d] = particle_cuda_real_fint(xi[d]);
    }
    
    int idx = (((pos[2] / 8) * (patch->ldims[1] / 8) + (pos[1] / 8)) << 6) |
      ((pos[2] & 4) << 3) |
      ((pos[2] & 2) << 2) |
      ((pos[2] & 1) << 1) |
      ((pos[1] & 4) << 2) |
      ((pos[1] & 2) << 1) |
      ((pos[1] & 1) << 0);
    h_cnis[i] = idx;
    h_ids[i] = i;
  }

  check(cudaMemcpy(d_cnis, h_cnis, n_part * sizeof(*h_cnis),
		   cudaMemcpyHostToDevice));
  check(cudaMemcpy(d_ids, h_ids, n_part * sizeof(*h_ids),
		   cudaMemcpyHostToDevice));

  free(h_xi4);
  free(h_cnis);
  free(h_ids);
}

static void __unused
sort_reorder_host(particles_cuda_t *pp, int *d_ids)
{
  int n_part = pp->n_part;
  particles_cuda_dev_t *d_part = &pp->d_part;
  int *h_ids = (int *) malloc(n_part);

  float4 *h_xi4 = (float4 *) malloc(n_part * sizeof(*h_xi4));
  float4 *h_pxi4 = (float4 *) malloc(n_part * sizeof(*h_pxi4));
  check(cudaMemcpy(h_xi4, d_part->xi4, n_part * sizeof(float4),
		   cudaMemcpyDeviceToHost));
  check(cudaMemcpy(h_pxi4, d_part->pxi4, n_part * sizeof(float4),
		   cudaMemcpyDeviceToHost));
  check(cudaMemcpy(h_ids, d_ids, n_part * sizeof(*h_ids),
		   cudaMemcpyDeviceToHost));

  // move into new position
  float4 *xi4 = (float4 *) malloc(n_part * sizeof(*xi4));
  float4 *pxi4 = (float4 *) malloc(n_part * sizeof(*pxi4));
  for (int i = 0; i < n_part; i++) {
    xi4[i] = h_xi4[h_ids[i]];
    pxi4[i] = h_pxi4[h_ids[i]];
  }
  
  check(cudaMemcpy(d_part->xi4, xi4, n_part * sizeof(float4),
		   cudaMemcpyHostToDevice));
  check(cudaMemcpy(d_part->pxi4, pxi4, n_part * sizeof(float4),
		   cudaMemcpyHostToDevice));
  free(xi4);
  free(pxi4);
  free(h_xi4);
  free(h_pxi4);
  free(h_ids);
}

__global__ static void
sort_reorder_by_cell(int n_part, particles_cuda_dev_t d_part, float4 *xi4, float4 *pxi4,
		     int *d_cnis, int *d_ids)
{
  int i = threadIdx.x + THREADS_PER_BLOCK * blockIdx.x;
  int blocksize = BLOCKSIZE_X * BLOCKSIZE_Y * BLOCKSIZE_Z;

  if (i > n_part)
    return;

  int block, prev_block;
  if (i < n_part) {
    xi4[i] = d_part.xi4[d_ids[i]];
    pxi4[i] = d_part.pxi4[d_ids[i]];
    
    block = d_cnis[i];
  } else if (i == n_part) { // needed if there is no particle in the last block
    block = d_consts.b_mx[0] * d_consts.b_mx[1] * d_consts.b_mx[2] * blocksize;
  }

  // create offsets per block into particle array
  prev_block = -1;
  if (i > 0) {
    prev_block = d_cnis[i-1];
  }
  for (int b = prev_block + 1; b <= block; b++) {
    d_part.c_offsets[b] = i;
  }
}

static void
sort_reorder_by_cell_device(particles_cuda_t *pp, int *d_cnis, int *d_ids)
{
  float4 *xi4, *pxi4;
  check(cudaMalloc((void **) &xi4, pp->n_alloced * sizeof(*xi4)));
  check(cudaMalloc((void **) &pxi4, pp->n_alloced * sizeof(*pxi4)));

  int dimBlock[2] = { THREADS_PER_BLOCK, 1 };
  int dimGrid[2]  = { (pp->n_part + 1 + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK, 1 };
  RUN_KERNEL(dimGrid, dimBlock,
	     sort_reorder_by_cell, (pp->n_part, pp->d_part, xi4, pxi4, d_cnis, d_ids));

  check(cudaFree(pp->d_part.xi4));
  check(cudaFree(pp->d_part.pxi4));
  pp->d_part.xi4 = xi4;
  pp->d_part.pxi4 = pxi4;
}

EXTERN_C void
sort_patch_prep(int p, particles_cuda_t *pp, int **d_cnis, int **d_ids)
{
  check(cudaMalloc((void **) d_cnis, pp->n_alloced * sizeof(*d_cnis)));
  check(cudaMalloc((void **) d_ids, pp->n_alloced * sizeof(*d_ids)));

  fields_cuda_t pf_dummy;
  sort_set_constants(pp, &pf_dummy);
}

EXTERN_C void
sort_patch_done(int p, particles_cuda_t *pp, int *d_cnis, int *d_ids)
{
  check(cudaFree(d_cnis));
  check(cudaFree(d_ids));
}

EXTERN_C void
sort_patch_by_cell(int p, particles_cuda_t *pp)
{
  struct psc_patch *patch = &ppsc->patch[p];

  int *d_cnis, *d_ids;
  check(cudaMalloc((void **) &d_cnis, pp->n_part * sizeof(*d_cnis)));
  check(cudaMalloc((void **) &d_ids, pp->n_part * sizeof(*d_ids)));

  fields_cuda_t pf_dummy;
  sort_set_constants(pp, &pf_dummy);

  sort_find_cell_indices_by_cell_device(pp, patch, d_cnis, d_ids);
  sort_pairs_device((unsigned int *) d_cnis, (unsigned int *) d_ids, pp->n_part);
  sort_reorder_by_cell_device(pp, d_cnis, d_ids);
  
  check(cudaFree(d_cnis));
  check(cudaFree(d_ids));
}

