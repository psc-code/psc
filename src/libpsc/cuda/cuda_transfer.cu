
#include "psc_cuda.h"

EXTERN_C void
__particles_cuda_alloc(particles_cuda_t *pp, bool need_block_offsets,
		       bool need_cell_offsets)
{
  int n_part = pp->n_part;
  particles_cuda_dev_t *d_part = &pp->d_part;

  const int cells_per_block = BLOCKSIZE_X * BLOCKSIZE_Y * BLOCKSIZE_Z;

  check(cudaMalloc((void **) &d_part->xi4, n_part * sizeof(float4)));
  check(cudaMalloc((void **) &d_part->pxi4, n_part * sizeof(float4)));

  if (need_block_offsets) {
    check(cudaMalloc((void **) &d_part->offsets, 
		     (pp->nr_blocks + 1) * sizeof(int)));
    check(cudaMemcpy(&d_part->offsets[pp->nr_blocks], &n_part, sizeof(int),
		     cudaMemcpyHostToDevice));
  }

  if (need_cell_offsets) {
    check(cudaMalloc((void **) &d_part->c_offsets, 
		     (pp->nr_blocks * cells_per_block + 1) * sizeof(int)));
  }

  check(cudaMalloc((void **) &d_part->c_pos, 
		   (pp->nr_blocks * cells_per_block * 3) * sizeof(int)));
}

EXTERN_C void
__particles_cuda_to_device(particles_cuda_t *pp, float4 *xi4, float4 *pxi4,
			   int *offsets, int *c_offsets, int *c_pos)
{
  int n_part = pp->n_part;
  particles_cuda_dev_t *d_part = &pp->d_part;

  const int cells_per_block = BLOCKSIZE_X * BLOCKSIZE_Y * BLOCKSIZE_Z;

  check(cudaMemcpy(d_part->xi4, xi4, n_part * sizeof(*xi4),
		   cudaMemcpyHostToDevice));
  check(cudaMemcpy(d_part->pxi4, pxi4, n_part * sizeof(*pxi4),
		   cudaMemcpyHostToDevice));
  if (offsets) {
    check(cudaMemcpy(d_part->offsets, offsets,
		     pp->nr_blocks * sizeof(int), cudaMemcpyHostToDevice));
  }
  if (c_offsets) {
    check(cudaMemcpy(d_part->c_offsets,c_offsets,
		     (pp->nr_blocks * cells_per_block + 1) * sizeof(int),
		     cudaMemcpyHostToDevice));
  }
  check(cudaMemcpy(d_part->c_pos, c_pos,
		   (pp->nr_blocks * cells_per_block * 3) * sizeof(int),
		   cudaMemcpyHostToDevice));
}

EXTERN_C void
__particles_cuda_from_device(particles_cuda_t *pp, float4 *xi4, float4 *pxi4)
{
  int n_part = pp->n_part;
  particles_cuda_dev_t *d_part = &pp->d_part;

  check(cudaMemcpy(xi4, d_part->xi4, n_part * sizeof(*xi4),
		   cudaMemcpyDeviceToHost));
  check(cudaMemcpy(pxi4, d_part->pxi4, n_part * sizeof(*pxi4),
		   cudaMemcpyDeviceToHost));
}

EXTERN_C void
__particles_cuda_free(particles_cuda_t *pp)
{
  particles_cuda_dev_t *d_part = &pp->d_part;

  check(cudaFree(d_part->xi4));
  check(cudaFree(d_part->pxi4));
  check(cudaFree(d_part->offsets));
  check(cudaFree(d_part->c_offsets));
  check(cudaFree(d_part->c_pos));
}

EXTERN_C void
__fields_cuda_to_device(fields_cuda_t *pf, real *h_flds, int mb, int me)
{
  assert(!ppsc->domain.use_pml);

  unsigned int size = pf->im[0] * pf->im[1] * pf->im[2];
  check(cudaMalloc((void **) &pf->d_flds, pf->nr_comp * size * sizeof(float)));
  check(cudaMemcpy(pf->d_flds + mb * size,
		   h_flds + mb * size,
		   (me - mb) * size * sizeof(float),
		   cudaMemcpyHostToDevice));
}

EXTERN_C void
__fields_cuda_from_device(fields_cuda_t *pf, real *h_flds, int mb, int me)
{
  unsigned int size = pf->im[0] * pf->im[1] * pf->im[2];
  check(cudaMemcpy(h_flds + mb * size,
		   pf->d_flds + mb * size,
		   (me - mb) * size * sizeof(float),
		   cudaMemcpyDeviceToHost));
  check(cudaFree(pf->d_flds));
}
