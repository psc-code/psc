
#include "psc_cuda.h"

EXTERN_C void
__particles_cuda_get(particles_cuda_t *pp)
{
  int n_part = pp->n_part;
  particles_cuda_dev_t *h_part = &pp->h_part;
  particles_cuda_dev_t *d_part = &pp->d_part;

  const int cells_per_block = BLOCKSIZE_X * BLOCKSIZE_Y * BLOCKSIZE_Z;
  check(cudaMalloc((void **) &d_part->xi4, n_part * sizeof(float4)));
  check(cudaMalloc((void **) &d_part->pxi4, n_part * sizeof(float4)));
  check(cudaMalloc((void **) &d_part->offsets, 
		   (pp->nr_blocks + 1) * sizeof(int)));
  check(cudaMalloc((void **) &d_part->c_pos, 
		   (pp->nr_blocks * cells_per_block * 3) * sizeof(int)));

  check(cudaMemcpy(d_part->xi4, h_part->xi4, n_part * sizeof(float4),
		   cudaMemcpyHostToDevice));
  check(cudaMemcpy(d_part->pxi4, h_part->pxi4, n_part * sizeof(float4),
		   cudaMemcpyHostToDevice));
  if (h_part->offsets) {
    check(cudaMemcpy(d_part->offsets, h_part->offsets,
		     pp->nr_blocks * sizeof(int), cudaMemcpyHostToDevice));
  }
  check(cudaMemcpy(&d_part->offsets[pp->nr_blocks], &n_part, sizeof(int),
		   cudaMemcpyHostToDevice));
  if (h_part->c_offsets) {
    check(cudaMalloc((void **) &d_part->c_offsets, 
		     (pp->nr_blocks * cells_per_block + 1) * sizeof(int)));
    check(cudaMemcpy(d_part->c_offsets, h_part->c_offsets,
		     (pp->nr_blocks * cells_per_block + 1) * sizeof(int),
		     cudaMemcpyHostToDevice));
  }
  check(cudaMemcpy(d_part->c_pos, h_part->c_pos,
		   (pp->nr_blocks * cells_per_block * 3) * sizeof(int),
		   cudaMemcpyHostToDevice));
}

EXTERN_C void
__particles_cuda_put(particles_cuda_t *pp)
{
  int n_part = pp->n_part;
  particles_cuda_dev_t *h_part = &pp->h_part;
  particles_cuda_dev_t *d_part = &pp->d_part;

  check(cudaMemcpy(h_part->xi4, d_part->xi4, n_part * sizeof(float4),
		   cudaMemcpyDeviceToHost));
  check(cudaMemcpy(h_part->pxi4, d_part->pxi4, n_part * sizeof(float4),
		   cudaMemcpyDeviceToHost));
  check(cudaFree(d_part->xi4));
  check(cudaFree(d_part->pxi4));
  check(cudaFree(d_part->offsets));
  check(cudaFree(d_part->c_offsets));
  check(cudaFree(d_part->c_pos));
}

EXTERN_C void
__fields_cuda_get(fields_cuda_t *pf)
{
  assert(!ppsc->domain.use_pml);

  unsigned int size = pf->im[0] * pf->im[1] * pf->im[2];
  check(cudaMalloc((void **) &pf->d_flds, pf->nr_comp * size * sizeof(float)));
  check(cudaMemcpy(pf->d_flds + EX * size,
		   pf->h_flds + EX * size,
		   6 * size * sizeof(float),
		   cudaMemcpyHostToDevice));
}

EXTERN_C void
__fields_cuda_put(fields_cuda_t *pf)
{
  unsigned int size = pf->im[0] * pf->im[1] * pf->im[2];
  check(cudaMemcpy(pf->h_flds + JXI * size,
		   pf->d_flds + JXI * size,
		   3 * size * sizeof(float),
		   cudaMemcpyDeviceToHost));
  check(cudaFree(pf->d_flds));
}
