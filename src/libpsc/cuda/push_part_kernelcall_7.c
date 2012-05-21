
// ----------------------------------------------------------------------
// cuda_push_part_p2
//
// particle push

EXTERN_C void
PFX(cuda_push_part_p2)(particles_cuda_t *pp, struct psc_fields *pf)
{
  struct psc_fields_cuda *pfc = psc_fields_cuda(pf);
  assert(pp->nr_blocks == pp->b_mx[1] * pp->b_mx[2]);

  int dimBlock[2] = { THREADS_PER_BLOCK, 1 };
  int dimGrid[2]  = { pp->b_mx[1], pp->b_mx[2] };

  RUN_KERNEL(dimGrid, dimBlock,
	     push_part_p1, (pp->n_part, pp->d_part, pfc->d_flds));
}

// ----------------------------------------------------------------------
// cuda_push_part_p3
//
// calculate currents

EXTERN_C void
PFX(cuda_push_part_p3)(particles_cuda_t *pp, struct psc_fields *pf, real *dummy,
		       int block_stride)
{
  struct psc_fields_cuda *pfc = psc_fields_cuda(pf);
  unsigned int size = pf->im[0] * pf->im[1] * pf->im[2];
  check(cudaMemset(pfx->d_flds + JXI * size, 0,
		   3 * size * sizeof(*pfc->d_flds)));

  assert(pp->nr_blocks % block_stride == 0);
  assert(block_stride == 4);
  assert(pp->nr_blocks == pp->b_mx[1] * pp->b_mx[2]);

  int dimBlock[2] = { THREADS_PER_BLOCK, 1 };
  int dimGrid[2]  = { (pp->b_mx[1] + 1) / 2, (pp->b_mx[2] + 1) / 2 };

  for (int block_start = 0; block_start < block_stride; block_start++) {
    RUN_KERNEL(dimGrid, dimBlock,
	       push_part_p2x, (pp->n_part, pp->d_part, pfc->d_flds,
			       block_stride, block_start));
  }
}

