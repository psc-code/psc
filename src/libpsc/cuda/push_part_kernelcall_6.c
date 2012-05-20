
// ----------------------------------------------------------------------
// cuda_push_part_p2
//
// particle push

EXTERN_C void
PFX(cuda_push_part_p2)(struct psc_particles *prts, fields_cuda_t *pf)
{
  struct psc_particles_cuda *cuda = psc_particles_cuda(prts);
  int dimBlock[2] = { THREADS_PER_BLOCK, 1 };
  int dimGrid[2]  = { cuda->nr_blocks, 1 };
  RUN_KERNEL(dimGrid, dimBlock,
	     push_part_p1, (prts->n_part, cuda->d_part, pf->d_flds));
}

// ----------------------------------------------------------------------
// cuda_push_part_p3
//
// calculate currents

EXTERN_C void
PFX(cuda_push_part_p3)(struct psc_particles *prts, fields_cuda_t *pf, real *dummy,
		       int block_stride)
{
  struct psc_particles_cuda *cuda = psc_particles_cuda(prts);
  unsigned int size = pf->im[0] * pf->im[1] * pf->im[2];
  check(cudaMemset(pf->d_flds + JXI * size, 0,
		   3 * size * sizeof(*pf->d_flds)));

  assert(cuda->nr_blocks % block_stride == 0);
  int dimBlock[2] = { THREADS_PER_BLOCK, 1 };
  int dimGrid[2]  = { cuda->nr_blocks / block_stride, 1 };

  for (int block_start = 0; block_start < block_stride; block_start++) {
    RUN_KERNEL(dimGrid, dimBlock,
	       push_part_p2x, (prts->n_part, cuda->d_part, pf->d_flds,
			       block_stride, block_start));

    RUN_KERNEL(dimGrid, dimBlock,
	       push_part_p2y, (prts->n_part, cuda->d_part, pf->d_flds,
			       block_stride, block_start));

    RUN_KERNEL(dimGrid, dimBlock,
	       push_part_p2z, (prts->n_part, cuda->d_part, pf->d_flds,
			       block_stride, block_start));
  }
}

