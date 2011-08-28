
// ----------------------------------------------------------------------
// cuda_push_part_p1
//
// alloc scratch

EXTERN_C void
PFX(cuda_push_part_p1)(particles_cuda_t *pp, fields_cuda_t *pf,
		       real **d_scratch)
{
  int size = pp->nr_blocks * 3 * BLOCKSTRIDE;
  check(cudaMalloc((void **)d_scratch, size * sizeof(real)));
  check(cudaMemset(*d_scratch, 0, size * sizeof(real)));
  cudaThreadSynchronize(); // FIXME
}

// ----------------------------------------------------------------------
// cuda_push_part_p2
//
// particle push

EXTERN_C void
PFX(cuda_push_part_p2)(particles_cuda_t *pp, fields_cuda_t *pf)
{
  int dimBlock[2] = { THREADS_PER_BLOCK, 1 };
  int dimGrid[2]  = { pp->nr_blocks, 1 };
  RUN_KERNEL(dimGrid, dimBlock,
	     push_part_p1, (pp->n_part, pp->d_part, pf->d_flds));
}

// ----------------------------------------------------------------------
// cuda_push_part_p3
//
// calculate currents locally

EXTERN_C void
PFX(cuda_push_part_p3)(particles_cuda_t *pp, fields_cuda_t *pf, real *d_scratch)
{
  int dimBlock[2] = { THREADS_PER_BLOCK, 1 };
  int dimGrid[2]  = { pp->nr_blocks, 1 };
  RUN_KERNEL(dimGrid, dimBlock,
	     push_part_p2, (pp->n_part, pp->d_part, pf->d_flds, d_scratch));
}

// ----------------------------------------------------------------------
// cuda_push_part_p4
//
// collect calculation

EXTERN_C void
PFX(cuda_push_part_p4)(particles_cuda_t *pp, fields_cuda_t *pf, real *d_scratch)
{
  check(cudaMemset(pf->d_flds + JXI * psc.fld_size, 0,
		   3 * psc.fld_size * sizeof(*pf->d_flds)));

#if DIM == DIM_Z
  int dimBlock[2] = { BLOCKSIZE_Z + 2*SW, 1 };
#elif DIM == DIM_YZ
  int dimBlock[2]  = { BLOCKSIZE_Y + 2*SW, BLOCKSIZE_Z + 2*SW };
#endif
  int dimGrid[2]  = { 1, 1 };
  RUN_KERNEL(dimGrid, dimBlock,
	     collect_currents, (pf->d_flds, d_scratch, pp->nr_blocks));
}

// ----------------------------------------------------------------------
// cuda_push_part_p5
//
// free

EXTERN_C void
PFX(cuda_push_part_p5)(particles_cuda_t *pp, fields_cuda_t *pf, real *d_scratch)
{
  check(cudaFree(d_scratch));
}

