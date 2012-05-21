
// ----------------------------------------------------------------------
// cuda_push_part_p1
//
// alloc scratch

EXTERN_C void
PFX(cuda_push_part_p1)(struct psc_particles *prts, struct psc_fields *pf,
		       real **d_scratch)
{
  struct psc_fields_cuda *pfc = psc_fields_cuda(pf);
  struct psc_particles_cuda *cuda = psc_particles_cuda(prts);
  unsigned int size = cuda->nr_blocks * 3 * BLOCKSTRIDE;
  check(cudaMalloc((void **)d_scratch, size * sizeof(real)));
  check(cudaMemset(*d_scratch, 0, size * sizeof(real)));
  size = pf->im[0] * pf->im[1] * pf->im[2];
  check(cudaMemset(pfc->d_flds + JXI * size, 0,
		   3 * size * sizeof(*pfc->d_flds)));

  cuda_sync_if_enabled();
}

// ----------------------------------------------------------------------
// cuda_push_part_p2
//
// particle push

EXTERN_C void
PFX(cuda_push_part_p2)(struct psc_particles *prts, struct psc_fields *pf)
{
  struct psc_fields_cuda *pfc = psc_fields_cuda(pf);
  struct psc_particles_cuda *cuda = psc_particles_cuda(prts);
  int dimBlock[2] = { THREADS_PER_BLOCK, 1 };
  int dimGrid[2]  = { cuda->nr_blocks, 1 };
  RUN_KERNEL(dimGrid, dimBlock,
	     push_part_p1, (prts->n_part, cuda->d_part, pfc->d_flds));
}

// ----------------------------------------------------------------------
// cuda_push_part_p3
//
// calculate currents locally

EXTERN_C void
PFX(cuda_push_part_p3)(struct psc_particles *prts, struct psc_fields *pf, real *d_scratch,
		       int block_stride)
{
  struct psc_fields_cuda *pfc = psc_fields_cuda(pf);
  struct psc_particles_cuda *cuda = psc_particles_cuda(prts);
  assert(cuda->nr_blocks % block_stride == 0);
  int dimBlock[2] = { THREADS_PER_BLOCK, 1 };
  int dimGrid[2]  = { cuda->nr_blocks / block_stride, 1 };
  for (int block_start = 0; block_start < block_stride; block_start++) {
    RUN_KERNEL(dimGrid, dimBlock,
	       push_part_p2, (prts->n_part, cuda->d_part, pfc->d_flds, d_scratch,
			      block_stride, block_start));
  }
}

// ----------------------------------------------------------------------
// cuda_push_part_p4
//
// collect calculation

EXTERN_C void
PFX(cuda_push_part_p4)(struct psc_particles *prts, struct psc_fields *pf, real *d_scratch)
{
  struct psc_fields_cuda *pfc = psc_fields_cuda(pf);
  struct psc_particles_cuda *cuda = psc_particles_cuda(prts);
#if DIM == DIM_Z
  int dimBlock[2] = { BLOCKSIZE_Z + 2*SW, 1 };
#elif DIM == DIM_YZ
  int dimBlock[2]  = { BLOCKSIZE_Y + 2*SW, BLOCKSIZE_Z + 2*SW };
#endif
  int dimGrid[2]  = { 1, 1 };
  RUN_KERNEL(dimGrid, dimBlock,
	     collect_currents, (pfc->d_flds, d_scratch, cuda->nr_blocks));

  int sz = cuda->nr_blocks * 3 * BLOCKSTRIDE;
  real *h_scratch = (real *) malloc(sz * sizeof(real));
  check(cudaMemcpy(h_scratch, d_scratch, sz * sizeof(real),
		   cudaMemcpyDeviceToHost));
#if 0
  // version which does the calculation on the host (z-only)
#define h_scratch(m,jy,jz) (p[(m)*BLOCKSTRIDE +		\
			      ((jz)+3) * (BLOCKSIZE_Y+6) +	\
			      (jy)+3])
  for (int m = 0; m < 3; m++) {
    printf("m %d\n", m);
    for (int b = 0; b < cuda->nr_blocks; b++) {
      real *p = h_scratch + b * 3 * BLOCKSTRIDE;
      for (int iz = -3; iz <= 3; iz++) {
	if (h_scratch(m, 0, iz) != 0.) {
	  printf("b %d iz %d : %g\n", b, iz, h_scratch(m, 0, iz));
	}
      }
    }
  }
  free(h_scratch);

  for (int b = 0; b < cuda->nr_blocks; b++) {
    real *p = h_scratch + b * 3 * BLOCKSTRIDE;
    for (int m = 0; m < 3; m++) {
      for (int iz = -3; iz <= 3; iz++) {
	F3_CUDA(pf, m, 0,0,b+iz) += h_scratch(m, 0,iz);
      }
    }
  }
  for (int iz = -3; iz < 10+3; iz++) {
    printf("iz %d: %g\n", iz, F3_CUDA(pf, 2, 0,0,iz));
  }
#endif
}

// ----------------------------------------------------------------------
// cuda_push_part_p5
//
// free

EXTERN_C void
PFX(cuda_push_part_p5)(struct psc_particles *prts, struct psc_fields *pf, real *d_scratch)
{
  check(cudaFree(d_scratch));
}

