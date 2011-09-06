
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
// calculate currents

EXTERN_C void
PFX(cuda_push_part_p3)(particles_cuda_t *pp, fields_cuda_t *pf, real *dummy,
		       int block_stride)
{
  struct shapeinfo_h *d_shapeinfo_h;
  struct shapeinfo_i *d_shapeinfo_i;
  check(cudaMalloc((void **)&d_shapeinfo_h, pp->n_part * sizeof(*d_shapeinfo_h)));
  check(cudaMalloc((void **)&d_shapeinfo_i, pp->n_part * sizeof(*d_shapeinfo_i)));
  real *d_vxi;
  check(cudaMalloc((void **)&d_vxi, pp->n_part * sizeof(*d_vxi)));
  real *d_qni;
  check(cudaMalloc((void **)&d_qni, pp->n_part * sizeof(*d_qni)));
  int *d_ci1;
  check(cudaMalloc((void **)&d_ci1, pp->n_part * sizeof(*d_ci1)));


  unsigned int size = pf->im[0] * pf->im[1] * pf->im[2];
  check(cudaMemset(pf->d_flds + JXI * size, 0,
		   3 * size * sizeof(*pf->d_flds)));

  assert(pp->nr_blocks % block_stride == 0);
  int dimBlock[2] = { THREADS_PER_BLOCK, 1 };
  int dimGrid[2]  = { pp->nr_blocks / block_stride, 1 };

  for (int block_start = 0; block_start < block_stride; block_start++) {
    RUN_KERNEL(dimGrid, dimBlock,
	       push_part_p1_5, (pp->n_part, pp->d_part,
				d_shapeinfo_h, d_shapeinfo_i, d_vxi, d_qni, d_ci1,
				block_stride, block_start));

    RUN_KERNEL(dimGrid, dimBlock,
	       push_part_p2x, (pp->n_part, pp->d_part,
			       d_shapeinfo_h, d_shapeinfo_i, d_vxi, d_qni, d_ci1,
			       pf->d_flds, block_stride, block_start));

    RUN_KERNEL(dimGrid, dimBlock,
	       push_part_p2y, (pp->n_part, pp->d_part,
			       d_shapeinfo_h, d_shapeinfo_i, d_vxi, d_qni, d_ci1,
			       pf->d_flds, block_stride, block_start));

    RUN_KERNEL(dimGrid, dimBlock,
	       push_part_p2z, (pp->n_part, pp->d_part,
			       d_shapeinfo_h, d_shapeinfo_i, d_vxi, d_qni, d_ci1,
			       pf->d_flds, block_stride, block_start));
  }

  check(cudaFree(d_shapeinfo_h));
  check(cudaFree(d_shapeinfo_i));
  check(cudaFree(d_vxi));
  check(cudaFree(d_qni));
  check(cudaFree(d_ci1));
}

