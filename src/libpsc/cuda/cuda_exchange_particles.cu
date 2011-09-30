
#include <psc_cuda.h>

#define PFX(x) xchg_##x
#include "constants.c"

// FIXME const mem for dims?
// FIXME probably should do our own loop rather than use blockIdx

__global__ static void
exchange_particles(int n_part, particles_cuda_dev_t d_part,
		   int ldimsx, int ldimsy, int ldimsz)
{
  int ldims[3] = { ldimsx, ldimsy, ldimsz };
  int xm[3];

  for (int d = 0; d < 3; d++) {
    xm[d] = ldims[d] / d_dxi[d];
  }

  int i = threadIdx.x + THREADS_PER_BLOCK * blockIdx.x;
  if (i < n_part) {
    particle_cuda_real_t xi[3] = {
      d_part.xi4[i].x * d_dxi[0],
      d_part.xi4[i].y * d_dxi[1],
      d_part.xi4[i].z * d_dxi[2] };
    int pos[3];
    for (int d = 0; d < 3; d++) {
      pos[d] = cuda_fint(xi[d]);
    }
    if (pos[1] < 0) {
      d_part.xi4[i].y += xm[1];
      if (d_part.xi4[i].y >= xm[1])
	d_part.xi4[i].y = 0.f;
    }
    if (pos[2] < 0) {
      d_part.xi4[i].z += xm[2];
      if (d_part.xi4[i].z >= xm[2])
	d_part.xi4[i].z = 0.f;
    }
    if (pos[1] >= ldims[1]) {
      d_part.xi4[i].y -= xm[1];
    }
    if (pos[2] >= ldims[2]) {
      d_part.xi4[i].z -= xm[2];
    }
  }
}

EXTERN_C void
cuda_exchange_particles(int p, particles_cuda_t *pp)
{
  struct psc_patch *patch = &ppsc->patch[p];

  fields_cuda_t pf_dummy;
  xchg_set_constants(pp, &pf_dummy);

  int dimBlock[2] = { THREADS_PER_BLOCK, 1 };
  int dimGrid[2]  = { (pp->n_part + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK, 1 };
  RUN_KERNEL(dimGrid, dimBlock,
	     exchange_particles, (pp->n_part, pp->d_part,
				  patch->ldims[0], patch->ldims[1], patch->ldims[2]));
  
}
