
#include <psc_cuda.h>

// FIXME, use const mem for dxi stuff (find_indices, too)
// probably should do our own loop rather than use blockIdx

__global__ static void
exchange_particles(int n_part, particles_cuda_dev_t d_part,
		   int ldimsx, int ldimsy, int ldimsz,
		   real _dxi, real _dyi, real _dzi)
{
  int ldims[3] = { ldimsx, ldimsy, ldimsz };
  real dxi[3] = { _dxi, _dyi, _dzi };
  int xm[3];

  for (int d = 0; d < 3; d++) {
    xm[d] = ldims[d] / dxi[d];
  }

  int i = threadIdx.x + THREADS_PER_BLOCK * blockIdx.x;
  if (i < n_part) {
    particle_cuda_real_t xi[3] = {
      d_part.xi4[i].x * dxi[0],
      d_part.xi4[i].y * dxi[1],
      d_part.xi4[i].z * dxi[2] };
    int pos[3];
    for (int d = 0; d < 3; d++) {
      pos[d] = cuda_fint(xi[d]);
    }
    if (pos[1] < 0) {
      d_part.xi4[i].y += xm[1];
    }
    if (pos[2] < 0) {
      d_part.xi4[i].z += xm[2];
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
  real dxi[3];
  for (int d = 0; d < 3; d++) {
    dxi[d] = 1.f / ppsc->dx[d];
  }

  int dimBlock[2] = { THREADS_PER_BLOCK, 1 };
  int dimGrid[2]  = { (pp->n_part + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK, 1 };
  RUN_KERNEL(dimGrid, dimBlock,
	     exchange_particles, (pp->n_part, pp->d_part,
				  patch->ldims[0], patch->ldims[1], patch->ldims[2],
				  dxi[0], dxi[1], dxi[2]));
  
}
