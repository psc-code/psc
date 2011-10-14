
#include "psc_cuda.h"

#define BLOCKSIZE_X 1
#define BLOCKSIZE_Y 4
#define BLOCKSIZE_Z 4

#define DIM DIM_YZ
#define PFX(x) yz_a_ ## x

#include "constants.c"
#include "common.c"
#include "common_push.c"
#include "common_fld_cache.c"

__device__ static void
push_part_yz_a_one(int n, particles_cuda_dev_t d_particles, real *d_flds)
{
  struct d_particle p;
  LOAD_PARTICLE(p, d_particles, n);
  real vxi[3];

  // x^n, p^n -> x^(n+0.5), p^n
  
  calc_vxi(vxi, p);
  push_xi(&p, vxi, .5f * d_dt);

  STORE_PARTICLE_POS(p, d_particles, n);
  //  STORE_PARTICLE_MOM(p, d_particles, n);
}

__global__ static void
push_part_yz_a(int n_part, particles_cuda_dev_t d_part, float *d_flds, int stride)
{
  int n = threadIdx.x + blockDim.x * blockIdx.x;
  while (n < n_part) {
    push_part_yz_a_one(n, d_part, d_flds);
    n += stride;
  }
}

EXTERN_C void
__cuda_push_part_yz_a(particles_cuda_t *pp, fields_cuda_t *pf)
{
  const int threadsPerBlock = 128;
  const int gridSize = 256;
  int dimBlock[2]  = { threadsPerBlock, 1 };
  int dimGrid[2] = { gridSize, 1 };
  RUN_KERNEL(dimGrid, dimBlock,
	     push_part_yz_a, (pp->n_part, pp->d_part, pf->d_flds,
			      gridSize * threadsPerBlock));
}

