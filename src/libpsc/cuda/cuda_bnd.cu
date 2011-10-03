
#include "psc_cuda.h"

#define PFX(x) cuda_bnd_##x
#include "constants.c"

#define SW (3)

__global__ static void
fill_ghosts(real *d_flds, int mb, int me)
{
  int iy = blockIdx.x * blockDim.x + threadIdx.x;
  int iz = blockIdx.y * blockDim.y + threadIdx.y;

  if (!(iy < d_mx[1] && iz < d_mx[2]))
    return;

  bool inside = true;
  int jy = iy, jz = iz;
  if (jy < SW           ) { jy += d_mx[1] - 2*SW; inside = false; }
  if (jy >= d_mx[1] - SW) { jy -= d_mx[1] - 2*SW; inside = false; }
  if (jz < SW           ) { jz += d_mx[2] - 2*SW; inside = false; }
  if (jz >= d_mx[2] - SW) { jz -= d_mx[2] - 2*SW; inside = false; }

  if (inside)
    return;

  for (int m = mb; m < me; m++) {
    F3_DEV(m, 0,iy-SW,iz-SW) = F3_DEV(m, 0,jy-SW,jz-SW);
  }
}

EXTERN_C void
cuda_fill_ghosts(int p, fields_cuda_t *pf, int mb, int me)
{
  particles_cuda_t pp;
  cuda_bnd_set_constants(&pp, pf);

  struct psc_patch *patch = &ppsc->patch[p];
  int dimBlock[2] = { BLOCKSIZE_Y, BLOCKSIZE_Z };
  int dimGrid[2]  = { (patch->ldims[1] + 2*SW + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y,
		      (patch->ldims[2] + 2*SW + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z };
  RUN_KERNEL(dimGrid, dimBlock,
	     fill_ghosts, (pf->d_flds, mb, me));
}

