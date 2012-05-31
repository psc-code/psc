
#include "psc_cuda.h"

// FIXME, merge with F3_DEV{,_YZ}, OPT (precalc offset)

#define X3_DEV_OFF_YZ(fldnr, jy,jz)					\
  ((((fldnr)								\
     *mz + ((jz)+3))							\
    *my + ((jy)+3))							\
   *7 + (0+3))

#undef F3_DEV

#define F3_DEV(fldnr,ix,jy,jz)			\
  (d_flds)[X3_DEV_OFF_YZ(fldnr, jy,jz)]


__global__ static void
push_fields_E_yz(real *d_flds, real dt, real cny, real cnz, int my, int mz)
{
  int iy = blockIdx.x * blockDim.x + threadIdx.x;
  int iz = blockIdx.y * blockDim.y + threadIdx.y;

  if (!(iy < my && iz < mz))
    return;

  F3_DEV(EX, 0,iy,iz) +=
    cny * (F3_DEV(HZ, 0,iy,iz) - F3_DEV(HZ, 0,iy-1,iz)) -
    cnz * (F3_DEV(HY, 0,iy,iz) - F3_DEV(HY, 0,iy,iz-1)) -
    .5f * dt * F3_DEV(JXI, 0,iy,iz);
  
  F3_DEV(EY, 0,iy,iz) +=
    cnz * (F3_DEV(HX, 0,iy,iz) - F3_DEV(HX, 0,iy,iz-1)) -
    0.f -
    .5f * dt * F3_DEV(JYI, 0,iy,iz);
  
  F3_DEV(EZ, 0,iy,iz) +=
    0.f -
    cny * (F3_DEV(HX, 0,iy,iz) - F3_DEV(HX, 0,iy-1,iz)) -
    .5f * dt * F3_DEV(JZI, 0,iy,iz);
}

__global__ static void
push_fields_H_yz(real *d_flds, real cny, real cnz, int my, int mz)
{
  int iy = blockIdx.x * blockDim.x + threadIdx.x;
  int iz = blockIdx.y * blockDim.y + threadIdx.y;

  if (!(iy < my && iz < mz))
    return;

  F3_DEV(HX, 0,iy,iz) -=
    cny * (F3_DEV(EZ, 0,iy+1,iz) - F3_DEV(EZ, 0,iy,iz)) -
    cnz * (F3_DEV(EY, 0,iy,iz+1) - F3_DEV(EY, 0,iy,iz));
  
  F3_DEV(HY, 0,iy,iz) -=
    cnz * (F3_DEV(EX, 0,iy,iz+1) - F3_DEV(EX, 0,iy,iz)) -
    0.f;
  
  F3_DEV(HZ, 0,iy,iz) -=
    0.f -
    cny * (F3_DEV(EX, 0,iy+1,iz) - F3_DEV(EX, 0,iy,iz));
}

#define BLOCKSIZE_X 1
#define BLOCKSIZE_Y 16
#define BLOCKSIZE_Z 16

EXTERN_C void
cuda_push_fields_E_yz(int p, struct psc_fields *pf)
{
  struct psc_fields_cuda *pfc = psc_fields_cuda(pf);
  struct psc_patch *patch = &ppsc->patch[p];

  real dt = ppsc->dt;
  real cny = .5f * ppsc->dt / ppsc->dx[1];
  real cnz = .5f * ppsc->dt / ppsc->dx[2];
  assert(pf->ib[0] == -3);
  assert(pf->ib[1] == -3);
  assert(pf->ib[2] == -3);
  assert(pf->im[0] == 7);
  int my = pf->im[1];
  int mz = pf->im[2];

  int dimBlock[2] = { BLOCKSIZE_Y, BLOCKSIZE_Z };
  int dimGrid[2]  = { (patch->ldims[1] + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y,
		      (patch->ldims[2] + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z };
  RUN_KERNEL(dimGrid, dimBlock,
	     push_fields_E_yz, (pfc->d_flds, dt, cny, cnz, my, mz));
}

EXTERN_C void
cuda_push_fields_H_yz(int p, struct psc_fields *pf)
{
  struct psc_fields_cuda *pfc = psc_fields_cuda(pf);
  struct psc_patch *patch = &ppsc->patch[p];

  real cny = .5f * ppsc->dt / ppsc->dx[1];
  real cnz = .5f * ppsc->dt / ppsc->dx[2];
  int my = pf->im[1];
  int mz = pf->im[2];

  int dimBlock[2] = { BLOCKSIZE_Y, BLOCKSIZE_Z };
  int dimGrid[2]  = { (patch->ldims[1] + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y,
		      (patch->ldims[2] + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z };
  RUN_KERNEL(dimGrid, dimBlock,
	     push_fields_H_yz, (pfc->d_flds, cny, cnz, my, mz));
}

