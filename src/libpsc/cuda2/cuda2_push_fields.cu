
#include "psc_cuda2.h"

#define BLOCKSIZE_X 1
#define BLOCKSIZE_Y 16
#define BLOCKSIZE_Z 16

#include "psc_fields_cuda2.h"

// ----------------------------------------------------------------------
// FIXME
#include "cuda_wrap.h"

#define BND (2)

#define X3_DEV_OFF_YZ(fldnr, jy,jz)					\
  ((((fldnr)								\
     *mz + ((jz)+2))							\
    *my + ((jy)+2))							\
   *1 + (0))

#undef F3_DEV

#define F3_DEV(fldnr,ix,jy,jz)			\
  (d_flds)[X3_DEV_OFF_YZ(fldnr, jy,jz)]

// FIXME end
// ----------------------------------------------------------------------

__global__ static void
push_fields_E_yz(real *d_flds0, real dt, real cny, real cnz, int my, int mz,
		 unsigned int size, int gridy)
{
  int bidx_y = blockIdx.y % gridy;
  int p = blockIdx.y / gridy;
  int iy = blockIdx.x * blockDim.x + threadIdx.x;
  int iz = bidx_y * blockDim.y + threadIdx.y;

  if (!(iy < my - 2 * (2-BND) && iz < mz - 2 * (2-BND)))
    return;
  iy -= BND;
  iz -= BND;

  real *d_flds = d_flds0 + p * size;

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
push_fields_H_yz(real *d_flds0, real cny, real cnz, int my, int mz,
		 unsigned int size, int gridy)
{
  int bidx_y = blockIdx.y % gridy;
  int p = blockIdx.y / gridy;
  int iy = blockIdx.x * blockDim.x + threadIdx.x;
  int iz = bidx_y * blockDim.y + threadIdx.y;

  if (!(iy < my - 2 * (2-BND) && iz < mz - 2 * (2-BND)))
    return;
  iy -= BND;
  iz -= BND;

  real *d_flds = d_flds0 + p * size;

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

// ----------------------------------------------------------------------
// cuda2_push_mflds_E_yz

void
cuda2_push_mflds_E_yz(struct psc_mfields *mflds)
{
  struct psc_mfields_cuda2 *sub = psc_mfields_cuda2(mflds);

  if (mflds->nr_patches == 0) {
    return;
  }

  struct psc_patch *patch = &ppsc->patch[0];

  real dt = ppsc->dt;
  real cny = .5f * ppsc->dt / patch->dx[1];
  real cnz = .5f * ppsc->dt / patch->dx[2];
  assert(patch->ldims[0] == 1);

  unsigned int size = mflds->nr_fields *
    sub->im[0] * sub->im[1] * sub->im[2];
  int my = sub->im[1];
  int mz = sub->im[2];

  int dimBlock[2] = { BLOCKSIZE_Y, BLOCKSIZE_Z };
  int grid[2]  = { (patch->ldims[1] + 2*BND + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y,
		   (patch->ldims[2] + 2*BND + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z };
  int dimGrid[2] = { grid[0], grid[1] * mflds->nr_patches };

  RUN_KERNEL(dimGrid, dimBlock,
	     push_fields_E_yz, (sub->d_flds, dt, cny, cnz, my, mz,
				size, grid[1]));
}

// ----------------------------------------------------------------------
// cuda2_push_mflds_H_yz

void
cuda2_push_mflds_H_yz(struct psc_mfields *mflds)
{
  struct psc_mfields_cuda2 *sub = psc_mfields_cuda2(mflds);

  if (mflds->nr_patches == 0) {
    return;
  }

  struct psc_patch *patch = &ppsc->patch[0];

  real cny = .5f * ppsc->dt / patch->dx[1];
  real cnz = .5f * ppsc->dt / patch->dx[2];
  assert(patch->ldims[0] == 1);

  unsigned int size = mflds->nr_fields *
    sub->im[0] * sub->im[1] * sub->im[2];
  int my = sub->im[1];
  int mz = sub->im[2];

  int dimBlock[2] = { BLOCKSIZE_Y, BLOCKSIZE_Z };
  int grid[2]  = { (patch->ldims[1] + 2*BND + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y,
		   (patch->ldims[2] + 2*BND + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z };
  int dimGrid[2] = { grid[0], grid[1] * mflds->nr_patches };

  RUN_KERNEL(dimGrid, dimBlock,
	     push_fields_H_yz, (sub->d_flds, cny, cnz, my, mz,
				size, grid[1]));
}

