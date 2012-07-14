
#include "psc_cuda.h"

// the loops include 2 levels of ghost cells
// they really only need -1:2 and -1:1, respectively (for 1st order)
// but always doing 2:2 seems cheap enough

#define BND 2

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
push_fields_E_yz(real *d_flds0, real dt, real cny, real cnz, int my, int mz,
		 unsigned int size, int gridy)
{
  int bidx_y = blockIdx.y % gridy;
  int p = blockIdx.y / gridy;
  int iy = blockIdx.x * blockDim.x + threadIdx.x;
  int iz = bidx_y * blockDim.y + threadIdx.y;

  if (!(iy < my - 2 * (3-BND) && iz < mz - 2 * (3-BND)))
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

  if (!(iy < my - 2 * (3-BND) && iz < mz - 2 * (3-BND)))
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

#define BLOCKSIZE_X 1
#define BLOCKSIZE_Y 16
#define BLOCKSIZE_Z 16

EXTERN_C void
cuda_push_fields_E_yz(struct psc_mfields *mflds)
{
  if (mflds->nr_patches == 0) {
    return;
  }

  struct psc_mfields_cuda *mflds_cuda = psc_mfields_cuda(mflds);
  struct psc_patch *patch = &ppsc->patch[0];

  real dt = ppsc->dt;
  real cny = .5f * ppsc->dt / patch->dx[1];
  real cnz = .5f * ppsc->dt / patch->dx[2];
  assert(patch->ldims[0] == 1);

  struct psc_fields *flds = psc_mfields_get_patch(mflds, 0);
  unsigned int size = flds->nr_comp * flds->im[0] * flds->im[1] * flds->im[2];
  int my = flds->im[1];
  int mz = flds->im[2];

  int dimBlock[2] = { BLOCKSIZE_Y, BLOCKSIZE_Z };
  int grid[2]  = { (patch->ldims[1] + 2*BND + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y,
		   (patch->ldims[2] + 2*BND + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z };
  int dimGrid[2] = { grid[0], grid[1] * mflds->nr_patches };

  RUN_KERNEL(dimGrid, dimBlock,
	     push_fields_E_yz, (mflds_cuda->d_flds, dt, cny, cnz, my, mz,
				size, grid[1]));
}

EXTERN_C void
cuda_push_fields_H_yz(struct psc_mfields *mflds)
{
  if (mflds->nr_patches == 0) {
    return;
  }

  struct psc_mfields_cuda *mflds_cuda = psc_mfields_cuda(mflds);
  struct psc_patch *patch = &ppsc->patch[0];

  real cny = .5f * ppsc->dt / patch->dx[1];
  real cnz = .5f * ppsc->dt / patch->dx[2];

  struct psc_fields *flds = psc_mfields_get_patch(mflds, 0);
  unsigned int size = flds->nr_comp * flds->im[0] * flds->im[1] * flds->im[2];
  int my = flds->im[1];
  int mz = flds->im[2];

  int dimBlock[2] = { BLOCKSIZE_Y, BLOCKSIZE_Z };
  int grid[2]  = { (patch->ldims[1] + 2*BND + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y,
		   (patch->ldims[2] + 2*BND + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z };
  int dimGrid[2] = { grid[0], grid[1] * mflds->nr_patches };

  RUN_KERNEL(dimGrid, dimBlock,
	     push_fields_H_yz, (mflds_cuda->d_flds, cny, cnz, my, mz,
				size, grid[1]));
}

