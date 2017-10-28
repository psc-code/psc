
#include "cuda_bits.h"
#include "cuda_mfields.h"
#include "cuda_mfields_const.h"

#include <psc.h>

// the loops include 2 levels of ghost cells
// they really only need -1:2 and -1:1, respectively (for 1st order)
// but always doing 2:2 seems cheap enough

#define BND 2

// FIXME, merge with F3_DEV{,_YZ}, OPT (precalc offset)

#define X3_DEV_OFF_YZ(fldnr, jy,jz)					\
  ((((fldnr)								\
     *d_cmflds_const.im[2] + ((jz)+2))					\
    *d_cmflds_const.im[1] + ((jy)+2))					\
   *1 + (0))

#undef F3_DEV

#define F3_DEV(fldnr,ix,jy,jz)			\
  (d_flds)[X3_DEV_OFF_YZ(fldnr, jy,jz)]

#define F3_DDEV(d_flds, fldnr,ix,jy,jz)		\
  (d_flds)[X3_DEV_OFF_YZ(fldnr, jy,jz)]

__global__ static void
push_fields_E_yz(float *d_flds0, float dt, float cny, float cnz,
		 unsigned int size, int gridy)
{
  int bidx_y = blockIdx.y % gridy;
  int p = blockIdx.y / gridy;
  int iy = blockIdx.x * blockDim.x + threadIdx.x;
  int iz = bidx_y * blockDim.y + threadIdx.y;

  if (!(iy < d_cmflds_const.im[1] - 2 * (2-BND) &&
	iz < d_cmflds_const.im[2] - 2 * (2-BND)))
    return;
  iy -= BND;
  iz -= BND;

  float *d_flds = d_flds0 + p * size;

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
push_fields_H_yz(float *d_flds0, float cny, float cnz,
		 unsigned int size, int gridy)
{
  int bidx_y = blockIdx.y % gridy;
  int p = blockIdx.y / gridy;
  int iy = blockIdx.x * blockDim.x + threadIdx.x;
  int iz = bidx_y * blockDim.y + threadIdx.y;

  if (!(iy < d_cmflds_const.im[1] - 2 * (2-BND) &&
	iz < d_cmflds_const.im[2] - 2 * (2-BND)))
    return;
  iy -= BND;
  iz -= BND;

  float *d_flds = d_flds0 + p * size;

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
cuda_push_fields_E_yz(struct cuda_mfields *cmflds, float dt)
{
  if (cmflds->n_patches == 0) {
    return;
  }

  cuda_mfields_const_set(cmflds);

  float cny = .5f * dt / cmflds->dx[1];
  float cnz = .5f * dt / cmflds->dx[2];
  assert(cmflds->ldims[0] == 1);

  unsigned int size = cmflds->n_fields * cmflds->n_cells_per_patch;

  int grid[2]  = { (cmflds->ldims[1] + 2*BND + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y,
		   (cmflds->ldims[2] + 2*BND + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z };
  dim3 dimBlock(BLOCKSIZE_Y, BLOCKSIZE_Z);
  dim3 dimGrid(grid[0], grid[1] * cmflds->n_patches);

  push_fields_E_yz<<<dimGrid, dimBlock>>>(cmflds->d_flds, dt, cny, cnz, size, grid[1]);
  cuda_sync_if_enabled();
}

EXTERN_C void
cuda_push_fields_H_yz(struct cuda_mfields *cmflds, float dt)
{
  if (cmflds->n_patches == 0) {
    return;
  }

  cuda_mfields_const_set(cmflds);

  float cny = .5f * dt / cmflds->dx[1];
  float cnz = .5f * dt / cmflds->dx[2];

  unsigned int size = cmflds->n_fields * cmflds->n_cells_per_patch;

  int grid[2]  = { (cmflds->ldims[1] + 2*BND + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y,
		   (cmflds->ldims[2] + 2*BND + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z };
  dim3 dimBlock(BLOCKSIZE_Y, BLOCKSIZE_Z);
  dim3 dimGrid(grid[0], grid[1] * cmflds->n_patches);

  push_fields_H_yz<<<dimGrid, dimBlock>>>(cmflds->d_flds, cny, cnz, size, grid[1]);
  cuda_sync_if_enabled();
}

void
cuda_marder_correct_yz_gold(struct cuda_mfields *cmflds, struct cuda_mfields *cmf,
			    int p, float fac[3],
			    int ly[3], int ry[3],
			    int lz[3], int rz[3])
{
  fields_single_t flds = cuda_mfields_get_host_fields(cmflds);
  fields_single_t f = cuda_mfields_get_host_fields(cmf);
  
  cuda_mfields_copy_from_device(cmflds, p, flds, EX, EX + 3);
  cuda_mfields_copy_from_device(cmf, p, f, 0, 1);
  
  for (int iz = -1; iz < cmflds->ldims[2]; iz++) {
    for (int iy = -1; iy < cmflds->ldims[1]; iy++) {
      if (iy >= -ly[1] && iy < ry[1] &&
	  iz >= -ly[2] && iz < ry[2]) {
	_F3_S(flds, EY, 0,iy,iz) += 
	  fac[1] * (_F3_S(f, 0, 0,iy+1,iz) - _F3_S(f, 0, 0,iy,iz));
	}
      
      if (iy >= -lz[1] && iy < rz[1] &&
	  iz >= -lz[2] && iz < rz[2]) {
	_F3_S(flds, EZ, 0,iy,iz) += 
	  fac[2] * (_F3_S(f, 0, 0,iy,iz+1) - _F3_S(f, 0, 0,iy,iz));
      }
    }
  }
  
  cuda_mfields_copy_to_device(cmflds, p, flds, EX, EX + 3);
  
  fields_single_t_dtor(&flds);
  fields_single_t_dtor(&f);
}

__global__ static void
marder_correct_yz(float *d_flds, float *d_f, float facy, float facz,
		  int lyy, int lyz, int ryy, int ryz,
		  int lzy, int lzz, int rzy, int rzz, int my, int mz)
{
  int iy = blockIdx.x * blockDim.x + threadIdx.x;
  int iz = blockIdx.y * blockDim.y + threadIdx.y;

  iy -= BND;
  iz -= BND;

  if (iy >= -lyy && iy < ryy &&
      iz >= -lyz && iz < ryz) {
    F3_DDEV(d_flds, EY, 0,iy,iz) += 
      facy * (F3_DDEV(d_f, 0, 0,iy+1,iz) - F3_DDEV(d_f, 0, 0,iy,iz));
  }
  
  if (iy >= -lzy && iy < rzy &&
      iz >= -lzz && iz < rzz) {
    F3_DDEV(d_flds, EZ, 0,iy,iz) += 
      facz * (F3_DDEV(d_f, 0, 0,iy,iz+1) - F3_DDEV(d_f, 0, 0,iy,iz));
  }
}

void
cuda_marder_correct_yz(struct cuda_mfields *cmflds, struct cuda_mfields *cmf,
		       int p, float fac[3],
		       int ly[3], int ry[3],
		       int lz[3], int rz[3])
{
#if 0
  cuda_marder_correct_yz_gold(mflds, mf, p, fac, ly, ry, lz, rz);
  return;
#endif

  if (cmflds->n_patches == 0) {
    return;
  }

  unsigned int size = cmflds->n_cells_per_patch;
  int my = cmflds->im[1];
  int mz = cmflds->im[2];

  int grid[2]  = { (cmflds->ldims[1] + 2*BND + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y,
		   (cmflds->ldims[2] + 2*BND + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z };
  dim3 dimBlock(BLOCKSIZE_Y, BLOCKSIZE_Z);
  dim3 dimGrid(grid[0], grid[1]);

  marder_correct_yz<<<dimGrid, dimBlock>>>(cmflds->d_flds + p * size * cmflds->n_fields,
					   cmf->d_flds + p * size * cmf->n_fields, fac[1], fac[2],
					   ly[1], ly[2], ry[1], ry[2],
					   lz[1], lz[2], rz[1], rz[2], my, mz);
  cuda_sync_if_enabled();
}

// ======================================================================

__global__ static void
calc_dive_yz(float *flds, float *f, float dy, float dz,
	     int ldimsy, int ldimsz, int my, int mz)
{
  int iy = blockIdx.x * blockDim.x + threadIdx.x;
  int iz = blockIdx.y * blockDim.y + threadIdx.y;

  if (iy >= ldimsy || iz >= ldimsz) {
    return;
  }

  F3_DDEV(f, 0, 0,iy,iz) = 
    ((F3_DDEV(flds, EY, 0,iy,iz) - F3_DDEV(flds, EY, 0,iy-1,iz)) / dy +
     (F3_DDEV(flds, EZ, 0,iy,iz) - F3_DDEV(flds, EZ, 0,iy,iz-1)) / dz);
}

void
cuda_mfields_calc_dive_yz(struct cuda_mfields *cmflds, struct cuda_mfields *cmf, int p)
{
  float dy = cmflds->dx[1];
  float dz = cmflds->dx[2];

  int my = cmflds->im[1];
  int mz = cmflds->im[2];

  int grid[2]  = { (cmflds->ldims[1] + 2*BND + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y,
		   (cmflds->ldims[2] + 2*BND + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z };
  dim3 dimBlock(BLOCKSIZE_Y, BLOCKSIZE_Z);
  dim3 dimGrid(grid[0], grid[1]);

  calc_dive_yz<<<dimGrid, dimBlock>>>(cmflds->d_flds_by_patch[p], cmf->d_flds_by_patch[p], dy, dz,
				      cmflds->ldims[1], cmflds->ldims[2], my, mz);
  cuda_sync_if_enabled();
}

