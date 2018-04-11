
#include "cuda_bits.h"
#include "cuda_mfields.h"

#include "fields.hxx"

#include <psc.h>

// the loops include 2 levels of ghost cells
// they really only need -1:2 and -1:1, respectively (for 1st order)
// but always doing 2:2 seems cheap enough

#define BND 2

// OPT: precalc offset

__global__ static void
push_fields_E_yz(DMFields dmflds, float dt, float cny, float cnz, int gridy)
{
  int bidx_y = blockIdx.y % gridy;
  int p = blockIdx.y / gridy;
  int iy = blockIdx.x * blockDim.x + threadIdx.x;
  int iz = bidx_y * blockDim.y + threadIdx.y;

  if (!(iy >= 1 && iy < dmflds.im(1) - 2 * (2-BND) &&
	iz >= 1 && iz < dmflds.im(2) - 2 * (2-BND)))
    return;
  
  iy -= BND;
  iz -= BND;

  DFields F = dmflds[p];

  F(EX, 0,iy,iz) +=
    cny * (F(HZ, 0,iy,iz) - F(HZ, 0,iy-1,iz)) -
    cnz * (F(HY, 0,iy,iz) - F(HY, 0,iy,iz-1)) -
    dt * F(JXI, 0,iy,iz);
  
  F(EY, 0,iy,iz) +=
    cnz * (F(HX, 0,iy,iz) - F(HX, 0,iy,iz-1)) -
    0.f -
    dt * F(JYI, 0,iy,iz);
  
  F(EZ, 0,iy,iz) +=
    0.f -
    cny * (F(HX, 0,iy,iz) - F(HX, 0,iy-1,iz)) -
    dt * F(JZI, 0,iy,iz);
}

__global__ static void
push_fields_H_yz(DMFields dmflds, float cny, float cnz, int gridy)
{
  int bidx_y = blockIdx.y % gridy;
  int p = blockIdx.y / gridy;
  int iy = blockIdx.x * blockDim.x + threadIdx.x;
  int iz = bidx_y * blockDim.y + threadIdx.y;

  if (!(iy < dmflds.im(1) - 2 * (2-BND) - 1 &&
	iz < dmflds.im(2) - 2 * (2-BND) - 1))
    return;
  iy -= BND;
  iz -= BND;

  DFields F = dmflds[p];

  F(HX, 0,iy,iz) -=
    cny * (F(EZ, 0,iy+1,iz) - F(EZ, 0,iy,iz)) -
    cnz * (F(EY, 0,iy,iz+1) - F(EY, 0,iy,iz));
  
  F(HY, 0,iy,iz) -=
    cnz * (F(EX, 0,iy,iz+1) - F(EX, 0,iy,iz)) -
    0.f;
  
  F(HZ, 0,iy,iz) -=
    0.f -
    cny * (F(EX, 0,iy+1,iz) - F(EX, 0,iy,iz));
}

#define BLOCKSIZE_X 1
#define BLOCKSIZE_Y 16
#define BLOCKSIZE_Z 16

void
cuda_push_fields_E_yz(struct cuda_mfields *cmflds, float dt)
{
  if (cmflds->n_patches == 0) {
    return;
  }

  assert(cmflds->n_fields == NR_FIELDS);

  float cny = dt / cmflds->grid().domain.dx[1];
  float cnz = dt / cmflds->grid().domain.dx[2];
  assert(cmflds->im[0] == 1);
  assert(cmflds->ib[1] == -BND);
  assert(cmflds->ib[2] == -BND);

  int grid[2]  = { (cmflds->im[1] + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y,
		   (cmflds->im[2] + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z };
  dim3 dimBlock(BLOCKSIZE_Y, BLOCKSIZE_Z);
  dim3 dimGrid(grid[0], grid[1] * cmflds->n_patches);

  push_fields_E_yz<<<dimGrid, dimBlock>>>(*cmflds, dt, cny, cnz, grid[1]);
  cuda_sync_if_enabled();
}

void
cuda_push_fields_H_yz(struct cuda_mfields *cmflds, float dt)
{
  if (cmflds->n_patches == 0) {
    return;
  }

  float cny = dt / cmflds->grid().domain.dx[1];
  float cnz = dt / cmflds->grid().domain.dx[2];

  int grid[2]  = { (cmflds->im[1] + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y,
		   (cmflds->im[2] + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z };
  dim3 dimBlock(BLOCKSIZE_Y, BLOCKSIZE_Z);
  dim3 dimGrid(grid[0], grid[1] * cmflds->n_patches);

  push_fields_H_yz<<<dimGrid, dimBlock>>>(*cmflds, cny, cnz, grid[1]);
  cuda_sync_if_enabled();
}

void
cuda_marder_correct_yz_gold(struct cuda_mfields *cmflds, struct cuda_mfields *cmf,
			    int p, float fac[3],
			    int ly[3], int ry[3],
			    int lz[3], int rz[3])
{
  fields_single_t flds = cmflds->get_host_fields();
  Fields3d<fields_single_t> Flds(flds);
  fields_single_t f = cmf->get_host_fields();
  Fields3d<fields_single_t> F(f);
  
  cmflds->copy_from_device(p, flds, EX, EX + 3);
  cmf->copy_from_device(p, f, 0, 1);

  Int3 ldims = cmflds->grid().ldims;
  for (int iz = -1; iz < ldims[2]; iz++) {
    for (int iy = -1; iy < ldims[1]; iy++) {
      if (iy >= -ly[1] && iy < ry[1] &&
	  iz >= -ly[2] && iz < ry[2]) {
	Flds(EY, 0,iy,iz) += fac[1] * (F(0, 0,iy+1,iz) - F(0, 0,iy,iz));
	}
      
      if (iy >= -lz[1] && iy < rz[1] &&
	  iz >= -lz[2] && iz < rz[2]) {
	Flds(EZ, 0,iy,iz) += fac[2] * (F(0, 0,iy,iz+1) - F(0, 0,iy,iz));
      }
    }
  }
  
  cmflds->copy_to_device(p, flds, EX, EX + 3);

  flds.dtor();
  f.dtor();
}

__global__ static void
marder_correct_yz(DFields d_flds, DFields d_f, float facy, float facz,
		  int lyy, int lyz, int ryy, int ryz,
		  int lzy, int lzz, int rzy, int rzz, int my, int mz)
{
  int iy = blockIdx.x * blockDim.x + threadIdx.x;
  int iz = blockIdx.y * blockDim.y + threadIdx.y;

  iy -= BND;
  iz -= BND;

  if (iy >= -lyy && iy < ryy &&
      iz >= -lyz && iz < ryz) {
    d_flds(EY, 0,iy,iz) += 
      facy * (d_f(0, 0,iy+1,iz) - d_f(0, 0,iy,iz));
  }
  
  if (iy >= -lzy && iy < rzy &&
      iz >= -lzz && iz < rzz) {
    d_flds(EZ, 0,iy,iz) += 
      facz * (d_f(0, 0,iy,iz+1) - d_f(0, 0,iy,iz));
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

  int my = cmflds->im[1];
  int mz = cmflds->im[2];

  int grid[2]  = { (my + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y,
		   (mz + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z };
  dim3 dimBlock(BLOCKSIZE_Y, BLOCKSIZE_Z);
  dim3 dimGrid(grid[0], grid[1]);

  marder_correct_yz<<<dimGrid, dimBlock>>>((*cmflds)[p], (*cmf)[p],
					   fac[1], fac[2],
					   ly[1], ly[2], ry[1], ry[2],
					   lz[1], lz[2], rz[1], rz[2], my, mz);
  cuda_sync_if_enabled();
}

// ======================================================================

__global__ static void
calc_dive_yz(DFields flds, DFields f, float dy, float dz,
	     int ldimsy, int ldimsz, int my, int mz)
{
  int iy = blockIdx.x * blockDim.x + threadIdx.x;
  int iz = blockIdx.y * blockDim.y + threadIdx.y;

  if (iy >= ldimsy || iz >= ldimsz) {
    return;
  }

  f(0, 0,iy,iz) = 
    ((flds(EY, 0,iy,iz) - flds(EY, 0,iy-1,iz)) / dy +
     (flds(EZ, 0,iy,iz) - flds(EZ, 0,iy,iz-1)) / dz);
}

void
cuda_mfields_calc_dive_yz(struct cuda_mfields *cmflds, struct cuda_mfields *cmf, int p)
{
  float dy = cmflds->grid().domain.dx[1];
  float dz = cmflds->grid().domain.dx[2];

  int my = cmflds->im[1];
  int mz = cmflds->im[2];

  int grid[2]  = { (my + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y,
		   (mz + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z };
  dim3 dimBlock(BLOCKSIZE_Y, BLOCKSIZE_Z);
  dim3 dimGrid(grid[0], grid[1]);

  calc_dive_yz<<<dimGrid, dimBlock>>>((*cmflds)[p], (*cmf)[p], dy, dz,
				      cmflds->grid().ldims[1], cmflds->grid().ldims[2], my, mz);
  cuda_sync_if_enabled();
}

