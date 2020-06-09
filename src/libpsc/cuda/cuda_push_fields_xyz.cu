
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
push_fields_E_xyz(DMFields dmflds, float dt, float cnx, float cny, float cnz, int gridz)
{
  int bidx_z = blockIdx.z % gridz;
  int p = blockIdx.z / gridz;
  int ix = blockIdx.x * blockDim.x + threadIdx.x;
  int iy = blockIdx.y * blockDim.y + threadIdx.y;
  int iz = bidx_z     * blockDim.z + threadIdx.z;

  if (!(ix < dmflds.im(0) - 2 * (2-BND) - 1 &&
	iy < dmflds.im(1) - 2 * (2-BND) - 1 &&
	iz < dmflds.im(2) - 2 * (2-BND) - 1))
    return;
  ix -= BND;
  iy -= BND;
  iz -= BND;

  DFields F = dmflds[p];

  F(EX, ix,iy,iz) +=
    cny * (F(HZ, ix,iy,iz) - F(HZ, ix,iy-1,iz)) -
    cnz * (F(HY, ix,iy,iz) - F(HY, ix,iy,iz-1)) -
    dt * F(JXI, ix,iy,iz);
  
  F(EY, ix,iy,iz) +=
    cnz * (F(HX, ix,iy,iz) - F(HX, ix,iy,iz-1)) -
    cnx * (F(HZ, ix,iy,iz) - F(HZ, ix-1,iy,iz)) -
    dt * F(JYI, ix,iy,iz);
  
  F(EZ, ix,iy,iz) +=
    cnx * (F(HY, ix,iy,iz) - F(HY, ix-1,iy,iz)) -
    cny * (F(HX, ix,iy,iz) - F(HX, ix,iy-1,iz)) -
    dt * F(JZI, ix,iy,iz);
}

__global__ static void
push_fields_H_xyz(DMFields dmflds, float cnx, float cny, float cnz, int gridz)
{
  int bidx_z = blockIdx.z % gridz;
  int p = blockIdx.z / gridz;
  int ix = blockIdx.x * blockDim.x + threadIdx.x;
  int iy = blockIdx.y * blockDim.y + threadIdx.y;
  int iz = bidx_z     * blockDim.z + threadIdx.z;

  if (!(ix < dmflds.im(0) - 2 * (2-BND) - 1 &&
	iy < dmflds.im(1) - 2 * (2-BND) - 1 &&
	iz < dmflds.im(2) - 2 * (2-BND) - 1))
    return;
  ix -= BND;
  iy -= BND;
  iz -= BND;

  DFields F = dmflds[p];

  F(HX, ix,iy,iz) -=
    cny * (F(EZ, ix,iy+1,iz) - F(EZ, ix,iy,iz)) -
    cnz * (F(EY, ix,iy,iz+1) - F(EY, ix,iy,iz));
  
  F(HY, ix,iy,iz) -=
    cnz * (F(EX, ix,iy,iz+1) - F(EX, ix,iy,iz)) -
    cnx * (F(EZ, ix+1,iy,iz) - F(EZ, ix,iy,iz));
  
  F(HZ, ix,iy,iz) -=
    cnx * (F(EY, ix+1,iy,iz) - F(EY, ix,iy,iz)) -
    cny * (F(EX, ix,iy+1,iz) - F(EX, ix,iy,iz));
}

#define BLOCKSIZE_X 8
#define BLOCKSIZE_Y 8
#define BLOCKSIZE_Z 8

void
cuda_push_fields_E_xyz(struct cuda_mfields *cmflds, float dt)
{
  if (cmflds->n_patches() == 0) {
    return;
  }

  assert(cmflds->n_comps() == NR_FIELDS);
  assert(cmflds->ib(0) == -BND);
  assert(cmflds->ib(1) == -BND);
  assert(cmflds->ib(2) == -BND);

  float cnx = dt / cmflds->grid().domain.dx[0];
  float cny = dt / cmflds->grid().domain.dx[1];
  float cnz = dt / cmflds->grid().domain.dx[2];

  int grid[3]  = { (cmflds->im(0) + BLOCKSIZE_X - 1) / BLOCKSIZE_X,
		   (cmflds->im(1) + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y,
		   (cmflds->im(2) + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z };
  dim3 dimBlock(BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z);
  dim3 dimGrid(grid[0], grid[1], grid[2] * cmflds->n_patches());


  push_fields_E_xyz<<<dimGrid, dimBlock>>>(*cmflds, dt, cnx, cny, cnz, grid[2]);
  cuda_sync_if_enabled();
}

void
cuda_push_fields_H_xyz(struct cuda_mfields *cmflds, float dt)
{
  if (cmflds->n_patches() == 0) {
    return;
  }

  float cnx = dt / cmflds->grid().domain.dx[0];
  float cny = dt / cmflds->grid().domain.dx[1];
  float cnz = dt / cmflds->grid().domain.dx[2];

  int grid[3]  = { (cmflds->im(0) + BLOCKSIZE_X - 1) / BLOCKSIZE_X,
		   (cmflds->im(1) + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y,
		   (cmflds->im(2) + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z };
  dim3 dimBlock(BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z);
  dim3 dimGrid(grid[0], grid[1], grid[2] * cmflds->n_patches());

  push_fields_H_xyz<<<dimGrid, dimBlock>>>(*cmflds, cnx, cny, cnz, grid[2]);
  cuda_sync_if_enabled();
}

__global__ static void
marder_correct_xyz(DFields d_flds, DFields d_f,
		   float facx, float facy, float facz,
		   int lxx, int lxy, int lxz,
		   int rxx, int rxy, int rxz,
		   int lyx, int lyy, int lyz,
		   int ryx, int ryy, int ryz,
		   int lzx, int lzy, int lzz,
		   int rzx, int rzy, int rzz,
		   int mx, int my, int mz)
{
  int ix = blockIdx.x * blockDim.x + threadIdx.x;
  int iy = blockIdx.y * blockDim.y + threadIdx.y;
  int iz = blockIdx.z * blockDim.z + threadIdx.z;

  ix -= BND;
  iy -= BND;
  iz -= BND;

  if (ix >= -lxx && ix < rxx &&
      iy >= -lxy && iy < rxy &&
      iz >= -lxz && iz < rxz) {
    d_flds(EX, ix,iy,iz) += 
      facx * (d_f(0, ix+1,iy,iz) - d_f(0, ix,iy,iz));
  }
  
  if (ix >= -lyx && ix < ryx &&
      iy >= -lyy && iy < ryy &&
      iz >= -lyz && iz < ryz) {
    d_flds(EY, ix,iy,iz) += 
      facy * (d_f(0, ix,iy+1,iz) - d_f(0, ix,iy,iz));
  }
  
  if (ix >= -lzx && ix < rzx &&
      iy >= -lzy && iy < rzy &&
      iz >= -lzz && iz < rzz) {
    d_flds(EZ, ix,iy,iz) += 
      facz * (d_f(0, ix,iy,iz+1) - d_f(0, ix,iy,iz));
  }
}

void
cuda_marder_correct_xyz(struct cuda_mfields *cmflds, struct cuda_mfields *cmf,
			int p, float fac[3],
			int lx[3], int rx[3],
			int ly[3], int ry[3],
			int lz[3], int rz[3])
{
  if (cmflds->n_patches() == 0) {
    return;
  }

  int mx = cmflds->im(0);
  int my = cmflds->im(1);
  int mz = cmflds->im(2);

  int grid[3]  = { (mx + BLOCKSIZE_X - 1) / BLOCKSIZE_X,
		   (my + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y,
		   (mz + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z };
  dim3 dimBlock(BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z);
  dim3 dimGrid(grid[0], grid[1], grid[2]);

  marder_correct_xyz<<<dimGrid, dimBlock>>>((*cmflds)[p], (*cmf)[p],
					    fac[0], fac[1], fac[2],
					    lx[0], lx[1], lx[2], rx[0], rx[1], rx[2],
					    ly[0], ly[1], ly[2], ry[0], ry[1], ry[2],
					    lz[0], lz[1], lz[2], rz[0], rz[1], rz[2],
					    mx, my, mz);
  cuda_sync_if_enabled();
}

