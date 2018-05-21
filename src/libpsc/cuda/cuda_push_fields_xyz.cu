
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
  if (cmflds->n_patches == 0) {
    return;
  }

  assert(cmflds->n_fields == NR_FIELDS);
  assert(cmflds->ib[0] == -BND);
  assert(cmflds->ib[1] == -BND);
  assert(cmflds->ib[2] == -BND);

  float cnx = dt / cmflds->grid().domain.dx[0];
  float cny = dt / cmflds->grid().domain.dx[1];
  float cnz = dt / cmflds->grid().domain.dx[2];

  int grid[3]  = { (cmflds->im[0] + BLOCKSIZE_X - 1) / BLOCKSIZE_X,
		   (cmflds->im[1] + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y,
		   (cmflds->im[2] + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z };
  dim3 dimBlock(BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z);
  dim3 dimGrid(grid[0], grid[1], grid[2] * cmflds->n_patches);


  push_fields_E_xyz<<<dimGrid, dimBlock>>>(*cmflds, dt, cnx, cny, cnz, grid[2]);
  cuda_sync_if_enabled();
}

void
cuda_push_fields_H_xyz(struct cuda_mfields *cmflds, float dt)
{
  if (cmflds->n_patches == 0) {
    return;
  }

  float cnx = dt / cmflds->grid().domain.dx[0];
  float cny = dt / cmflds->grid().domain.dx[1];
  float cnz = dt / cmflds->grid().domain.dx[2];

  int grid[3]  = { (cmflds->im[0] + BLOCKSIZE_X - 1) / BLOCKSIZE_X,
		   (cmflds->im[1] + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y,
		   (cmflds->im[2] + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z };
  dim3 dimBlock(BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z);
  dim3 dimGrid(grid[0], grid[1], grid[2] * cmflds->n_patches);

  push_fields_H_xyz<<<dimGrid, dimBlock>>>(*cmflds, cnx, cny, cnz, grid[2]);
  cuda_sync_if_enabled();
}

