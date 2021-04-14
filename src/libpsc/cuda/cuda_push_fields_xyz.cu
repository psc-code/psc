
#include "cuda_bits.h"
#include "cuda_mfields.h"

#include "fields.hxx"

#include <psc.h>

#define BND 2
#define BLOCKSIZE_X 8
#define BLOCKSIZE_Y 8
#define BLOCKSIZE_Z 8

__global__ static void marder_correct_xyz(
  DFields d_flds, DFields d_f, float facx, float facy, float facz, int lxx,
  int lxy, int lxz, int rxx, int rxy, int rxz, int lyx, int lyy, int lyz,
  int ryx, int ryy, int ryz, int lzx, int lzy, int lzz, int rzx, int rzy,
  int rzz, int mx, int my, int mz)
{
  int ix = blockIdx.x * blockDim.x + threadIdx.x;
  int iy = blockIdx.y * blockDim.y + threadIdx.y;
  int iz = blockIdx.z * blockDim.z + threadIdx.z;

  ix -= BND;
  iy -= BND;
  iz -= BND;

  auto _d_flds = make_Fields3d<dim_xyz>(d_flds);
  auto _d_f = make_Fields3d<dim_xyz>(d_f);

  if (ix >= -lxx && ix < rxx && iy >= -lxy && iy < rxy && iz >= -lxz &&
      iz < rxz) {
    _d_flds(EX, ix, iy, iz) +=
      facx * (_d_f(0, ix + 1, iy, iz) - _d_f(0, ix, iy, iz));
  }

  if (ix >= -lyx && ix < ryx && iy >= -lyy && iy < ryy && iz >= -lyz &&
      iz < ryz) {
    _d_flds(EY, ix, iy, iz) +=
      facy * (_d_f(0, ix, iy + 1, iz) - _d_f(0, ix, iy, iz));
  }

  if (ix >= -lzx && ix < rzx && iy >= -lzy && iy < rzy && iz >= -lzz &&
      iz < rzz) {
    _d_flds(EZ, ix, iy, iz) +=
      facz * (_d_f(0, ix, iy, iz + 1) - _d_f(0, ix, iy, iz));
  }
}

void cuda_marder_correct_xyz(struct cuda_mfields* cmflds,
                             struct cuda_mfields* cmf, int p, float fac[3],
                             int lx[3], int rx[3], int ly[3], int ry[3],
                             int lz[3], int rz[3])
{
  if (cmflds->n_patches() == 0) {
    return;
  }

  int mx = cmflds->im(0);
  int my = cmflds->im(1);
  int mz = cmflds->im(2);

  int grid[3] = {(mx + BLOCKSIZE_X - 1) / BLOCKSIZE_X,
                 (my + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y,
                 (mz + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z};
  dim3 dimBlock(BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z);
  dim3 dimGrid(grid[0], grid[1], grid[2]);

  marder_correct_xyz<<<dimGrid, dimBlock>>>(
    (*cmflds)[p], (*cmf)[p], fac[0], fac[1], fac[2], lx[0], lx[1], lx[2], rx[0],
    rx[1], rx[2], ly[0], ly[1], ly[2], ry[0], ry[1], ry[2], lz[0], lz[1], lz[2],
    rz[0], rz[1], rz[2], mx, my, mz);
  cuda_sync_if_enabled();
}
