
#include "cuda_bits.h"
#include "cuda_mfields.h"

#include "fields.hxx"

#include <psc.h>

#define BND 2
#define BLOCKSIZE_X 1
#define BLOCKSIZE_Y 16
#define BLOCKSIZE_Z 16

void cuda_marder_correct_yz_gold(struct cuda_mfields* cmflds,
                                 struct cuda_mfields* cmf, int p, float fac[3],
                                 int ly[3], int ry[3], int lz[3], int rz[3])
{
  auto mflds = hostMirror(*cmflds);
  auto mf = hostMirror(*cmf);

  copy(*cmflds, mflds);
  copy(*cmf, mf);

  auto flds = make_Fields3d<dim_xyz>(mflds[p]);
  auto f = make_Fields3d<dim_xyz>(mf[p]);

  Int3 ldims = cmflds->grid().ldims;
  for (int iz = -1; iz < ldims[2]; iz++) {
    for (int iy = -1; iy < ldims[1]; iy++) {
      if (iy >= -ly[1] && iy < ry[1] && iz >= -ly[2] && iz < ry[2]) {
        flds(EY, 0, iy, iz) += fac[1] * (f(0, 0, iy + 1, iz) - f(0, 0, iy, iz));
      }

      if (iy >= -lz[1] && iy < rz[1] && iz >= -lz[2] && iz < rz[2]) {
        flds(EZ, 0, iy, iz) += fac[2] * (f(0, 0, iy, iz + 1) - f(0, 0, iy, iz));
      }
    }
  }

  copy(mflds, *cmflds);
}

__global__ static void marder_correct_yz(DFields d_flds, DFields d_f,
                                         float facy, float facz, int lyy,
                                         int lyz, int ryy, int ryz, int lzy,
                                         int lzz, int rzy, int rzz, int my,
                                         int mz)
{
  int iy = blockIdx.x * blockDim.x + threadIdx.x;
  int iz = blockIdx.y * blockDim.y + threadIdx.y;

  iy -= BND;
  iz -= BND;

  auto _d_flds = make_Fields3d<dim_xyz>(d_flds);
  auto _d_f = make_Fields3d<dim_xyz>(d_f);

  if (iy >= -lyy && iy < ryy && iz >= -lyz && iz < ryz) {
    _d_flds(EY, 0, iy, iz) +=
      facy * (_d_f(0, 0, iy + 1, iz) - _d_f(0, 0, iy, iz));
  }

  if (iy >= -lzy && iy < rzy && iz >= -lzz && iz < rzz) {
    _d_flds(EZ, 0, iy, iz) +=
      facz * (_d_f(0, 0, iy, iz + 1) - _d_f(0, 0, iy, iz));
  }
}

void cuda_marder_correct_yz(struct cuda_mfields* cmflds,
                            struct cuda_mfields* cmf, int p, float fac[3],
                            int ly[3], int ry[3], int lz[3], int rz[3])
{
#if 0
  cuda_marder_correct_yz_gold(mflds, mf, p, fac, ly, ry, lz, rz);
  return;
#endif

  if (cmflds->n_patches() == 0) {
    return;
  }

  int my = cmflds->im(1);
  int mz = cmflds->im(2);

  int grid[2] = {(my + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y,
                 (mz + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z};
  dim3 dimBlock(BLOCKSIZE_Y, BLOCKSIZE_Z);
  dim3 dimGrid(grid[0], grid[1]);

  marder_correct_yz<<<dimGrid, dimBlock>>>((*cmflds)[p], (*cmf)[p], fac[1],
                                           fac[2], ly[1], ly[2], ry[1], ry[2],
                                           lz[1], lz[2], rz[1], rz[2], my, mz);
  cuda_sync_if_enabled();
}

// ======================================================================

__global__ static void calc_dive_yz(DFields flds, DFields f, float dy, float dz,
                                    int ldimsy, int ldimsz, int my, int mz)
{
  int iy = blockIdx.x * blockDim.x + threadIdx.x;
  int iz = blockIdx.y * blockDim.y + threadIdx.y;

  if (iy >= ldimsy || iz >= ldimsz) {
    return;
  }

  auto _flds = make_Fields3d<dim_xyz>(flds);
  auto _f = make_Fields3d<dim_xyz>(f);
  _f(0, 0, iy, iz) = ((_flds(EY, 0, iy, iz) - _flds(EY, 0, iy - 1, iz)) / dy +
                      (_flds(EZ, 0, iy, iz) - _flds(EZ, 0, iy, iz - 1)) / dz);
}

void cuda_mfields_calc_dive_yz(struct cuda_mfields* cmflds,
                               struct cuda_mfields* cmf, int p)
{
  float dy = cmflds->grid().domain.dx[1];
  float dz = cmflds->grid().domain.dx[2];

  int my = cmflds->im(1);
  int mz = cmflds->im(2);

  int grid[2] = {(my + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y,
                 (mz + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z};
  dim3 dimBlock(BLOCKSIZE_Y, BLOCKSIZE_Z);
  dim3 dimGrid(grid[0], grid[1]);

  calc_dive_yz<<<dimGrid, dimBlock>>>((*cmflds)[p], (*cmf)[p], dy, dz,
                                      cmflds->grid().ldims[1],
                                      cmflds->grid().ldims[2], my, mz);
  cuda_sync_if_enabled();
}
