
#include "cuda_bits.h"
#include "cuda_mfields.h"

#include "fields.hxx"

#include <psc.h>

#define BND 2
#define BLOCKSIZE_X 1
#define BLOCKSIZE_Y 16
#define BLOCKSIZE_Z 16

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
