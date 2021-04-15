
#include "../libpsc/psc_output_fields/fields_item_fields.hxx"
#include "fields_item_dive_cuda.hxx"

#define BND 2
#define BLOCKSIZE_X 1
#define BLOCKSIZE_Y 16
#define BLOCKSIZE_Z 16

// ======================================================================

__global__ static void calc_dive_yz(DFields flds, DFields f, float dy, float dz,
                                    int ldimsy, int ldimsz)
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

void cuda_mfields_calc_dive_yz(MfieldsStateCuda& mflds, MfieldsCuda& mf, int p)
{
  auto cmflds = mflds.cmflds();
  auto cmf = mf.cmflds();
  auto dx = mflds.grid().domain.dx;

  int my = mflds.gt().shape(1);
  int mz = mflds.gt().shape(2);

  int grid[2] = {(my + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y,
                 (mz + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z};
  dim3 dimBlock(BLOCKSIZE_Y, BLOCKSIZE_Z);
  dim3 dimGrid(grid[0], grid[1]);

  calc_dive_yz<<<dimGrid, dimBlock>>>((*cmflds)[p], (*cmf)[p], dx[1], dx[2],
                                      cmflds->grid().ldims[1],
                                      cmflds->grid().ldims[2]);
  cuda_sync_if_enabled();
}
