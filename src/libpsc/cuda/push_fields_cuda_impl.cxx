
#include "push_fields_cuda_impl.hxx"

#include "fields.hxx"

// the loops include 2 levels of ghost cells
// they really only need -1:2 and -1:1, respectively (for 1st order)
// but always doing 2:2 seems cheap enough

#define BND 2
#define BLOCKSIZE_X 1
#define BLOCKSIZE_Y 16
#define BLOCKSIZE_Z 16

// OPT: precalc offset

// ======================================================================
// dim_yz

__global__ static void push_fields_E_yz(DMFields dmflds, float dt, float cny,
                                        float cnz, int gridy)
{
  int bidx_y = blockIdx.y % gridy;
  int p = blockIdx.y / gridy;
  int iy = blockIdx.x * blockDim.x + threadIdx.x;
  int iz = bidx_y * blockDim.y + threadIdx.y;

  if (!(iy >= 1 && iy < dmflds.im(1) - 2 * (2 - BND) && iz >= 1 &&
        iz < dmflds.im(2) - 2 * (2 - BND)))
    return;

  iy -= BND;
  iz -= BND;

  auto F = make_Fields3d<dim_xyz>(dmflds[p]);

  F(EX, 0, iy, iz) += cny * (F(HZ, 0, iy, iz) - F(HZ, 0, iy - 1, iz)) -
                      cnz * (F(HY, 0, iy, iz) - F(HY, 0, iy, iz - 1)) -
                      dt * F(JXI, 0, iy, iz);

  F(EY, 0, iy, iz) += cnz * (F(HX, 0, iy, iz) - F(HX, 0, iy, iz - 1)) - 0.f -
                      dt * F(JYI, 0, iy, iz);

  F(EZ, 0, iy, iz) += 0.f - cny * (F(HX, 0, iy, iz) - F(HX, 0, iy - 1, iz)) -
                      dt * F(JZI, 0, iy, iz);
}

__global__ static void push_fields_H_yz(DMFields dmflds, float cny, float cnz,
                                        int gridy)
{
  int bidx_y = blockIdx.y % gridy;
  int p = blockIdx.y / gridy;
  int iy = blockIdx.x * blockDim.x + threadIdx.x;
  int iz = bidx_y * blockDim.y + threadIdx.y;

  if (!(iy < dmflds.im(1) - 2 * (2 - BND) - 1 &&
        iz < dmflds.im(2) - 2 * (2 - BND) - 1))
    return;
  iy -= BND;
  iz -= BND;

  auto F = make_Fields3d<dim_xyz>(dmflds[p]);

  F(HX, 0, iy, iz) -= cny * (F(EZ, 0, iy + 1, iz) - F(EZ, 0, iy, iz)) -
                      cnz * (F(EY, 0, iy, iz + 1) - F(EY, 0, iy, iz));

  F(HY, 0, iy, iz) -= cnz * (F(EX, 0, iy, iz + 1) - F(EX, 0, iy, iz)) - 0.f;

  F(HZ, 0, iy, iz) -= 0.f - cny * (F(EX, 0, iy + 1, iz) - F(EX, 0, iy, iz));
}

void cuda_push_fields_E_yz(struct cuda_mfields* cmflds, float dt)
{
  if (cmflds->n_patches() == 0) {
    return;
  }

  assert(cmflds->n_comps() == NR_FIELDS);

  float cny = dt / cmflds->grid().domain.dx[1];
  float cnz = dt / cmflds->grid().domain.dx[2];
  assert(cmflds->im(0) == 1);
  assert(cmflds->ib(1) == -BND);
  assert(cmflds->ib(2) == -BND);

  int grid[2] = {(cmflds->im(1) + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y,
                 (cmflds->im(2) + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z};
  dim3 dimBlock(BLOCKSIZE_Y, BLOCKSIZE_Z);
  dim3 dimGrid(grid[0], grid[1] * cmflds->n_patches());

  push_fields_E_yz<<<dimGrid, dimBlock>>>(*cmflds, dt, cny, cnz, grid[1]);
  cuda_sync_if_enabled();
}

void cuda_push_fields_H_yz(struct cuda_mfields* cmflds, float dt)
{
  if (cmflds->n_patches() == 0) {
    return;
  }

  float cny = dt / cmflds->grid().domain.dx[1];
  float cnz = dt / cmflds->grid().domain.dx[2];

  int grid[2] = {(cmflds->im(1) + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y,
                 (cmflds->im(2) + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z};
  dim3 dimBlock(BLOCKSIZE_Y, BLOCKSIZE_Z);
  dim3 dimGrid(grid[0], grid[1] * cmflds->n_patches());

  push_fields_H_yz<<<dimGrid, dimBlock>>>(*cmflds, cny, cnz, grid[1]);
  cuda_sync_if_enabled();
}

// ----------------------------------------------------------------------
// push_E

void PushFieldsCuda::push_E(MfieldsStateCuda& mflds, double dt_fac, dim_yz tag)
{
  cuda_push_fields_E_yz(mflds.cmflds(), dt_fac * mflds.grid().dt);
}

void PushFieldsCuda::push_E(MfieldsStateCuda& mflds, double dt_fac, dim_xyz tag)
{
  cuda_push_fields_E_xyz(mflds.cmflds(), dt_fac * mflds.grid().dt);
}

// ----------------------------------------------------------------------
// push_H

void PushFieldsCuda::push_H(MfieldsStateCuda& mflds, double dt_fac, dim_yz tag)
{
  cuda_push_fields_H_yz(mflds.cmflds(), dt_fac * mflds.grid().dt);
}

void PushFieldsCuda::push_H(MfieldsStateCuda& mflds, double dt_fac, dim_xyz tag)
{
  cuda_push_fields_H_xyz(mflds.cmflds(), dt_fac * mflds.grid().dt);
}
