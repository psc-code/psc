
#include "push_fields_cuda_impl.hxx"

#include "fields.hxx"

// OPT: precalc offset

// ======================================================================
// dim_yz

#if 0
template <typename E>
GT_INLINE static void push_fields_E_yz(E gt, float dt, float cny, float cnz,
                                       int iy, int iz, int p)
{
  if (!(iy >= 1 && iz >= 1))
    return;

  gt(0, iy, iz, EX, p) =
    gt(0, iy, iz, EX, p) +
    cny * (gt(0, iy, iz, HZ, p) - gt(0, iy - 1, iz, HZ, p)) -
    cnz * (gt(0, iy, iz, HY, p) - gt(0, iy, iz - 1, HY, p)) -
    dt * gt(0, iy, iz, JXI, p);

  gt(0, iy, iz, EY, p) =
    gt(0, iy, iz, EY, p) +
    cnz * (gt(0, iy, iz, HX, p) - gt(0, iy, iz - 1, HX, p)) - 0.f -
    dt * gt(0, iy, iz, JYI, p);

  gt(0, iy, iz, EZ, p) =
    gt(0, iy, iz, EZ, p) + 0.f -
    cny * (gt(0, iy, iz, HX, p) - gt(0, iy - 1, iz, HX, p)) -
    dt * gt(0, iy, iz, JZI, p);
}
#endif

void PushFieldsCuda::push_E(MfieldsStateCuda& mflds, double dt_fac, dim_yz tag)
{
  if (mflds.n_patches() == 0) {
    return;
  }

  assert(mflds.n_comps() == NR_FIELDS);
  assert(mflds.ibn() == Int3({0, 2, 2}));

  double dt = dt_fac * mflds.grid().dt;
  float cny = dt / mflds.grid().domain.dx[1];
  float cnz = dt / mflds.grid().domain.dx[2];

#if 0
  auto shape = mflds.gt().shape();
  auto gt = mflds.gt().to_kernel();
  gt::launch<3>({shape[1], shape[2], mflds.n_patches()},
                GT_LAMBDA(int iy, int iz, int p) {
                  push_fields_E_yz(gt, dt, cny, cnz, iy, iz, p);
                });
#endif

  auto gt = mflds.gt();

  gt.view(0, _s(1, _), _s(1, _), EX) =
    gt.view(0, _s(1, _), _s(1, _), EX) +
    (cny * (gt.view(0, _s(1, _), _s(1, _), HZ) -
            gt.view(0, _s(_, -1), _s(1, _), HZ)) -
     cnz * (gt.view(0, _s(1, _), _s(1, _), HY) -
            gt.view(0, _s(1, _), _s(_, -1), HY)) -
     0.f - dt * gt.view(0, _s(1, _), _s(1, _), JXI));

  gt.view(0, _s(1, _), _s(1, _), EY) =
    gt.view(0, _s(1, _), _s(1, _), EY) +
    (cnz * (gt.view(0, _s(1, _), _s(1, _), HX) -
            gt.view(0, _s(1, _), _s(_, -1), HX)) -
     0.f - dt * gt.view(0, _s(1, _), _s(1, _), JYI));

  gt.view(0, _s(1, _), _s(1, _), EZ) =
    gt.view(0, _s(1, _), _s(1, _), EZ) +
    (0.f -
     cny * (gt.view(0, _s(1, _), _s(1, _), HX) -
            gt.view(0, _s(_, -1), _s(1, _), HX)) -
     dt * gt.view(0, _s(1, _), _s(1, _), JZI));

  cuda_sync_if_enabled();
}

void PushFieldsCuda::push_H(MfieldsStateCuda& mflds, double dt_fac, dim_yz tag)
{
  if (mflds.n_patches() == 0) {
    return;
  }

  double dt = dt_fac * mflds.grid().dt;
  float cny = dt / mflds.grid().domain.dx[1];
  float cnz = dt / mflds.grid().domain.dx[2];

  auto gt = mflds.gt();

  gt.view(0, _s(_, -1), _s(_, -1), HX) =
    gt.view(0, _s(_, -1), _s(_, -1), HX) -
    (cny * (gt.view(0, _s(1, _), _s(_, -1), EZ) -
            gt.view(0, _s(_, -1), _s(_, -1), EZ)) -
     cnz * (gt.view(0, _s(_, -1), _s(1, _), EY) -
            gt.view(0, _s(_, -1), _s(_, -1), EY)));

  gt.view(0, _s(_, -1), _s(_, -1), HY) =
    gt.view(0, _s(_, -1), _s(_, -1), HY) -
    (cnz * (gt.view(0, _s(_, -1), _s(1, _), EX) -
            gt.view(0, _s(_, -1), _s(_, -1), EX)));

  gt.view(0, _s(_, -1), _s(_, -1), HZ) =
    gt.view(0, _s(_, -1), _s(_, -1), HZ) -
    (0.f - cny * (gt.view(0, _s(1, _), _s(_, -1), EX) -
                  gt.view(0, _s(_, -1), _s(_, -1), EX)));

  cuda_sync_if_enabled();
}

// ======================================================================
// dim_xyz

#define BND 2
#define BLOCKSIZE_X 8
#define BLOCKSIZE_Y 8
#define BLOCKSIZE_Z 8

__global__ static void push_fields_E_xyz(DMFields dmflds, float dt, float cnx,
                                         float cny, float cnz, int gridz)
{
  int bidx_z = blockIdx.z % gridz;
  int p = blockIdx.z / gridz;
  int ix = blockIdx.x * blockDim.x + threadIdx.x;
  int iy = blockIdx.y * blockDim.y + threadIdx.y;
  int iz = bidx_z * blockDim.z + threadIdx.z;

  if (!((ix >= 1 && ix < dmflds.im(0) - 2 * (2 - BND)) &&
        (iy >= 1 && iy < dmflds.im(1) - 2 * (2 - BND)) &&
        (iz >= 1 && iz < dmflds.im(2) - 2 * (2 - BND))))
    return;
  ix -= BND;
  iy -= BND;
  iz -= BND;

  auto F = make_Fields3d<dim_xyz>(dmflds[p]);

  F(EX, ix, iy, iz) += cny * (F(HZ, ix, iy, iz) - F(HZ, ix, iy - 1, iz)) -
                       cnz * (F(HY, ix, iy, iz) - F(HY, ix, iy, iz - 1)) -
                       dt * F(JXI, ix, iy, iz);

  F(EY, ix, iy, iz) += cnz * (F(HX, ix, iy, iz) - F(HX, ix, iy, iz - 1)) -
                       cnx * (F(HZ, ix, iy, iz) - F(HZ, ix - 1, iy, iz)) -
                       dt * F(JYI, ix, iy, iz);

  F(EZ, ix, iy, iz) += cnx * (F(HY, ix, iy, iz) - F(HY, ix - 1, iy, iz)) -
                       cny * (F(HX, ix, iy, iz) - F(HX, ix, iy - 1, iz)) -
                       dt * F(JZI, ix, iy, iz);
}

__global__ static void push_fields_H_xyz(DMFields dmflds, float cnx, float cny,
                                         float cnz, int gridz)
{
  int bidx_z = blockIdx.z % gridz;
  int p = blockIdx.z / gridz;
  int ix = blockIdx.x * blockDim.x + threadIdx.x;
  int iy = blockIdx.y * blockDim.y + threadIdx.y;
  int iz = bidx_z * blockDim.z + threadIdx.z;

  if (!(ix < dmflds.im(0) - 2 * (2 - BND) - 1 &&
        iy < dmflds.im(1) - 2 * (2 - BND) - 1 &&
        iz < dmflds.im(2) - 2 * (2 - BND) - 1))
    return;
  ix -= BND;
  iy -= BND;
  iz -= BND;

  auto F = make_Fields3d<dim_xyz>(dmflds[p]);

  F(HX, ix, iy, iz) -= cny * (F(EZ, ix, iy + 1, iz) - F(EZ, ix, iy, iz)) -
                       cnz * (F(EY, ix, iy, iz + 1) - F(EY, ix, iy, iz));

  F(HY, ix, iy, iz) -= cnz * (F(EX, ix, iy, iz + 1) - F(EX, ix, iy, iz)) -
                       cnx * (F(EZ, ix + 1, iy, iz) - F(EZ, ix, iy, iz));

  F(HZ, ix, iy, iz) -= cnx * (F(EY, ix + 1, iy, iz) - F(EY, ix, iy, iz)) -
                       cny * (F(EX, ix, iy + 1, iz) - F(EX, ix, iy, iz));
}

void cuda_push_fields_E_xyz(struct cuda_mfields* cmflds, float dt)
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

  int grid[3] = {(cmflds->im(0) + BLOCKSIZE_X - 1) / BLOCKSIZE_X,
                 (cmflds->im(1) + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y,
                 (cmflds->im(2) + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z};
  dim3 dimBlock(BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z);
  dim3 dimGrid(grid[0], grid[1], grid[2] * cmflds->n_patches());

  push_fields_E_xyz<<<dimGrid, dimBlock>>>(*cmflds, dt, cnx, cny, cnz, grid[2]);
  cuda_sync_if_enabled();
}

void cuda_push_fields_H_xyz(struct cuda_mfields* cmflds, float dt)
{
  if (cmflds->n_patches() == 0) {
    return;
  }

  float cnx = dt / cmflds->grid().domain.dx[0];
  float cny = dt / cmflds->grid().domain.dx[1];
  float cnz = dt / cmflds->grid().domain.dx[2];

  int grid[3] = {(cmflds->im(0) + BLOCKSIZE_X - 1) / BLOCKSIZE_X,
                 (cmflds->im(1) + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y,
                 (cmflds->im(2) + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z};
  dim3 dimBlock(BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z);
  dim3 dimGrid(grid[0], grid[1], grid[2] * cmflds->n_patches());

  push_fields_H_xyz<<<dimGrid, dimBlock>>>(*cmflds, cnx, cny, cnz, grid[2]);
  cuda_sync_if_enabled();
}

// ----------------------------------------------------------------------
// push_E

void PushFieldsCuda::push_E(MfieldsStateCuda& mflds, double dt_fac, dim_xyz tag)
{
  cuda_push_fields_E_xyz(mflds.cmflds(), dt_fac * mflds.grid().dt);
}

// ----------------------------------------------------------------------
// push_H

void PushFieldsCuda::push_H(MfieldsStateCuda& mflds, double dt_fac, dim_xyz tag)
{
  cuda_push_fields_H_xyz(mflds.cmflds(), dt_fac * mflds.grid().dt);
}
