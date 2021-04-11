
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

void PushFieldsCuda::push_E(MfieldsStateCuda& mflds, double dt_fac, dim_yz tag)
{
  if (mflds.n_patches() == 0) {
    return;
  }

  assert(mflds.n_comps() == NR_FIELDS);
  assert(mflds.ibn() == Int3({0, BND, BND}));

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

  auto cmflds = mflds.cmflds();
  double dt = dt_fac * mflds.grid().dt;

  float cny = dt / mflds.grid().domain.dx[1];
  float cnz = dt / mflds.grid().domain.dx[2];

  auto shape = mflds.gt().shape();
  int grid[2] = {(shape[1] + BLOCKSIZE_Y - 1) / BLOCKSIZE_Y,
                 (shape[2] + BLOCKSIZE_Z - 1) / BLOCKSIZE_Z};
  dim3 dimBlock(BLOCKSIZE_Y, BLOCKSIZE_Z);
  dim3 dimGrid(grid[0], grid[1] * mflds.n_patches());

  push_fields_H_yz<<<dimGrid, dimBlock>>>(*cmflds, cny, cnz, grid[1]);
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
