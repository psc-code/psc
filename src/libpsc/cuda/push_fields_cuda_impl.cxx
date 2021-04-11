
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
