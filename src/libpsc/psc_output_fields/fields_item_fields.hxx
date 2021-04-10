
#pragma once

#include "fields_item.hxx"
#include "psc_fields_c.h"
#ifdef USE_CUDA
#include "psc_fields_cuda.h"
#endif

template <typename E>
inline auto to_gt(const E& e)
{
  auto& grid = e.grid();
  assert(e.ibn() == Int3{});
  auto res = gt::empty<typename E::Real>(
    {grid.ldims[0], grid.ldims[1], grid.ldims[2], e.n_comps(), e.n_patches()});
  auto k_res = res.to_kernel();

  gt::launch<5, gt::space::host>(res.shape(),
                                 [=](int i, int j, int k, int m, int p) {
                                   k_res(i, j, k, m, p) = e(m, {i, j, k}, p);
                                 });

  return res;
}

// ======================================================================

using MfieldsState_t = MfieldsStateDouble;
using Mfields_t = MfieldsC;

// ======================================================================

#define define_dxdydz(dx, dy, dz)                                              \
  int dx _mrc_unused = (grid.isInvar(0)) ? 0 : 1;                              \
  int dy _mrc_unused = (grid.isInvar(1)) ? 0 : 1;                              \
  int dz _mrc_unused = (grid.isInvar(2)) ? 0 : 1

// ======================================================================
// Item_jeh
//
// Main fields in their natural staggering

template <typename MfieldsState>
class Item_jeh : public MFexpression<Item_jeh<MfieldsState>>
{
public:
  using Real = typename MfieldsState::real_t;
  using value_type = Real;
  using space = typename MfieldsState::space;

  static char const* name() { return "jeh"; }
  static int n_comps() { return 9; }
  static std::vector<std::string> comp_names()
  {
    return {"jx_ec", "jy_ec", "jz_ec", "ex_ec", "ey_ec",
            "ez_ec", "hx_fc", "hy_fc", "hz_fc"};
  }

  Item_jeh(MfieldsState& mflds) : mflds_{mflds} {}

  auto gt() const
  {
    auto bnd = mflds_.ibn();
    return mflds_.gt().view(_s(bnd[0], -bnd[0]), _s(bnd[1], -bnd[1]),
                            _s(bnd[2], -bnd[2]));
  }

  const Grid_t& grid() const { return mflds_.grid(); }
  Int3 ibn() const { return {}; }
  int n_patches() const { return grid().n_patches(); }

  MfieldsState& result() const { return mflds_; }

private:
  MfieldsState& mflds_;
};

// ======================================================================
// Item_dive

template <typename MfieldsState>
class Item_dive : public MFexpression<Item_dive<MfieldsState>>
{
public:
  using Real = typename MfieldsState::real_t;

  static char const* name() { return "dive"; }
  static int n_comps() { return 1; }
  static std::vector<std::string> comp_names() { return {"dive"}; }

  Item_dive(MfieldsState& mflds) : mflds_{mflds} {}

  Real operator()(int m, Int3 ijk, int p) const
  {
    const auto& grid = mflds_.grid();
    define_dxdydz(dx, dy, dz);

    auto flds = make_Fields3d<dim_xyz>(mflds_[p]);
    return ((flds(EX, ijk[0], ijk[1], ijk[2]) -
             flds(EX, ijk[0] - dx, ijk[1], ijk[2])) /
              grid.domain.dx[0] +
            (flds(EY, ijk[0], ijk[1], ijk[2]) -
             flds(EY, ijk[0], ijk[1] - dy, ijk[2])) /
              grid.domain.dx[1] +
            (flds(EZ, ijk[0], ijk[1], ijk[2]) -
             flds(EZ, ijk[0], ijk[1], ijk[2] - dz)) /
              grid.domain.dx[2]);
  }

  auto gt() const { return to_gt(*this); }

  const Grid_t& grid() const { return mflds_.grid(); }
  Int3 ibn() const { return {}; }
  int n_patches() const { return grid().n_patches(); }

private:
  MfieldsState& mflds_;
};

// ======================================================================
// Item_divj

// FIXME, almost same as dive

template <typename MfieldsState>
class Item_divj : public MFexpression<Item_divj<MfieldsState>>
{
public:
  using Real = typename MfieldsState::real_t;

  static char const* name() { return "divj"; }
  static int n_comps() { return 1; }
  static std::vector<std::string> comp_names() { return {"divj"}; }

  Item_divj(MfieldsState& mflds) : mflds_{mflds} {}

  Real operator()(int m, Int3 ijk, int p) const
  {
    const auto& grid = mflds_.grid();
    define_dxdydz(dx, dy, dz);
    auto flds = make_Fields3d<dim_xyz>(mflds_[p]);

    return ((flds(JXI, ijk[0], ijk[1], ijk[2]) -
             flds(JXI, ijk[0] - dx, ijk[1], ijk[2])) /
              grid.domain.dx[0] +
            (flds(JYI, ijk[0], ijk[1], ijk[2]) -
             flds(JYI, ijk[0], ijk[1] - dy, ijk[2])) /
              grid.domain.dx[1] +
            (flds(JZI, ijk[0], ijk[1], ijk[2]) -
             flds(JZI, ijk[0], ijk[1], ijk[2] - dz)) /
              grid.domain.dx[2]);
  }

  const Grid_t& grid() const { return mflds_.grid(); }
  Int3 ibn() const { return {}; }
  int n_patches() const { return grid().n_patches(); }

private:
  MfieldsState& mflds_;
};

#undef define_dxdydz
