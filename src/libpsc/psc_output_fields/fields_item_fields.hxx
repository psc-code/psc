
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

namespace psc
{
namespace item
{

template <typename E>
inline auto div_yz(const E& flds, const Grid_t& grid)
{
  Int3 bnd = {(flds.shape(0) - grid.domain.ldims[0]) / 2,
              (flds.shape(1) - grid.domain.ldims[1]) / 2,
              (flds.shape(2) - grid.domain.ldims[2]) / 2};

  auto dxyz = grid.domain.dx;

  auto s0 = _s(1, _);
  auto sm = _s(_, -1);

  auto res = gt::empty<gt::expr_value_type<E>, gt::expr_space_type<E>>(
    {grid.ldims[0], grid.ldims[1], grid.ldims[2], 1, grid.n_patches()});

  auto _flds =
    flds.view(_all, _s(-1 + bnd[1], -bnd[1]), _s(-1 + bnd[2], -bnd[2]));

  res.view(_all, _all, _all, 0) =
    (_flds.view(_all, s0, s0, 1) - _flds.view(_all, sm, s0, 1)) / dxyz[1] +
    (_flds.view(_all, s0, s0, 2) - _flds.view(_all, s0, sm, 2)) / dxyz[2];

  return res;
}

template <typename E>
inline auto div_xyz(const E& flds, const Grid_t& grid)
{
  Int3 bnd = {(flds.shape(0) - grid.domain.ldims[0]) / 2,
              (flds.shape(1) - grid.domain.ldims[1]) / 2,
              (flds.shape(2) - grid.domain.ldims[2]) / 2};
  auto dxyz = grid.domain.dx;

  auto s0 = _s(1, _);
  auto sm = _s(_, -1);

  auto res = gt::empty<gt::expr_value_type<E>, gt::expr_space_type<E>>(
    {grid.ldims[0], grid.ldims[1], grid.ldims[2], 1, grid.n_patches()});

  auto _flds = flds.view(_s(-1 + bnd[0], -bnd[0]), _s(-1 + bnd[1], -bnd[1]),
                         _s(-1 + bnd[2], -bnd[2]));

  res.view(_all, _all, _all, 0) =
    (_flds.view(s0, s0, s0, 0) - _flds.view(sm, s0, s0, 0)) / dxyz[0] +
    (_flds.view(s0, s0, s0, 1) - _flds.view(s0, sm, s0, 1)) / dxyz[1] +
    (_flds.view(s0, s0, s0, 2) - _flds.view(s0, s0, sm, 2)) / dxyz[2];

  return res;
}

template <typename E>
static auto div_nc(const E& flds, const Grid_t& grid)
{
  if (grid.isInvar(0)) {
    return psc::item::div_yz(flds, grid);
  } else {
    return psc::item::div_xyz(flds, grid);
  }
}

template <typename E>
inline auto grad_yz(const E& fld, const Grid_t& grid)
{
  Int3 bnd = {(fld.shape(0) - grid.domain.ldims[0]) / 2,
              (fld.shape(1) - grid.domain.ldims[1]) / 2,
              (fld.shape(2) - grid.domain.ldims[2]) / 2};
  auto dxyz = grid.domain.dx;

  auto s0 = _s(1, _);
  auto sm = _s(_, -1);

  auto res = gt::empty<gt::expr_value_type<E>, gt::expr_space_type<E>>(
    {grid.ldims[0], grid.ldims[1], grid.ldims[2], 3, grid.n_patches()});

  // FIXME? each field has to be shifted in a different way...
  auto _fldy = fld.view(_s(-1 + bnd[0], -bnd[0]), _s(bnd[1], 1-bnd[1]),
                       _s(-1 + bnd[2], -bnd[2]));
  auto _fldz = fld.view(_s(-1 + bnd[0], -bnd[0]), _s(-1 + bnd[1], -bnd[1]),
                       _s(bnd[2], 1-bnd[2]));

  res.view(_all, _all, _all, 1) =
    (_fldy.view(_all, s0, s0, 0) - _fldy.view(_all, sm, s0, 0)) / dxyz[1];
  res.view(_all, _all, _all, 2) =
    (_fldz.view(_all, s0, s0, 0) - _fldz.view(_all, s0, sm, 0)) / dxyz[2];

  return res;
}

template <typename E>
static auto grad_ec(const E& fld, const Grid_t& grid)
{
  if (grid.isInvar(0)) {
    return psc::item::grad_yz(fld, grid);
  } else {
    // FIXME implement the following
    // return psc::item::grad_xyz(fld, grid);
    assert(false);
  }
}

} // namespace item
} // namespace psc

template <typename MfieldsState>
class Item_dive : public MFexpression<Item_dive<MfieldsState>>
{
public:
  using Real = typename MfieldsState::real_t;

  static char const* name() { return "dive"; }
  static int n_comps() { return 1; }
  static std::vector<std::string> comp_names() { return {"dive"}; }

  Item_dive(MfieldsState& mflds) : mflds_{mflds} {}

  const Grid_t& grid() const { return mflds_.grid(); }
  Int3 ibn() const { return {}; }
  int n_patches() const { return grid().n_patches(); }

  auto gt() const
  {
    return psc::item::div_nc(mflds_.gt().view(_all, _all, _all, _s(EX, EX + 3)),
                             mflds_.grid());
  }

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

  const Grid_t& grid() const { return mflds_.grid(); }
  Int3 ibn() const { return {}; }
  int n_patches() const { return grid().n_patches(); }

  auto gt() const
  {
    return psc::item::div_nc(
      mflds_.gt().view(_all, _all, _all, _s(JXI, JXI + 3)), mflds_.grid());
  }

private:
  MfieldsState& mflds_;
};

// ======================================================================
// Item_grad

template <typename Mfields>
class Item_grad : public MFexpression<Item_grad<Mfields>>
{
public:
  using Real = typename Mfields::real_t;

  static char const* name() { return "grad"; }
  static int n_comps() { return 3; }
  static std::vector<std::string> comp_names()
  {
    return {"gradx", "grady", "gradz"};
  }

  Item_grad(Mfields& mflds) : mflds_{mflds} {}

  const Grid_t& grid() const { return mflds_.grid(); }
  Int3 ibn() const { return {}; }
  int n_patches() const { return grid().n_patches(); }

  auto gt() const
  {
    return psc::item::grad_ec(mflds_.gt().view(_all, _all, _all, _all),
                              mflds_.grid());
  }

private:
  Mfields& mflds_;
};

#undef define_dxdydz
