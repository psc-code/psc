
#pragma once

#include "fields_item.hxx"
#include "psc_fields_c.h"
#ifdef USE_CUDA
#include "psc_fields_cuda.h"
#endif

// ======================================================================
// Item_jeh
//
// Main fields in their natural staggering

template <typename MfieldsState>
class Item_jeh
{
public:
  using value_type = typename MfieldsState::real_t;
  using space_type = typename MfieldsState::space;

  static std::string name() { return "jeh"; }
  static int n_comps() { return 9; }
  static std::vector<std::string> comp_names()
  {
    return {"jx_ec", "jy_ec", "jz_ec", "ex_ec", "ey_ec",
            "ez_ec", "hx_fc", "hy_fc", "hz_fc"};
  }

  auto operator()(MfieldsState& mflds) const { return mflds.storage().view(); };
};

// ======================================================================
// psc::item::

namespace psc
{
namespace item
{

namespace
{
// FIXME is this namespace appropriate?
auto s0 = _s(1, _);
auto sm = _s(_, -1);

template <typename GT, typename... SLICES>
auto viewFromSlices(GT& gt, gt::gslice slices[], SLICES... moreSlices)
{
  return gt.view(slices[0], slices[1], slices[2], moreSlices...);
}

template <typename E>
Int3 getBnd(const E& flds, const Grid_t& grid)
{
  return {(flds.shape(0) - grid.domain.ldims[0]) / 2,
          (flds.shape(1) - grid.domain.ldims[1]) / 2,
          (flds.shape(2) - grid.domain.ldims[2]) / 2};
}
} // namespace

// ======================================================================
// div_nc

template <typename E>
static auto div_nc(const Grid_t& grid, const E& flds)
{
  Int3 bnd = getBnd(flds, grid);
  auto dxyz = grid.domain.dx;

  auto res = gt::full<gt::expr_value_type<E>, gt::expr_space_type<E>>(
    {grid.ldims[0], grid.ldims[1], grid.ldims[2], 1, grid.n_patches()}, 0);

  // initial values for slices
  gt::gslice trims[3] = {_all, _all, _all};
  gt::gslice lhs[3] = {s0, s0, s0};
  gt::gslice rhs[3][3] = {{sm, s0, s0}, {s0, sm, s0}, {s0, s0, sm}};

  // modify slices to accommodate invariant axes
  for (int a = 0; a < 3; ++a) {
    if (!grid.isInvar(a)) {
      trims[a] = _s(bnd[a] - 1, -bnd[a]);
    } else {
      lhs[a] = _all;
      for (int a2 = 0; a2 < 3; ++a2)
        rhs[a2][a] = _all;
    }
  }

  auto&& trimmed_flds = viewFromSlices(flds, trims);

  for (int a = 0; a < 3; ++a) {
    if (!grid.isInvar(a)) {
      // FIXME use +=
      res.view(_all, _all, _all, 0) =
        res.view(_all, _all, _all, 0) +
        (viewFromSlices(trimmed_flds, lhs, a) -
         viewFromSlices(trimmed_flds, rhs[a], a)) /
          dxyz[a];
    }
  }

  return res;
}

// ----------------------------------------------------------------------
// grad_ec

template <typename E>
static auto grad_ec(const E& fld, const Grid_t& grid)
{
  Int3 bnd = getBnd(fld, grid);
  auto dxyz = grid.domain.dx;

  auto res = gt::full<gt::expr_value_type<E>, gt::expr_space_type<E>>(
    {grid.ldims[0], grid.ldims[1], grid.ldims[2], 3, grid.n_patches()}, 0);

  // initial values for slices
  gt::gslice trims[3] = {_all, _all, _all};
  gt::gslice lhs[3][3] = {{s0, sm, sm}, {sm, s0, sm}, {sm, sm, s0}};
  gt::gslice rhs[3] = {sm, sm, sm};

  // modify slices to accommodate invariant axes
  for (int a = 0; a < 3; ++a) {
    if (!grid.isInvar(a)) {
      trims[a] = _s(bnd[a], 1 - bnd[a]);
    } else {
      rhs[a] = _all;
      for (int a2 = 0; a2 < 3; ++a2)
        lhs[a2][a] = _all;
    }
  }

  auto&& trimmed_flds = viewFromSlices(fld, trims);

  for (int a = 0; a < 3; ++a) {
    if (!grid.isInvar(a)) {
      res.view(_all, _all, _all, a) = (viewFromSlices(trimmed_flds, lhs[a], a) -
                                       viewFromSlices(trimmed_flds, rhs, a)) /
                                      dxyz[a];
    }
  }

  return res;
}

} // namespace item
} // namespace psc

// ======================================================================
// Item_dive

template <typename MfieldsState>
class Item_dive
{
public:
  static char const* name() { return "dive"; }
  static int n_comps() { return 1; }
  static std::vector<std::string> comp_names() { return {"dive"}; }

  auto operator()(MfieldsState& mflds) const
  {
    return psc::item::div_nc(
      mflds.grid(), mflds.storage().view(_all, _all, _all, _s(EX, EX + 3)));
  }
};

// ======================================================================
// Item_divj

template <typename MfieldsState>
class Item_divj
{
public:
  static char const* name() { return "divj"; }
  static int n_comps() { return 1; }
  static std::vector<std::string> comp_names() { return {"divj"}; }

  auto operator()(MfieldsState& mflds) const
  {
    return psc::item::div_nc(
      mflds.grid(), mflds.storage().view(_all, _all, _all, _s(JXI, JXI + 3)));
  }
};

// ======================================================================
// Item_grad

template <typename Mfields>
class Item_grad
{
public:
  using value_type = typename Mfields::real_t;

  static char const* name() { return "grad"; }
  static int n_comps() { return 3; }
  static std::vector<std::string> comp_names()
  {
    return {"gradx", "grady", "gradz"};
  }

  Item_grad(Mfields& mflds) : mflds_{mflds} {}

  const Grid_t& grid() const { return mflds_.grid(); }

  auto gt() const { return psc::item::grad_ec(mflds_.gt(), mflds_.grid()); }

private:
  Mfields& mflds_;
};
