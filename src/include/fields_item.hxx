
#pragma once

#include "../libpsc/psc_bnd/psc_bnd_impl.hxx"
#include "bnd.hxx"
#include "fields.hxx"
#include "fields3d.hxx"
#include "particles.hxx"
#include "psc_fields_c.h"

#include <mrc_profile.h>

#include <array>
#include <map>
#include <string>

enum
{
  POFI_BY_KIND = 2, // this item needs to be replicated by kind
};

// ======================================================================
// FieldsItemFields

template <typename Item>
struct FieldsItemFields
{
  using Mfields = typename Item::Mfields;

  FieldsItemFields(const Grid_t& grid) : mres_{grid, Item::n_comps, grid.ibn} {}

  template <typename MfieldsState>
  void operator()(const Grid_t& grid, MfieldsState& mflds)
  {
    Item::run(grid, mflds, mres_);
  }

  Mfields& result() { return mres_; }

  static const char* name() { return Item::name; }
  static int n_comps(const Grid_t& grid) { return Item::n_comps; }

  std::vector<std::string> comp_names(const Grid_t& grid)
  {
    return Item::fld_names();
  }

private:
  Mfields mres_;
};

template <template <typename> class Item>
struct _FieldsItemFields
{
  using Mfields = MfieldsC;

  _FieldsItemFields(const Grid_t& grid)
    : mres_{grid, Item<MfieldsFake>::n_comps, grid.ibn}
  {}

  template <typename MfieldsState>
  void operator()(const Grid_t& grid, MfieldsState& mflds)
  {
    Item<MfieldsState> item{mflds};

    for (int p = 0; p < mres_.n_patches(); p++) {
      auto R = mres_[p];
      for (int m = 0; m < item.n_comps; m++) {
        mres_.Foreach_3d(0, 0, [&](int i, int j, int k) {
          R(m, i, j, k) = item(m, {i, j, k}, p);
        });
      }
    }
  }

  Mfields& result() { return mres_; }

  using MfieldsFake = MfieldsC;

  static const char* name() { return Item<MfieldsFake>::name; }
  static int n_comps(const Grid_t& grid) { return Item<MfieldsFake>::n_comps; }

  static std::vector<std::string> comp_names(const Grid_t& grid)
  {
    return Item<MfieldsFake>::fld_names();
  }

private:
  Mfields mres_;
};

// ======================================================================
// ItemLoopPatches
//
// Adapter from per-patch Item with ::set

template <typename ItemPatch>
struct ItemLoopPatches : ItemPatch
{
  using MfieldsState = typename ItemPatch::MfieldsState;
  using Mfields = typename ItemPatch::Mfields;

  static void run(const Grid_t& grid, MfieldsState& mflds, Mfields& mres)
  {
    for (int p = 0; p < mres.n_patches(); p++) {
      auto F = mflds[p];
      auto R = mres[p];
      mres.Foreach_3d(0, 0, [&](int i, int j, int k) {
        ItemPatch::set(grid, R, F, i, j, k);
      });
    }
  }
};

// ======================================================================
// addKindSuffix

inline std::vector<std::string> addKindSuffix(
  const std::vector<std::string>& names, const Grid_t::Kinds& kinds)
{
  std::vector<std::string> result;
  for (int k = 0; k < kinds.size(); k++) {
    for (int m = 0; m < names.size(); m++) {
      result.emplace_back(names[m] + "_" + kinds[k].name);
    }
  }
  return result;
}

// ======================================================================
// ItemMomentCRTP
//
// deriving from this class adds the result field mres_

template <typename Derived, typename MF>
struct ItemMomentCRTP
{
  using Mfields = MF;

  ItemMomentCRTP(const Grid_t& grid)
    : mres_{grid, Derived::n_comps(grid), grid.ibn}
  {}

  Mfields& result() { return mres_; }

protected:
  Mfields mres_;
};

// ======================================================================
// ItemMomentAddBnd

template <typename Moment_t, typename Bnd = Bnd_<typename Moment_t::Mfields>>
struct ItemMomentAddBnd
{
  using Mfields = typename Moment_t::Mfields;

  static const char* name() { return Moment_t::name; }

  static int n_comps(const Grid_t& grid) { return Moment_t::n_comps(grid); }
  static std::vector<std::string> comp_names(const Grid_t& grid)
  {
    return Moment_t::comp_names(grid);
  }

  ItemMomentAddBnd(const Grid_t& grid)
    : moment_{grid},
      bnd_{grid, grid.ibn}
  {}

  template <typename Mparticles>
  void operator()(Mparticles& mprts)
  {
    auto& mres = moment_.result();
    mres.zero();
    Moment_t::run(mres, mprts);
    for (int p = 0; p < mprts.n_patches(); p++) {
      add_ghosts_boundary(mprts.grid(), mres[p], p, 0, mres.n_comps());
    }

    bnd_.add_ghosts(mres, 0, mres.n_comps());
  }

  Mfields& result() { return moment_.result(); }

  // ----------------------------------------------------------------------
  // boundary stuff FIXME, should go elsewhere...

  template <typename FE>
  void add_ghosts_reflecting_lo(const Grid_t& grid, FE flds, int p, int d,
                                int mb, int me)
  {
    auto ldims = grid.ldims;

    int bx = ldims[0] == 1 ? 0 : 1;
    if (d == 1) {
      for (int iz = -1; iz < ldims[2] + 1; iz++) {
        for (int ix = -bx; ix < ldims[0] + bx; ix++) {
          int iy = 0;
          {
            for (int m = mb; m < me; m++) {
              flds(m, ix, iy, iz) += flds(m, ix, iy - 1, iz);
            }
          }
        }
      }
    } else if (d == 2) {
      for (int iy = 0 * -1; iy < ldims[1] + 0 * 1; iy++) {
        for (int ix = -bx; ix < ldims[0] + bx; ix++) {
          int iz = 0;
          {
            for (int m = mb; m < me; m++) {
              flds(m, ix, iy, iz) += flds(m, ix, iy, iz - 1);
            }
          }
        }
      }
    } else {
      assert(0);
    }
  }

  template <typename FE>
  void add_ghosts_reflecting_hi(const Grid_t& grid, FE flds, int p, int d,
                                int mb, int me)
  {
    auto ldims = grid.ldims;

    int bx = ldims[0] == 1 ? 0 : 1;
    if (d == 1) {
      for (int iz = -1; iz < ldims[2] + 1; iz++) {
        for (int ix = -bx; ix < ldims[0] + bx; ix++) {
          int iy = ldims[1] - 1;
          {
            for (int m = mb; m < me; m++) {
              flds(m, ix, iy, iz) += flds(m, ix, iy + 1, iz);
            }
          }
        }
      }
    } else if (d == 2) {
      for (int iy = 0 * -1; iy < ldims[1] + 0 * 1; iy++) {
        for (int ix = -bx; ix < ldims[0] + bx; ix++) {
          int iz = ldims[2] - 1;
          {
            for (int m = mb; m < me; m++) {
              flds(m, ix, iy, iz) += flds(m, ix, iy, iz + 1);
            }
          }
        }
      }
    } else {
      assert(0);
    }
  }

  template <typename FE>
  void add_ghosts_boundary(const Grid_t& grid, FE res, int p, int mb, int me)
  {
    // lo
    for (int d = 0; d < 3; d++) {
      if (grid.atBoundaryLo(p, d)) {
        if (grid.bc.prt_lo[d] == BND_PRT_REFLECTING ||
            grid.bc.prt_lo[d] == BND_PRT_OPEN) {
          add_ghosts_reflecting_lo(grid, res, p, d, mb, me);
        }
      }
    }
    // hi
    for (int d = 0; d < 3; d++) {
      if (grid.atBoundaryHi(p, d)) {
        if (grid.bc.prt_hi[d] == BND_PRT_REFLECTING ||
            grid.bc.prt_hi[d] == BND_PRT_OPEN) {
          add_ghosts_reflecting_hi(grid, res, p, d, mb, me);
        }
      }
    }
  }

private:
  Moment_t moment_;
  Bnd bnd_;
};
