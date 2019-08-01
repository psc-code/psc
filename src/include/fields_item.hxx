
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

#define POFI_MAX_COMPS (16)

using fld_names_t = std::array<const char*, POFI_MAX_COMPS>;

// ======================================================================
// FieldsItemBase

struct FieldsItemBase
{
  virtual ~FieldsItemBase(){};

  virtual void run(const Grid_t& grid, MfieldsStateBase& mflds_base,
                   MparticlesBase& mprts_base)
  {
    assert(0);
  }

  virtual MfieldsBase& mres() = 0;

  virtual const char* name() const = 0;

  virtual int n_comps(const Grid_t& grid) const = 0;

  virtual std::vector<std::string> comp_names() = 0;

  bool inited = true; // FIXME hack to avoid dtor call when not yet constructed
};

// ======================================================================
// FieldsItemFields

template <typename Item>
struct FieldsItemFields : FieldsItemBase
{
  using Mfields = typename Item::Mfields;

  const char* name() const override { return Item::name; }

  int n_comps(const Grid_t& grid) const override { return Item::n_comps; }

  FieldsItemFields(const Grid_t& grid) : mres_{grid, Item::n_comps, grid.ibn} {}

  template <typename MfieldsState>
  void operator()(const Grid_t& grid, MfieldsState& mflds)
  {
    Item::run(grid, mflds, mres_);
  }

  void run(const Grid_t& grid, MfieldsStateBase& mflds_base,
           MparticlesBase& mprts_base) override
  {
    using MfieldsState = typename Item::MfieldsState;
    auto& mflds = mflds_base.get_as<MfieldsState>(0, mflds_base._n_comps());
    (*this)(grid, mflds);
    mflds_base.put_as(mflds, 0, 0);
  }

  virtual std::vector<std::string> comp_names() override
  {
    std::vector<std::string> comp_names;
    for (int m = 0; m < Item::n_comps; m++) {
      comp_names.emplace_back(Item::fld_names()[m]);
    }
    return comp_names;
  }

  virtual MfieldsBase& mres() override { return mres_; }

  Mfields& result() { return mres_; }

private:
  Mfields mres_;
};

template <template <typename> class Item>
struct _FieldsItemFields
{
  using Mfields = MfieldsC;

  _FieldsItemFields(const Grid_t& grid) : mres_{grid, Item<MfieldsFake>::n_comps, grid.ibn}
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

  using MfieldsFake = MfieldsC;

  static const char* name() { return Item<MfieldsFake>::name; }
  static int n_comps(const Grid_t& grid) { return Item<MfieldsFake>::n_comps; }

  static std::vector<std::string> comp_names()
  {
    std::vector<std::string> comp_names;
    for (int m = 0; m < Item<MfieldsFake>::n_comps; m++) {
      comp_names.emplace_back(Item<MfieldsFake>::fld_names()[m]);
    }
    return comp_names;
  }

  MfieldsBase& mres() { return mres_; }
  Mfields& result() { return mres_; }

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
// ItemMomentCRTP
//
// deriving from this class adds the result field mres_

template <typename Derived, typename MF>
struct ItemMomentCRTP
{
  using Mfields = MF;

  ItemMomentCRTP(const Grid_t& grid)
    : mres_{grid,
            int(Derived::n_comps *
                ((Derived::flags & POFI_BY_KIND) ? grid.kinds.size() : 1)),
            grid.ibn}
  {
    auto n_comps = Derived::n_comps;
    auto fld_names = Derived::fld_names();
    auto& kinds = grid.kinds;
    assert(n_comps <= POFI_MAX_COMPS);

    if (!(Derived::flags & POFI_BY_KIND)) {
      for (int m = 0; m < n_comps; m++) {
        comp_names_.emplace_back(fld_names[m]);
      }
    } else {
      for (int k = 0; k < kinds.size(); k++) {
        for (int m = 0; m < n_comps; m++) {
          comp_names_.emplace_back(std::string(fld_names[m]) + "_" +
                                   kinds[k].name);
        }
      }
    }
  }

  Mfields& result() { return mres_; }
  std::vector<std::string> comp_names() { return comp_names_; }

protected:
  Mfields mres_;
  std::vector<std::string> comp_names_;
};

// ======================================================================
// ItemMomentAddBnd

template <typename Moment_t, typename Bnd = Bnd_<typename Moment_t::Mfields>>
struct ItemMomentAddBnd
  : ItemMomentCRTP<ItemMomentAddBnd<Moment_t>, typename Moment_t::Mfields>
{
  using Base =
    ItemMomentCRTP<ItemMomentAddBnd<Moment_t>, typename Moment_t::Mfields>;
  using Mfields = typename Moment_t::Mfields;

  constexpr static const char* name = Moment_t::name;
  constexpr static int n_comps = Moment_t::n_comps;
  constexpr static fld_names_t fld_names() { return Moment_t::fld_names(); }
  constexpr static int flags = Moment_t::flags;

  ItemMomentAddBnd(const Grid_t& grid) : Base{grid}, bnd_{grid, grid.ibn} {}

  template <typename Mparticles>
  void run(Mparticles& mprts)
  {
    auto& mres = this->mres_;
    mres.zero();
    Moment_t::run(mres, mprts);
    for (int p = 0; p < mprts.n_patches(); p++) {
      add_ghosts_boundary(mprts.grid(), mres[p], p, 0, mres.n_comps());
    }

    bnd_.add_ghosts(mres, 0, mres.n_comps());
  }

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
  Bnd bnd_;
};

// ======================================================================
// FieldsItemMoment

template <typename Moment_t>
struct FieldsItemMoment : FieldsItemBase
{
  const char* name() const override { return Moment_t::name; }

  int n_comps(const Grid_t& grid) const override
  {
    return Moment_t::n_comps *
           ((Moment_t::flags & POFI_BY_KIND) ? grid.kinds.size() : 1);
  }

  FieldsItemMoment(const Grid_t& grid) : moment_(grid) {}

  template <typename MfieldsState, typename Mparticles>
  void operator()(const Grid_t& grid, MfieldsState& mflds, Mparticles& mprts)
  {
    moment_.run(mprts);
  }

  virtual MfieldsBase& mres() override { return moment_.result(); }

  virtual std::vector<std::string> comp_names() override
  {
    return moment_.comp_names();
  }

private:
  Moment_t moment_;
};
