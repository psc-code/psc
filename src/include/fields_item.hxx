
#pragma once

#include "particles.hxx"
#include "fields3d.hxx"
#include "fields.hxx"
#include "bnd.hxx"
#include "psc_fields_c.h"
#include "../libpsc/psc_bnd/psc_bnd_impl.hxx"

#include <mrc_profile.h>

#include <map>
#include <string>
#include <array>

enum {
  POFI_BY_KIND    = 2, // this item needs to be replicated by kind
};

#define POFI_MAX_COMPS (16)

using fld_names_t = std::array<const char*, POFI_MAX_COMPS>;

// ======================================================================
// FieldsItemBase

struct FieldsItemBase
{
  virtual ~FieldsItemBase() {};

  virtual void run(MfieldsStateBase& mflds_base, MparticlesBase& mprts_base) = 0;

  virtual MfieldsBase& mres() = 0;

  virtual const char* _name() const = 0;

  virtual int n_comps(const Grid_t& grid) const = 0;

  virtual std::vector<std::string> comp_names() = 0;

  bool inited = true; // FIXME hack to avoid dtor call when not yet constructed
};

// ======================================================================
// FieldsItemFactory

class FieldsItemFactory
{
  using CreateMethod = std::unique_ptr<FieldsItemBase>(*)(const Grid_t& grid);
  using Map = std::map<std::string, CreateMethod>;
  
public:
  FieldsItemFactory() = delete;

  static std::unique_ptr<FieldsItemBase> create(const char *type, const Grid_t& grid)
  {
    auto it = map().find(type);    
    if (it != map().end()) 
        return it->second(grid); // call the create method

    mprintf("unknown type %s\n", type);
    assert(0);
  }

  static void registerType(const char* name, CreateMethod create)
  {
    map()[name] = create;
  }

private:
  static Map& map()
  {
    static std::map<std::string, CreateMethod> s_map;
    return s_map;
  }
  
};

void registerFieldsItemFields();

// ======================================================================
// FieldsItemOps

template<typename Item>
struct FieldsItemOps
{
  FieldsItemOps() {
    FieldsItemFactory::registerType(Item::name(), create);
  }

  static std::unique_ptr<FieldsItemBase> create(const Grid_t& grid)
  {
    return std::unique_ptr<FieldsItemBase>{new Item{grid}};
  }
};

// ======================================================================
// FieldsItemFields

template<typename Item>
struct FieldsItemFields : FieldsItemBase
{
  using MfieldsState = typename Item::MfieldsState;
  using Mfields = typename Item::Mfields;
  
  static const char* name()
  {
    if (std::is_same<Mfields, MfieldsC>::value && strcmp(Item::name, "dive") != 0) {
      return Item::name;
    } else {
      return strdup((std::string{Item::name} + "_" + Mfields_traits<Mfields>::name).c_str());
    }
  }

  const char* _name() const override { return name(); }

  int n_comps(const Grid_t& grid) const override { return Item::n_comps; }
 
  FieldsItemFields(const Grid_t& grid)
    : mres_{grid, Item::n_comps, grid.ibn}
  {}

  void operator()(MfieldsState& mflds)
  {
    Item::run(mflds, mres_);
  }

  void run(MfieldsStateBase& mflds_base, MparticlesBase& mprts_base) override
  {
    auto& mflds = mflds_base.get_as<MfieldsState>(0, mflds_base.n_comps());
    (*this)(mflds);
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

// ======================================================================
// FieldsItemFieldsOps
//
// Adapter from per-patch Item with ::set

template<typename ItemPatch>
struct ItemLoopPatches : ItemPatch
{
  using MfieldsState = typename ItemPatch::MfieldsState;
  using Mfields = typename ItemPatch::Mfields;
  using fields_t = typename Mfields::fields_t;
  using Fields = Fields3d<fields_t>;
  
  static void run(MfieldsState& mflds, Mfields& mres)
  {
    auto& grid = mres.grid();
    
    for (int p = 0; p < mres.n_patches(); p++) {
      Fields F(mflds[p]), R(mres[p]);
      grid.Foreach_3d(0, 0, [&](int i, int j, int k) {
	  ItemPatch::set(grid, R, F, i,j,k);
	});
    }
  }
};

template<typename Item_t>
using FieldsItemFieldsOps = FieldsItemOps<FieldsItemFields<ItemLoopPatches<Item_t>>>;

// ======================================================================
// FieldsItemMomentOps

// ----------------------------------------------------------------------
// ItemMomentCRTP
//
// deriving from this class adds the result field mres_

template<typename Derived, typename MF>
struct ItemMomentCRTP
{
  using Mfields = MF;
  
  ItemMomentCRTP(const Grid_t& grid, MPI_Comm comm)
    : mres_{grid, int(Derived::n_comps * ((Derived::flags & POFI_BY_KIND) ? grid.kinds.size() : 1)), grid.ibn}
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
	  comp_names_.emplace_back(std::string(fld_names[m]) + "_" + kinds[k].name);
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

// ----------------------------------------------------------------------
// ItemMomentAddBnd

template<typename Moment_t, typename Bnd = Bnd_<typename Moment_t::Mfields>>
struct ItemMomentAddBnd : ItemMomentCRTP<ItemMomentAddBnd<Moment_t>, typename Moment_t::Mfields>
{
  using Base = ItemMomentCRTP<ItemMomentAddBnd<Moment_t>, typename Moment_t::Mfields>;
  using Mfields = typename Moment_t::Mfields;
  using Mparticles = typename Moment_t::Mparticles;
  using fields_t = typename Mfields::fields_t;
  using Fields = Fields3d<fields_t>;

  constexpr static const char* name = Moment_t::name;
  constexpr static int n_comps = Moment_t::n_comps;
  constexpr static fld_names_t fld_names() { return Moment_t::fld_names(); }
  constexpr static int flags = Moment_t::flags;

  ItemMomentAddBnd(const Grid_t& grid, MPI_Comm comm)
    : Base{grid, comm},
      bnd_{grid, grid.ibn}
  {}

  void run(Mparticles& mprts)
  {
    auto& mres = this->mres_;
    mres.zero();
    Moment_t::run(mres, mprts);
    for (int p = 0; p < mprts.n_patches(); p++) {
      add_ghosts_boundary(mres[p], p, 0, mres.n_comps());
    }

    bnd_.add_ghosts(mres, 0, mres.n_comps());
  }

  // ----------------------------------------------------------------------
  // boundary stuff FIXME, should go elsewhere...

  static void add_ghosts_reflecting_lo(fields_t flds, int p, int d, int mb, int me)
  {
    Fields F(flds);
    const int *ldims = flds.grid().ldims;

    int bx = ldims[0] == 1 ? 0 : 1;
    if (d == 1) {
      for (int iz = -1; iz < ldims[2] + 1; iz++) {
	for (int ix = -bx; ix < ldims[0] + bx; ix++) {
	  int iy = 0; {
	    for (int m = mb; m < me; m++) {
	      F(m, ix,iy,iz) += F(m, ix,iy-1,iz);
	    }
	  }
	}
      }
    } else if (d == 2) {
      for (int iy = 0*-1; iy < ldims[1] + 0*1; iy++) {
	for (int ix = -bx; ix < ldims[0] + bx; ix++) {
	  int iz = 0; {
	    for (int m = mb; m < me; m++) {
	      F(m, ix,iy,iz) += F(m, ix,iy,iz-1);
	    }
	  }
	}
      }
    } else {
      assert(0);
    }
  }

  static void add_ghosts_reflecting_hi(fields_t flds, int p, int d, int mb, int me)
  {
    Fields F(flds);
    const int *ldims = flds.grid().ldims;

    int bx = ldims[0] == 1 ? 0 : 1;
    if (d == 1) {
      for (int iz = -1; iz < ldims[2] + 1; iz++) {
	for (int ix = -bx; ix < ldims[0] + bx; ix++) {
	  int iy = ldims[1] - 1; {
	    for (int m = mb; m < me; m++) {
	      F(m, ix,iy,iz) += F(m, ix,iy+1,iz);
	    }
	  }
	}
      }
    } else if (d == 2) {
      for (int iy = 0*-1; iy < ldims[1] + 0*1; iy++) {
	for (int ix = -bx; ix < ldims[0] + bx; ix++) {
	  int iz = ldims[2] - 1; {
	    for (int m = mb; m < me; m++) {
	      F(m, ix,iy,iz) += F(m, ix,iy,iz+1);
	    }
	  }
	}
      }
    } else {
      assert(0);
    }
  }

  static void add_ghosts_boundary(fields_t res, int p, int mb, int me)
  {
    const auto& grid = res.grid();
    // lo
    for (int d = 0; d < 3; d++) {
      if (grid.atBoundaryLo(p, d)) {
	if (grid.bc.prt_lo[d] == BND_PRT_REFLECTING ||
	    grid.bc.prt_lo[d] == BND_PRT_OPEN) {
	  add_ghosts_reflecting_lo(res, p, d, mb, me);
	}
      }
    }
    // hi
    for (int d = 0; d < 3; d++) {
      if (grid.atBoundaryHi(p, d)) {
	if (grid.bc.prt_hi[d] == BND_PRT_REFLECTING ||
	    grid.bc.prt_hi[d] == BND_PRT_OPEN) {
	  add_ghosts_reflecting_hi(res, p, d, mb, me);
	}
      }
    }
  }

private:
  Bnd bnd_;
};

// ----------------------------------------------------------------------
// FieldsItemMoment

template<typename Moment_t>
struct FieldsItemMoment : FieldsItemBase
{
  using Mparticles = typename Moment_t::Mparticles;
  
  static const char* name()
  {
    return strdup((std::string(Moment_t::name) + "_" +
		   Mparticles_traits<Mparticles>::name).c_str());
  }

  const char* _name() const override { return name(); }

  int n_comps(const Grid_t& grid) const override
  {
    return Moment_t::n_comps * ((Moment_t::flags & POFI_BY_KIND) ? grid.kinds.size() : 1);
  }
  
  FieldsItemMoment(const Grid_t& grid)
    : moment_(grid, grid.comm())
  {}

  void run(MfieldsStateBase& mflds_base, MparticlesBase& mprts_base) override
  {
    auto& mprts = mprts_base.get_as<Mparticles>();
    moment_.run(mprts);
    mprts_base.put_as(mprts, MP_DONT_COPY);
  }

  virtual MfieldsBase& mres() override { return moment_.result(); }

  virtual std::vector<std::string> comp_names()  override { return moment_.comp_names(); }
  
private:
  Moment_t moment_;
};

// ----------------------------------------------------------------------
// FieldsItemMomentOps
  
template<typename Moment_t>
using FieldsItemMomentOps = FieldsItemOps<FieldsItemMoment<ItemMomentAddBnd<Moment_t>>>;

