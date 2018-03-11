
#pragma once

#include "psc_output_fields_item_private.h"

#include "particles.hxx"
#include "fields3d.hxx"
#include "fields.hxx"
#include "bnd.hxx"
#include "psc_fields_c.h"

#include <mrc_profile.h>

#include <string>

// ======================================================================
// FieldsItemBase

struct FieldsItemBase
{
  virtual void run(PscMfieldsBase mflds_base, PscMparticlesBase mprts_base) = 0;

  virtual psc_mfields* mres() = 0;

  bool inited = true; // FIXME hack to avoid dtor call when not yet constructed
};

// ======================================================================
// PscFieldsItem

template<typename S>
struct PscFieldsItem
{
  using sub_t = S;
  
  static_assert(std::is_convertible<sub_t*, FieldsItemBase*>::value,
  		"sub classes used in PscFieldsItemParticles must derive from FieldsItemBase");
  
  explicit PscFieldsItem(psc_output_fields_item *item)
    : item_(item)
  {}

  void operator()(PscMfieldsBase mflds, PscMparticlesBase mprts, PscMfieldsBase mres)
  {
    // FIXME, we shouldn't even take the mres arg
    assert(!mres.mflds());
    sub()->run(mflds, mprts);
  }

  sub_t* sub() { return mrc_to_subobj(item_, sub_t); }
  sub_t* operator->() { return sub(); }

private:
  psc_output_fields_item* item_;
};

using PscFieldsItemBase = PscFieldsItem<FieldsItemBase>;

// ======================================================================
// PscFieldsItemWrapper

template<typename FieldsItem>
class PscFieldsItemWrapper
{
public:
  const static size_t size = sizeof(FieldsItem);
  
  static void setup(psc_output_fields_item* _item)
  {
    PscFieldsItem<FieldsItem> item(_item);
    new(item.sub()) FieldsItem(psc_output_fields_item_comm(_item),
			       PscBndBase{_item->bnd});
  }

  static void destroy(psc_output_fields_item* _item)
  {
    PscFieldsItem<FieldsItem> item(_item);

    if (!item->inited) return; // ctor hadn't run yet FIXME

    item->~FieldsItem();
  }
};

// ======================================================================
// FieldsItemOps

template<typename Item_t>
struct FieldsItemOps : psc_output_fields_item_ops {
  using Wrapper_t = PscFieldsItemWrapper<Item_t>;
  FieldsItemOps() {
    name      = Item_t::name();
    size      = Wrapper_t::size;
    setup     = Wrapper_t::setup;
    destroy   = Wrapper_t::destroy;
  }
};

// ======================================================================
// FieldsItemFields

template<typename Item>
struct FieldsItemFields : FieldsItemBase
{
  using mfields_t = typename Item::mfields_t;
  
  static char const* name()
  {
    if (std::is_same<mfields_t, PscMfieldsC>::value && strcmp(Item::name, "dive") != 0) {
      return Item::name;
    } else {
      return strdup((std::string{Item::name} + "_" + fields_traits<typename mfields_t::fields_t>::name).c_str());
    }
  }
 
  FieldsItemFields(MPI_Comm comm, PscBndBase bnd)
  {
    mres_base_ = mfields_t::create(comm, ppsc->grid(), Item::n_comps).mflds();
    for (int m = 0; m < Item::n_comps; m++) {
      psc_mfields_set_comp_name(mres_base_, m, Item::fld_names()[m]);
    }
    psc_mfields_list_add(&psc_mfields_base_list, &mres_base_);
  }

  ~FieldsItemFields()
  {
    psc_mfields_list_del(&psc_mfields_base_list, &mres_base_);
    psc_mfields_destroy(mres_base_);
  }

  void run(PscMfieldsBase mflds_base, PscMparticlesBase mprts_base) override
  {
    mfields_t mflds = mflds_base.get_as<mfields_t>(JXI, JXI + 3);
    mfields_t mres(mres_base_);
    Item::run(mflds, mres);
    mflds.put_as(mflds_base, 0, 0);
  }

  virtual psc_mfields* mres() override { return mres_base_; }

private:  
  psc_mfields* mres_base_;
};

// ======================================================================
// FieldsItemFieldsOps
//
// Adapter from per-patch Item with ::set

template<typename ItemPatch>
struct ItemLoopPatches : ItemPatch
{
  using mfields_t = typename ItemPatch::mfields_t;
  using fields_t = typename mfields_t::fields_t;
  using Fields = Fields3d<fields_t>;
  
  static void run(mfields_t mflds, mfields_t mres)
  {
    for (int p = 0; p < mres->n_patches(); p++) {
      Fields F(mflds[p]), R(mres[p]);
      psc_foreach_3d(ppsc, p, i,j,k, 0, 0) {
	ItemPatch::set(R, F, i,j,k);
      } foreach_3d_end;
    }
  }
};

template<typename Item_t>
using FieldsItemFieldsOps = FieldsItemOps<FieldsItemFields<ItemLoopPatches<Item_t>>>;

// ======================================================================
// FieldsItemMomentOps

// ----------------------------------------------------------------------
// ItemMomentLoopPatches

template<typename Moment_t>
struct ItemMomentLoopPatches
{
  using mfields_t = typename Moment_t::mfields_t;
  using mparticles_t = typename Moment_t::mparticles_t;
  using fields_t = typename mfields_t::fields_t;
  using Fields = Fields3d<fields_t>;

  constexpr static const char* name = Moment_t::name;
  constexpr static int n_comps = Moment_t::n_comps;
  constexpr static fld_names_t fld_names() { return Moment_t::fld_names(); }
  constexpr static int flags = Moment_t::flags;

  ItemMomentLoopPatches(PscBndBase bnd)
    : bnd_(bnd)
  {}
  
  void run(mfields_t mres, mparticles_t mprts)
  {
    for (int p = 0; p < mprts->n_patches(); p++) {
      mres[p].zero();
      Moment_t::run(mres[p], mprts[p]);
      add_ghosts_boundary(mres[p], p, 0, mres->n_comps());
    }

    bnd_.add_ghosts(mres.mflds(), 0, mres->n_comps());
  }

  // ----------------------------------------------------------------------
  // boundary stuff FIXME, should go elsewhere...

  static void add_ghosts_reflecting_lo(fields_t flds, int p, int d, int mb, int me)
  {
    Fields F(flds);
    const int *ldims = ppsc->grid().ldims;

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
    const int *ldims = ppsc->grid().ldims;

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
    // lo
    for (int d = 0; d < 3; d++) {
      if (psc_at_boundary_lo(ppsc, p, d)) {
	if (ppsc->domain.bnd_part_lo[d] == BND_PART_REFLECTING ||
	    ppsc->domain.bnd_part_lo[d] == BND_PART_OPEN) {
	  add_ghosts_reflecting_lo(res, p, d, mb, me);
	}
      }
    }
    // hi
    for (int d = 0; d < 3; d++) {
      if (psc_at_boundary_hi(ppsc, p, d)) {
	if (ppsc->domain.bnd_part_hi[d] == BND_PART_REFLECTING ||
	    ppsc->domain.bnd_part_hi[d] == BND_PART_OPEN) {
	  add_ghosts_reflecting_hi(res, p, d, mb, me);
	}
      }
    }
  }

private:
  PscBndBase bnd_;
};

// ----------------------------------------------------------------------
// ItemMoment

template<typename Moment_t>
struct ItemMoment : FieldsItemBase
{
  using mparticles_t = typename Moment_t::mparticles_t;
  using mfields_t = typename Moment_t::mfields_t;
  using fields_t = typename mfields_t::fields_t;
  using Fields = Fields3d<fields_t>;
  
  static const char* name()
  {
    return strdup((std::string(Moment_t::name) + "_" +
		   mparticles_traits<mparticles_t>::name).c_str());
  }

  ItemMoment(MPI_Comm comm, PscBndBase bnd)
    : moment_(bnd)
  {
    auto n_comps = Moment_t::n_comps;
    auto fld_names = Moment_t::fld_names();
    assert(n_comps <= POFI_MAX_COMPS);

    if (!Moment_t::flags & POFI_BY_KIND) {
      mres_base_ = mfields_t::create(comm, ppsc->grid(), n_comps).mflds();
      for (int m = 0; m < n_comps; m++) {
	psc_mfields_set_comp_name(mres_base_, m, fld_names[m]);
      }
    } else {
      mres_base_ = mfields_t::create(comm, ppsc->grid(), n_comps * ppsc->nr_kinds).mflds();
      for (int k = 0; k < ppsc->nr_kinds; k++) {
	for (int m = 0; m < n_comps; m++) {
	  auto s = std::string(fld_names[m]) + "_" + ppsc->kinds[k].name;
	  psc_mfields_set_comp_name(mres_base_, k * n_comps + m, s.c_str());
	}
      }
    }
    psc_mfields_list_add(&psc_mfields_base_list, &mres_base_);
  }

  ~ItemMoment()
  {
    psc_mfields_list_del(&psc_mfields_base_list, &mres_base_);
    psc_mfields_destroy(mres_base_);
  }

  void run(PscMfieldsBase mflds_base, PscMparticlesBase mprts_base) override
  {
    mparticles_t mprts = mprts_base.get_as<mparticles_t>();

    moment_.run(mfields_t{mres_base_}, mprts);
    
    mprts.put_as(mprts_base, MP_DONT_COPY);
  }

  virtual psc_mfields* mres() override { return mres_base_; }
  
private:
  Moment_t moment_;
  psc_mfields* mres_base_;
};

// ----------------------------------------------------------------------
// FieldsItemMomentOps
  
template<typename Moment_t>
using FieldsItemMomentOps = FieldsItemOps<ItemMoment<ItemMomentLoopPatches<Moment_t>>>;

