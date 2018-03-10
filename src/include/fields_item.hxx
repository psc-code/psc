
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
  virtual void run2(PscMfieldsBase mflds_base, PscMparticlesBase mprts_base) = 0;

  template<typename MF>
  void create_mres_(MPI_Comm comm, int n_comps)
  {
    mres_base_ = psc_mfields_create(comm);
    psc_mfields_set_type(mres_base_, fields_traits<typename MF::fields_t>::name);
    psc_mfields_set_param_int(mres_base_, "nr_fields", n_comps);
    psc_mfields_set_param_int3(mres_base_, "ibn", ppsc->ibn);
    mres_base_->grid = &ppsc->grid();
    psc_mfields_setup(mres_base_);
    psc_mfields_list_add(&psc_mfields_base_list, &mres_base_);
  }
  
  ~FieldsItemBase()
  {
    psc_mfields_list_del(&psc_mfields_base_list, &mres_base_);
    psc_mfields_destroy(mres_base_);
  }
  
  psc_mfields* mres_base_;
};

template<class Derived>
struct FieldsItemCRTP : FieldsItemBase
{
  using mres_t = MfieldsC; // default (FIXME, get rid of?)

  FieldsItemCRTP(MPI_Comm comm, PscBndBase bnd)
    : bnd_(bnd)
  {
    auto d = static_cast<Derived*>(this);

    create_mres_<typename Derived::mres_t>(comm, d->n_comps());
    auto fld_names = d->fld_names();
    for (int m = 0; m < d->n_comps(); m++) {
      psc_mfields_set_comp_name(mres_base_, m, fld_names[m]);
    }
  }

  void run2(PscMfieldsBase mflds_base, PscMparticlesBase mprts_base) override
  {
    auto d = static_cast<Derived*>(this);
    d->run(mflds_base, mprts_base);

    if (d->flags & POFI_ADD_GHOSTS) {
      bnd_.add_ghosts(mres_base_, 0, mres_base_->nr_fields);
    }
  }

private:
  PscBndBase bnd_;
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
    struct psc_output_fields_item_ops *ops = psc_output_fields_item_ops(item_);
    
    assert(mres.mflds() == sub()->mres_base_);
    sub()->run2(mflds, mprts);
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

    if (!item->mres_base_) return; // ctor hadn't run yet FIXME

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
    flags     = Item_t::flags;
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
    if (std::is_same<mfields_t, PscMfieldsC>::value) {
      return Item::name;
    } else {
      return strdup((std::string{Item::name} + "_" + fields_traits<typename mfields_t::fields_t>::name).c_str());
    }
  }
  constexpr static int flags = 0;
 
  FieldsItemFields(MPI_Comm comm, PscBndBase bnd)
  {
    create_mres_<mfields_t>(comm, Item::n_comps);
    for (int m = 0; m < Item::n_comps; m++) {
      psc_mfields_set_comp_name(mres_base_, m, Item::fld_names()[m]);
    }
  }

  void run(PscMfieldsBase mflds_base, PscMparticlesBase mprts_base) override
  {
    mfields_t mflds = mflds_base.get_as<mfields_t>(JXI, JXI + 3);
    mfields_t mres(mres_base_);
    Item::run(mflds, mres);
    mflds.put_as(mflds_base, 0, 0);
  }

  void run2(PscMfieldsBase mflds_base, PscMparticlesBase mprts_base) override
  {
    run(mflds_base, mprts_base);
  }
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
// ItemMomentWrap

template<typename Moment_t, typename mparticles_t, typename mfields_t>
struct ItemMomentWrap
{
  using fields_t = typename mfields_t::fields_t;
  using Fields = Fields3d<fields_t>;
  
  static void run(mparticles_t mprts, mfields_t mflds_res)
  {
    for (int p = 0; p < mprts->n_patches(); p++) {
      mflds_res[p].zero();
      Moment_t::run(mflds_res[p], mprts[p]);
      add_ghosts_boundary(mflds_res[p], p, 0, mflds_res->n_comps());
    }
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
};

// ----------------------------------------------------------------------
// ItemMoment

template<typename Moment_t, typename mparticles_t, typename mfields_t>
struct ItemMoment : FieldsItemCRTP<ItemMoment<Moment_t, mparticles_t, mfields_t>>
{
  using Base = FieldsItemCRTP<ItemMoment<Moment_t, mparticles_t, mfields_t>>;
  using Base::Base;
  
  static const char* name()
  {
    return strdup((std::string(Moment_t::name) + "_" +
		   mparticles_traits<mparticles_t>::name).c_str());
  }
  constexpr static int flags = Moment_t::flags;

  static int n_comps()
  {
    int n_comps = Moment_t::n_comps;
    if (flags & POFI_BY_KIND) {
      n_comps *= POFI_BY_KIND;
    }
    assert(n_comps <= POFI_MAX_COMPS);
    return n_comps;
  }
  static fld_names_t fld_names()
  {
    auto fld_names = Moment_t::fld_names();
    if (!(flags & POFI_BY_KIND)) {
      return fld_names;
    }

    static fld_names_t names;
    if (!names[0]) {
      for (int k = 0; k < ppsc->nr_kinds; k++) {
	for (int m = 0; m < Moment_t::n_comps; m++) {
	  auto s = std::string(fld_names[m]) + "_" + ppsc->kinds[k].name;
	  names[k * Moment_t::n_comps + m] = strdup(s.c_str());
	}
      }
    }
    return names;
  }

  void run(PscMfieldsBase mflds_base, PscMparticlesBase mprts_base) override
  {
    mparticles_t mprts = mprts_base.get_as<mparticles_t>();
    mfields_t mres = this->mres_base_->template get_as<mfields_t>(0, 0);

    ItemMomentWrap<Moment_t, mparticles_t, mfields_t>::run(mprts, mres);
    
    mres.put_as(this->mres_base_, 0, this->mres_base_->nr_fields);
    mprts.put_as(mprts_base, MP_DONT_COPY);
  }
};

// ----------------------------------------------------------------------
// FieldsItemMomentOps
  
template<typename Item_t, typename mparticles_t, typename mfields_t>
using FieldsItemMomentOps = FieldsItemOps<ItemMoment<Item_t, mparticles_t, mfields_t>>;

