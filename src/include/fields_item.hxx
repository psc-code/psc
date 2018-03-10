
#pragma once

#include "psc_output_fields_item_private.h"

#include "particles.hxx"
#include "fields3d.hxx"
#include "bnd.hxx"
#include "psc_fields_c.h"

#include <mrc_profile.h>

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
  
  constexpr static char const* name() { return Item::name; }
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

