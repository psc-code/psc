
#pragma once

#include "psc_output_fields_item_private.h"

#include "particles.hxx"
#include "fields3d.hxx"
#include "bnd.hxx"

#include <mrc_profile.h>

// ======================================================================
// FieldsItemBase

struct FieldsItemBase
{
  virtual void run(PscMfieldsBase mflds_base, PscMparticlesBase mprts_base) = 0;
  virtual void run2(PscMfieldsBase mflds_base, PscMparticlesBase mprts_base) = 0;

  psc_mfields* mres_base_;
};

template<class Derived>
struct FieldsItemCRTP : FieldsItemBase
{
  void run2(PscMfieldsBase mflds_base, PscMparticlesBase mprts_base) override
  {
    auto d = static_cast<Derived*>(this);
    d->run(mflds_base, mprts_base);
  }
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

    if (ops->flags & POFI_ADD_GHOSTS) {
      assert(item_->bnd);
      auto bnd = PscBndBase(item_->bnd);
      bnd.add_ghosts(mres, 0, mres.n_fields());
    }
  }

  sub_t* sub() { return mrc_to_subobj(item_, sub_t); }
  sub_t* operator->() { return sub(); }

private:
  psc_output_fields_item* item_;
};

using PscFieldsItemBase = PscFieldsItem<FieldsItemBase>;

// ----------------------------------------------------------------------
// psc_output_fields_item_create_mfields

static struct psc_mfields *
psc_output_fields_item_create_mfields(struct psc_output_fields_item *_item)
{
  struct psc_output_fields_item_ops *ops = psc_output_fields_item_ops(_item);
  MPI_Comm comm = psc_output_fields_item_comm(_item);
  const char* type = "c";
  if (strcmp(psc_output_fields_item_type(_item), "n_1st_cuda") == 0) { // FIXME
    type = "cuda";
  }
  int n_comps = ops->nr_comp;
  int n_comps_total = n_comps;
  if (ops->flags & POFI_BY_KIND) {
    n_comps_total *= ppsc->nr_kinds;
  }

  struct psc_mfields *flds = psc_mfields_create(comm);
  psc_mfields_set_type(flds, type);
  psc_mfields_set_param_int(flds, "nr_fields", n_comps_total);
  psc_mfields_set_param_int3(flds, "ibn", ppsc->ibn);
  flds->grid = &ppsc->grid();
  psc_mfields_setup(flds);
  assert(ops->nr_comp <= POFI_MAX_COMPS);
  for (int m = 0; m < n_comps_total; m++) {
    if (ops->flags & POFI_BY_KIND) {
      int mm = m % ops->nr_comp;
      int k = m / ops->nr_comp;
      char s[strlen(ops->fld_names[mm]) + strlen(ppsc->kinds[k].name) + 2];
      sprintf(s, "%s_%s", ops->fld_names[mm], ppsc->kinds[k].name);
      psc_mfields_set_comp_name(flds, m, s);
    } else {
      psc_mfields_set_comp_name(flds, m, ops->fld_names[m]);
    }
  }

  return flds;
}

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
    new(item.sub()) FieldsItem();

    item->mres_base_ = psc_output_fields_item_create_mfields(_item);
    psc_mfields_list_add(&psc_mfields_base_list, &item->mres_base_);
  }

  static void destroy(psc_output_fields_item* _item)
  {
    PscFieldsItem<FieldsItem> item(_item);

    if (!item->mres_base_) return; // ctor hadn't run yet FIXME

    item->~FieldsItem();

    psc_mfields_list_del(&psc_mfields_base_list, &item->mres_base_);
    psc_mfields_destroy(item->mres_base_);
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
    nr_comp   = Item_t::n_comps;
    fld_names = Item_t::fld_names();
    flags     = Item_t::flags;
  }
};

