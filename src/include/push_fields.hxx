
#pragma once

#include "psc_push_fields_private.h"
#include "fields3d.hxx"

#include "psc_bnd_fields.h"
#include "psc_bnd.h"

// ======================================================================
// PscPushFieldsBase

template<typename S>
struct PscPushFieldsBase
{
  using sub_t = S;

  explicit PscPushFieldsBase(psc_push_fields *pushf)
    : pushf_(pushf),
      sub_(mrc_to_subobj(pushf, sub_t))
  {}

  void advance_H(mfields_base_t mflds, double frac)
  {
    psc_push_fields_push_H(pushf_, mflds.mflds(), frac);
  }

  void advance_b2(mfields_base_t mflds)
  {
    psc_push_fields_step_b2(pushf_, mflds.mflds());
  }

  void advance_a(mfields_base_t mflds)
  {
    if (pushf_->variant == 0) {
      psc_bnd_fields_fill_ghosts_E(pushf_->bnd_fields, mflds.mflds());
      psc_bnd_fill_ghosts(ppsc->bnd, mflds.mflds(), EX, EX + 3);
    }
    
    // push H
    psc_push_fields_push_H(pushf_, mflds.mflds(), .5);
    psc_bnd_fields_fill_ghosts_H(pushf_->bnd_fields, mflds.mflds());
    if (pushf_->variant == 0) {
      psc_bnd_fill_ghosts(ppsc->bnd, mflds.mflds(), HX, HX + 3);
    }
  }

  psc_push_fields *pushf() { return pushf_; }
  
  sub_t* operator->() { return sub_; }

  sub_t* sub() { return sub_; }

private:
  psc_push_fields *pushf_;
  sub_t *sub_;
};

