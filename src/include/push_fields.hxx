
#pragma once

#include "psc_push_fields_private.h"
#include "fields3d.hxx"

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
    psc_push_fields_step_a(pushf_, mflds.mflds());
  }

  psc_push_fields *pushf() { return pushf_; }
  
  sub_t* operator->() { return sub_; }

  sub_t* sub() { return sub_; }

private:
  psc_push_fields *pushf_;
  sub_t *sub_;
};

