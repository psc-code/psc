
#pragma once

#include "psc_push_fields_private.h"

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

  void advance_H(struct psc_mfields *mfields, double frac)
  {
    psc_push_fields_push_H(pushf_, mfields, frac);
  }

  void advance_b2(struct psc_mfields *mfields)
  {
    psc_push_fields_step_b2(pushf_, mfields);
  }

  void advance_a(struct psc_mfields *mfields)
  {
    psc_push_fields_step_a(pushf_, mfields);
  }

  psc_push_fields *pushf() { return pushf_; }
  
  sub_t* operator->() { return sub_; }

  sub_t* sub() { return sub_; }

private:
  psc_push_fields *pushf_;
  sub_t *sub_;
};

