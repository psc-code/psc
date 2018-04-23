
#include "psc_push_fields_private.h"

#include "psc_fields_single.h"
#include "psc_fields_c.h"

#include "push_fields.hxx"
#include "psc_push_fields_impl.hxx"

// ======================================================================
// psc_push_fields: subclass "single"

struct psc_push_fields_ops_single : psc_push_fields_ops {
  using Wrapper = PscPushFieldsWrapper<PushFields<MfieldsSingle>>;
  psc_push_fields_ops_single() {
    name                  = "single";
    size                  = Wrapper::size;
    setup                 = Wrapper::setup;
    destroy               = Wrapper::destroy;
  }
} psc_push_fields_single_ops;

// ======================================================================
// psc_push_fields: subclass "c"

struct psc_push_fields_ops_c : psc_push_fields_ops {
  using Wrapper = PscPushFieldsWrapper<PushFields<MfieldsC>>;
  psc_push_fields_ops_c() {
    name                  = "single";
    size                  = Wrapper::size;
    setup                 = Wrapper::setup;
    destroy               = Wrapper::destroy;
  }
} psc_push_fields_c_ops;
