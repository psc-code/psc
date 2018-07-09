
#include "psc_push_fields_private.h"

#include "push_fields_vpic.hxx"

// ----------------------------------------------------------------------
// psc_push_fields: subclass "vpic"

struct psc_push_fields_ops_vpic : psc_push_fields_ops {
  using Wrapper = PscPushFieldsWrapper<PushFieldsVpic>;
  psc_push_fields_ops_vpic() {
    name                  = "vpic";
    size                  = Wrapper::size;
    setup                 = Wrapper::setup;
    destroy               = Wrapper::destroy;
  }
} psc_push_fields_vpic_ops;
