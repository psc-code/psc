
#include "bnd_fields_vpic.hxx"

// ----------------------------------------------------------------------
// psc_bnd_fields: subclass "vpic"

struct psc_bnd_fields_ops_vpic : psc_bnd_fields_ops {
  using Wrapper = PscBndFieldsWrapper<BndFieldsVpic>;
  psc_bnd_fields_ops_vpic() {
    name                  = "vpic";
    size                  = Wrapper::size;
    setup                 = Wrapper::setup;
    destroy               = Wrapper::destroy;
  }
} psc_bnd_fields_vpic_ops;


