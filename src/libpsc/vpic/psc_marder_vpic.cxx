
#include "marder_vpic.hxx"

// ----------------------------------------------------------------------
// psc_marder: subclass "vpic"

struct psc_marder_ops_vpic : psc_marder_ops {
  using Wrapper = MarderWrapper<MarderVpic>;
  psc_marder_ops_vpic() {
    name                  = "vpic";
    size                  = Wrapper::size;
    setup                 = Wrapper::setup;
    destroy               = Wrapper::destroy;
  }
} psc_marder_vpic_ops;

