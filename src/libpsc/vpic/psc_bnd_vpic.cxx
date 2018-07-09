
#include "bnd_vpic.hxx"

#include "psc_bnd_private.h"

// ----------------------------------------------------------------------
// psc_bnd: subclass "vpic"

struct psc_bnd_ops_vpic : psc_bnd_ops {
  using PscBnd_t = PscBndWrapper<BndVpic>;
  psc_bnd_ops_vpic() {
    name                    = "vpic";
    size                    = PscBnd_t::size;
    setup                   = PscBnd_t::setup;
    destroy                 = PscBnd_t::destroy;
  }
} psc_bnd_vpic_ops;

