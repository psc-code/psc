
#include "psc_bnd_private.h"

#include "bnd.hxx"

struct BndVpic : BndBase
{
  BndVpic(const Grid_t& grid, mrc_domain *domain, int ibn[3])
  {}

  void reset() override
  {}
  
  void fill_ghosts(PscMfieldsBase mflds_base, int mb, int me) override
  {}
  
  void add_ghosts(PscMfieldsBase mflds_base, int mb, int me) override
  {}
};

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

