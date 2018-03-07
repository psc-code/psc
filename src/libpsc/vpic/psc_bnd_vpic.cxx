
#include "psc_bnd_private.h"

#include "bnd.hxx"

struct BndVpic : BndBase
{
  BndVpic(const Grid_t& grid, mrc_domain *domain, int ibn[3])
  {}

  void reset() override
  {}
  
  void fill_ghosts(struct psc_mfields *mflds_base, int mb, int me) override
  {}
  
  void add_ghosts(struct psc_mfields *mflds_base, int mb, int me) override
  {}
};

// ----------------------------------------------------------------------
// psc_bnd_vpic_fill_ghosts

static void
psc_bnd_vpic_fill_ghosts(struct psc_bnd *bnd, struct psc_mfields *mflds,
			 int mb, int me)
{  
}

// ----------------------------------------------------------------------
// psc_bnd_vpic_add_ghosts

static void
psc_bnd_vpic_add_ghosts(struct psc_bnd *bnd, struct psc_mfields *mflds,
			int mb, int me)
{  
}

// ----------------------------------------------------------------------
// psc_bnd: subclass "vpic"

struct psc_bnd_ops_vpic : psc_bnd_ops {
  using PscBnd_t = PscBndWrapper<BndVpic>;
  psc_bnd_ops_vpic() {
    name                    = "vpic";
    size                    = PscBnd_t::size;
    reset                   = PscBnd_t::reset;
    setup                   = PscBnd_t::setup;
    destroy                 = PscBnd_t::destroy;
    add_ghosts              = PscBnd_t::add_ghosts;
    fill_ghosts             = PscBnd_t::fill_ghosts;
  }
} psc_bnd_vpic_ops;

