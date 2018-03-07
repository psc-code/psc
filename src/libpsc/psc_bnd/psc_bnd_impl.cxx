
#include "psc_bnd_private.h"
#include "psc_fields_c.h"
#include "psc_fields_single.h"

#include "psc_bnd_impl.hxx"

// ======================================================================
// psc_bnd: subclass "c"

struct psc_bnd_ops_c : psc_bnd_ops {
  using Bnd_t = Bnd_<PscMfieldsC>;
  using PscBnd_t = PscBndWrapper<Bnd_t>;
  psc_bnd_ops_c() {
    name                    = "c";
    size                    = PscBnd_t::size;
    setup                   = PscBnd_t::setup;
    destroy                 = PscBnd_t::destroy;
    reset                   = PscBnd_t::reset;
    create_ddc              = Bnd_t::create;
    add_ghosts              = Bnd_t::add_ghosts;
    fill_ghosts             = Bnd_t::fill_ghosts;
  }
} psc_bnd_c_ops;

// ======================================================================
// psc_bnd: subclass "single"


struct psc_bnd_ops_single : psc_bnd_ops {
  using Bnd_t = Bnd_<PscMfieldsSingle>;
  using PscBnd_t = PscBndWrapper<Bnd_t>;
  psc_bnd_ops_single() {
    name                    = "single";
    size                    = PscBnd_t::size;
    setup                   = PscBnd_t::setup;
    destroy                 = PscBnd_t::destroy;
    reset                   = PscBnd_t::reset;
    create_ddc              = Bnd_t::create;
    add_ghosts              = Bnd_t::add_ghosts;
    fill_ghosts             = Bnd_t::fill_ghosts;
  }
} psc_bnd_single_ops;

