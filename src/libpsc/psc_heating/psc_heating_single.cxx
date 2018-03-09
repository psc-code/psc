
#include "psc_particles_single.h"

#include "psc_heating_common.cxx"

// ----------------------------------------------------------------------
// psc_heating "single"

struct psc_heating_ops_single : psc_heating_ops {
  using PscHeating_t = PscHeatingWrapper<Heating_<MparticlesSingle>>;
  psc_heating_ops_single() {
    name                = "single";
    size                = PscHeating_t::size;
    setup               = PscHeating_t::setup;
    destroy             = PscHeating_t::destroy;
  }
} psc_heating_ops_single;


