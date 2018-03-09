
#include "psc_particles_double.h"

#include "psc_heating_common.cxx"

// ----------------------------------------------------------------------
// psc_heating "double"

struct psc_heating_ops_double : psc_heating_ops {
  using PscHeating_t = PscHeatingWrapper<Heating_<MparticlesDouble>>;
  psc_heating_ops_double() {
    name                = "double";
    size                = PscHeating_t::size;
    setup               = PscHeating_t::setup;
    destroy             = PscHeating_t::destroy;
  }
} psc_heating_ops_double;



