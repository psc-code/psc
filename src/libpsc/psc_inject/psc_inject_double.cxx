
#include <psc_particles_double.h>
#include <psc_fields_c.h>

#include "psc_inject_common.cxx"

// ----------------------------------------------------------------------
// psc_inject "double"

struct psc_inject_ops_double : psc_inject_ops {
  using Inject_t = Inject_<MparticlesDouble, MfieldsC>;
  psc_inject_ops_double() {
    name                = "double";
    create              = Inject_t::create;
    run                 = Inject_t::run;
  }
} psc_inject_ops_double;

