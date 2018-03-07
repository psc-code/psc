
#include <psc_particles_single.h>
#include <psc_particles_double.h>
#include <psc_fields_c.h>

#include "psc_inject_impl.hxx"

// ----------------------------------------------------------------------
// psc_inject "single"

struct psc_inject_ops_single : psc_inject_ops {
  using Inject_t = Inject_<MparticlesSingle, MfieldsC>;
  psc_inject_ops_single() {
    name                = "single";
    create              = Inject_t::create;
    run                 = Inject_t::run;
  }
} psc_inject_ops_single;

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

