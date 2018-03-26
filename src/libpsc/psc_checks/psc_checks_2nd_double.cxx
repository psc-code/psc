
#include "psc_fields_as_c.h"
#include "psc_particles_as_double.h"

#define PSC_CHECKS_ORDER "2nd"

#include "psc_checks_common.cxx"

// ----------------------------------------------------------------------
// psc_checks_2nd_double_ops

struct psc_checks_2nd_double_ops : psc_checks_ops {
  using Wrapper_t = ChecksWrapper<Checks_<MparticlesDouble, MfieldsC>>;
  psc_checks_2nd_double_ops() {
    name                            = PSC_CHECKS_ORDER "_" PARTICLE_TYPE;
    size                            = Wrapper_t::size;
    setup                           = Wrapper_t::setup;
    destroy                         = Wrapper_t::destroy;
    //read                            = psc_checks_sub_read;
  }
} psc_checks_2nd_double_ops;

