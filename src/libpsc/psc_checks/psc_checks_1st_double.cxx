
#include "psc_fields_c.h"
#include "psc_particles_double.h"

#define PSC_CHECKS_ORDER "1st"

#include "psc_checks_common.cxx"

// ----------------------------------------------------------------------
// psc_checks_1st_double_ops

struct psc_checks_1st_double_ops : psc_checks_ops {
  using Wrapper_t = ChecksWrapper<Checks_<MparticlesDouble, MfieldsC, checks_order_1st>>;
  psc_checks_1st_double_ops() {
    name                            = "1st_double";
    size                            = Wrapper_t::size;
    setup                           = Wrapper_t::setup;
    destroy                         = Wrapper_t::destroy;
    //read                            = psc_checks_sub_read;
  }
} psc_checks_1st_double_ops;



