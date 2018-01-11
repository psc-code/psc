
#include "psc_fields_c.h"
#include "psc_particles_as_double.h"

#include "psc_marder_common.cxx"

// ======================================================================
// psc_marder: c

using marder_ops_c = marder_ops<mfields_c_t>;

struct psc_marder_ops_c : psc_marder_ops {
  psc_marder_ops_c() {
    name                  = "c";
    setup                 = marder_ops_c::setup;
    destroy               = marder_ops_c::destroy;
    correct               = marder_ops_c::correct;
  }
} psc_marder_c_ops;
