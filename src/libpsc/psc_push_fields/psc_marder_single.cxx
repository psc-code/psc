
#include "psc_fields_single.h"
#include "psc_particles_as_single.h"

#include "psc_marder_common.cxx"

// ======================================================================
// psc_marder: single

using marder_ops_single = marder_ops<PscMfieldsSingle>;

struct psc_marder_ops_single : psc_marder_ops {
  psc_marder_ops_single() {
    name                  = "single";
    setup                 = marder_ops_single::setup;
    destroy               = marder_ops_single::destroy;
    correct               = marder_ops_single::correct;
  }
} psc_marder_single_ops;

