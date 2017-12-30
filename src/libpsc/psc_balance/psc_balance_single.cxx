
#include "psc_balance_private.h"

#include "psc_particles_as_single.h"
#include "psc_fields_as_single.h"

#include "psc_balance_common.cxx"

// ======================================================================
// psc_balance subclass "single"

struct psc_balance_ops psc_balance_single_ops = {
  .name                  = "single",
  .mprts_type            = "single",
  .mflds_type            = "single",
  .communicate_particles = psc_balance_sub_communicate_particles,
  .communicate_fields    = psc_balance_sub_communicate_fields,
};
