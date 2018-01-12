
#include "psc_balance_private.h"

#include "psc_particles_as_double.h"
#include "psc_fields_as_c.h"

#include "psc_balance_common.cxx"

// ======================================================================
// psc_balance subclass "double"

struct psc_balance_ops psc_balance_double_ops = {
  name                  : "double",
  mprts_type            : "double",
  mflds_type            : "c",
  communicate_particles : psc_balance_sub_communicate_particles,
  communicate_fields    : psc_balance_sub_communicate_fields
};
