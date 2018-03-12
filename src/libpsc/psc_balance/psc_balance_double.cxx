
#include "psc_balance_private.h"

#include "psc_particles_as_double.h"
#include "psc_fields_as_c.h"

#include "psc_balance_common.cxx"

// ======================================================================
// psc_balance subclass "double"

struct psc_balance_ops_double : psc_balance_ops
{
  using PscBalance_t = PscBalance_<PscMparticlesDouble, PscMfieldsC>;
  psc_balance_ops_double() {
    name                  = "double";
    mprts_type            = "double";
    mflds_type            = "c";
    communicate_particles = PscBalance_t::communicate_particles;
    communicate_fields    = PscBalance_t::communicate_fields;
  }
} psc_balance_double_ops;
