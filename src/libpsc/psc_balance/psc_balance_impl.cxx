
#include "psc_balance_private.h"

#include "psc_balance_impl.hxx"

#include "psc_particles_single.h"
#include "psc_particles_as_double.h"
#include "psc_fields_single.h"
#include "psc_fields_as_c.h"

// ======================================================================
// psc_balance subclass "single"

struct psc_balance_ops_single : psc_balance_ops
{
  using PscBalance_t = PscBalance_<PscMparticlesSingle, PscMfieldsSingle>;
  psc_balance_ops_single() {
    name                  = "single";
    mprts_type            = "single";
    mflds_type            = "single";
    communicate_particles = PscBalance_t::communicate_particles;
    communicate_fields    = PscBalance_t::communicate_fields;
  }
} psc_balance_single_ops;

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
