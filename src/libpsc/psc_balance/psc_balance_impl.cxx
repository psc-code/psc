
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
  using Balance_t = Balance_<PscMparticlesSingle, PscMfieldsSingle>;
  using Wrapper_t = BalanceWrapper<Balance_t>;
  psc_balance_ops_single() {
    name                  = "single";
    mprts_type            = "single";
    mflds_type            = "single";
    size                  = Wrapper_t::size;
    setup                 = Wrapper_t::setup;
    destroy               = Wrapper_t::destroy;
    communicate_particles = Balance_t::communicate_particles;
    communicate_fields    = Balance_t::communicate_fields;
  }
} psc_balance_single_ops;

// ======================================================================
// psc_balance subclass "double"

struct psc_balance_ops_double : psc_balance_ops
{
  using Balance_t = Balance_<PscMparticlesDouble, PscMfieldsC>;
  using Wrapper_t = BalanceWrapper<Balance_t>;
  psc_balance_ops_double() {
    name                  = "double";
    mprts_type            = "double";
    mflds_type            = "c";
    size                  = Wrapper_t::size;
    setup                 = Wrapper_t::setup;
    destroy               = Wrapper_t::destroy;
    communicate_particles = Balance_t::communicate_particles;
    communicate_fields    = Balance_t::communicate_fields;
  }
} psc_balance_double_ops;
