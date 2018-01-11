
#include "psc_push_particles_private.h"
#include "psc_push_particles_2nd.h"

// ======================================================================
// psc_push_particles: subclass "2nd_double"
//
// 2nd order, Esirkepov, like the old generic_c/fortran, but on "double"
// particles

struct psc_push_particles_ops_2nd_double : psc_push_particles_ops {
  psc_push_particles_ops_2nd_double() {
    name                  = "2nd_double";
    push_mprts_yz         = psc_push_particles_2nd_double_push_mprts_yz;
    particles_type        = "double";
  }
} psc_push_particles_2nd_double_ops;

