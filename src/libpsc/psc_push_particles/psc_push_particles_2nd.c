
#include "psc_push_particles_private.h"
#include "psc_push_particles_2nd.h"

// ======================================================================
// psc_push_particles: subclass "2nd_double"
//
// 2nd order, Esirkepov, like the old generic_c/fortran, but on "double"
// particles

struct psc_push_particles_ops psc_push_particles_2nd_double_ops = {
  .name                  = "2nd_double",
  .push_a_yz             = psc_push_particles_2nd_double_push_a_yz,
  .particles_type        = "double",
  .fields_type           = "c",
};

