
#include "psc_push_particles_private.h"
#include "psc_1p5_double.h"

// ======================================================================
// psc_push_particles: subclass "1p5_double"

struct psc_push_particles_ops psc_push_particles_1p5_double_ops = {
  .name                  = "1p5_double",
  .push_mprts_yz         = psc_push_particles_1p5_double_push_mprts_yz,
  .particles_type        = "double",
  .fields_type           = "c",
};
