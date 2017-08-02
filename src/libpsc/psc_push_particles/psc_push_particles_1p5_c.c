
#include "psc_push_particles_private.h"
#include "psc_1p5_c.h"

// ======================================================================
// psc_push_particles: subclass "1p5_c"

struct psc_push_particles_ops psc_push_particles_1p5_c_ops = {
  .name                  = "1p5_c",
  .push_mprts_yz         = psc_push_particles_1p5_c_push_mprts_yz,
  .particles_type        = "c",
  .fields_type           = "c",
};
