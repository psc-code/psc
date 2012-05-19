
#include "psc_push_particles_private.h"
#include "psc_push_particles_double.h"

// ======================================================================
// psc_push_particles: subclass "double_1vb"

struct psc_push_particles_ops psc_push_particles_double_1vb_ops = {
  .name                  = "double_1vb",
  .push_yz               = psc_push_particles_double_1vb_push_yz,
};

