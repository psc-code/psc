
#include "psc_push_particles_private.h"
#include "psc_push_particles_single.h"

// ======================================================================
// psc_push_particles: subclass "single_1vb"

struct psc_push_particles_ops psc_push_particles_single_1vb_ops = {
  .name                  = "single_1vb",
  .push_a_yz             = psc_push_particles_single_1vb_push_a_yz,
};

