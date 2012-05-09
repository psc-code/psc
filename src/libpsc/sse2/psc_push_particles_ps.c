
#include "psc_push_particles_private.h"
#include "psc_push_particles_ps.h"

// ======================================================================
// psc_push_particles: subclass "ps_1vb"

struct psc_push_particles_ops psc_push_particles_ps_1vb_ops = {
  .name                  = "ps_1vb",
  .push_yz               = psc_push_particles_ps_1vb_push_yz,
};

