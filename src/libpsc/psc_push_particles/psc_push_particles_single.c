
#include "psc_push_particles_private.h"
#include "psc_push_particles_single.h"

// ======================================================================
// psc_push_particles: subclass "1vb_single"

struct psc_push_particles_ops psc_push_particles_1vb_single_ops = {
  .name                  = "1vb_single",
  .push_a_yz             = psc_push_particles_1vb_single_push_a_yz,
};

// ======================================================================
// psc_push_particles: subclass "1vb2_single"

struct psc_push_particles_ops psc_push_particles_1vb2_single_ops = {
  .name                  = "1vb2_single2",
  .push_a_yz             = psc_push_particles_1vb2_single_push_a_yz,
};

