
#include "psc_push_particles_private.h"
#include "psc_push_particles_double.h"

// ======================================================================
// psc_push_particles: subclass "1vb_double"

struct psc_push_particles_ops psc_push_particles_1vb_double_ops = {
  .name                  = "1vb_double",
  .push_a_yz             = psc_push_particles_1vb_double_push_a_yz,
};

// ======================================================================
// psc_push_particles: subclass "1vbec_double"

struct psc_push_particles_ops psc_push_particles_1vbec_double_ops = {
  .name                  = "1vbec_double",
  .push_a_yz             = psc_push_particles_1vbec_double_push_a_yz,
};

// ======================================================================
// psc_push_particles: subclass "1vbec3d_double"

struct psc_push_particles_ops psc_push_particles_1vbec3d_double_ops = {
  .name                  = "1vbec3d_double",
  .push_a_yz             = psc_push_particles_1vbec3d_double_push_a_yz,
};

