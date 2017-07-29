
#include "psc_push_particles_private.h"
#include "psc_generic_c.h"

// ======================================================================
// psc_push_particles: subclass "generic_c"

struct psc_push_particles_ops psc_push_particles_generic_c_ops = {
  .name                  = "generic_c",
  .push_a_y              = psc_push_particles_generic_c_push_a_y,
  .push_a_z              = psc_push_particles_generic_c_push_a_z,
  .push_a_xy             = psc_push_particles_generic_c_push_a_xy,
  .push_a_xz             = psc_push_particles_generic_c_push_a_xz,
  .push_a_yz             = psc_push_particles_generic_c_push_a_yz,
  .push_a_xyz            = psc_push_particles_generic_c_push_a_xyz,
  .particles_type        = "c",
  .fields_type           = "c",
};
