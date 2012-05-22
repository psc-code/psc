
#include "psc_push_particles_private.h"
#include "psc_generic_c.h"

// ======================================================================
// psc_push_particles: subclass "generic_c"

struct psc_push_particles_ops psc_push_particles_generic_c_ops = {
  .name                  = "generic_c",
  .push_z                = psc_push_particles_generic_c_push_z,
  .push_xy               = psc_push_particles_generic_c_push_xy,
  .push_xz               = psc_push_particles_generic_c_push_xz,
  .push_a_yz             = psc_push_particles_generic_c_push_a_yz,
  .push_xyz              = psc_push_particles_generic_c_push_xyz,
  .calc_j_z              = psc_push_particles_generic_c_calc_j_z,
};
