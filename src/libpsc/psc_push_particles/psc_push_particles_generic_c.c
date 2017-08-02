
#include "psc_push_particles_private.h"
#include "psc_generic_c.h"

// ======================================================================
// psc_push_particles: subclass "generic_c"

struct psc_push_particles_ops psc_push_particles_generic_c_ops = {
  .name                  = "generic_c",
  .push_mprts_y          = psc_push_particles_generic_c_push_mprts_y,
  .push_mprts_z          = psc_push_particles_generic_c_push_mprts_z,
  .push_mprts_xy         = psc_push_particles_generic_c_push_mprts_xy,
  .push_mprts_xz         = psc_push_particles_generic_c_push_mprts_xz,
  .push_mprts_yz         = psc_push_particles_generic_c_push_mprts_yz,
  .push_mprts_xyz        = psc_push_particles_generic_c_push_mprts_xyz,
  .particles_type        = "c",
  .fields_type           = "c",
};
