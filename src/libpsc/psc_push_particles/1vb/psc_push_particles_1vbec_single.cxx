
#include "psc_push_particles_private.h"

#include "psc_particles_as_single.h"
#include "psc_fields_as_single.h"
#include "fields.hxx"
#include "1vb/psc_push_particles_1vb.h"

// ======================================================================
// psc_push_particles: subclass "1vbec_single"

struct psc_push_particles_ops psc_push_particles_1vbec_single_ops = {
  .name                  = "1vbec_single",
  .push_mprts_xyz        = psc_push_particles_1vbec_single_push_mprts_xyz,
  .push_mprts_xz         = psc_push_particles_1vbec_single_push_mprts_xyz_xz,
  .push_mprts_yz         = psc_push_particles_1vbec_single_push_mprts_yz,
  .push_mprts_1          = psc_push_particles_1vbec_single_push_mprts_1,
  .stagger_mprts_yz      = push_particles_ops<mfields_single_t, dim_yz>::stagger_mprts,
  .particles_type        = PARTICLE_TYPE,
};

