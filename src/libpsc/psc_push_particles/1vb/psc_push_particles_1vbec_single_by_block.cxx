
#include "psc_push_particles_private.h"

#include "psc_particles_as_single_by_block.h"
#include "psc_fields_as_single.h"
#include "1vb/psc_push_particles_1vb.h"

// ======================================================================
// psc_push_particles: subclass "1vbec_single_by_block"

struct psc_push_particles_ops psc_push_particles_1vbec_single_by_block_ops = {
  .name                  = "1vbec_single_by_block",
  .push_mprts_yz         = psc_push_particles_1vbec_single_by_block_push_mprts_yz,
  .push_mprts_xyz        = psc_push_particles_1vbec_single_by_block_push_mprts_xyz,
  .particles_type        = PARTICLE_TYPE,
  .fields_type           = FIELDS_TYPE,
};

