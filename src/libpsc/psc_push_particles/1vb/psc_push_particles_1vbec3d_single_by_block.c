
#include "psc_push_particles_private.h"

#include "psc_particles_as_single_by_block.h"
#include "psc_fields_as_single.h"

#include "../inc_defs.h"

#define DIM DIM_YZ
#define CALC_J CALC_J_1VB_VAR1
#define F3_CACHE F3_S
#define F3_CACHE_TYPE "single"
#define INTERPOLATE_1ST INTERPOLATE_1ST_EC
#define PUSHER_BY_BLOCK

#include "1vb_yz.c"

// ======================================================================
// psc_push_particles: subclass "1vbec3d_single_by_block"

struct psc_push_particles_ops psc_push_particles_1vbec3d_single_by_block_ops = {
  .name                  = "1vbec3d_single_by_block",
  .push_a_yz             = psc_push_particles_push_a_yz,
  .particles_type        = PARTICLE_TYPE,
  .fields_type           = FIELDS_TYPE,
};

