
#include "psc_push_particles_private.h"

#include "psc_particles_as_single.h"
#include "psc_fields_as_single.h"

#include "../inc_defs.h"

#define DIM DIM_YZ
#define CALC_J CALC_J_1VB_VAR1
#define F3_CURR F3_S
#define F3_CACHE F3_S
#define F3_CACHE_TYPE "single"
#define INTERPOLATE_1ST INTERPOLATE_1ST_EC
#define VB_2D

#include "1vb_yz.c"

// ======================================================================
// psc_push_particles: subclass "1vbec_single"

struct psc_push_particles_ops psc_push_particles_1vbec_single_ops = {
  .name                  = "1vbec_single",
  .push_a_yz             = psc_push_particles_push_a_yz,
  .particles_type        = PARTICLE_TYPE,
  .fields_type           = FIELDS_TYPE,
};

