
#include "psc_push_particles_private.h"

#include "psc_particles_as_double.h"
#include "psc_fields_as_c.h"

#include "../inc_defs.h"

#define DIM DIM_YZ
#define CALC_J CALC_J_1VB_2D
#define INTERPOLATE_1ST INTERPOLATE_1ST_STD

#define psc_push_particles_push_mprts_yz psc_push_particles_1vb_double_push_mprts_yz

#include "1vb.c"

// ======================================================================
// psc_push_particles: subclass "1vb_double"

struct psc_push_particles_ops psc_push_particles_1vb_double_ops = {
  .name                  = "1vb_double",
  .push_mprts_yz         = psc_push_particles_push_mprts_yz,
  .particles_type        = PARTICLE_TYPE,
  .fields_type           = FIELDS_TYPE,
};

