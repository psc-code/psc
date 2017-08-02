
#include "psc_push_particles_private.h"

#include "psc_particles_as_single_by_block.h"
#include "psc_fields_as_single.h"

#include "../inc_defs.h"

#define DIM DIM_XYZ
#define CALC_J CALC_J_1VB_SPLIT
#define INTERPOLATE_1ST INTERPOLATE_1ST_EC
#define PUSHER_BY_BLOCK

#define psc_push_particles_push_mprts_xyz psc_push_particles_1vbec_single_by_block_push_mprts_xyz

#include "1vb.c"

