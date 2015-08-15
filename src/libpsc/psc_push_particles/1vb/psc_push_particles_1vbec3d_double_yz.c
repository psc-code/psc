
#include "psc_push_particles_private.h"

#include "psc_particles_as_double.h"
#include "psc_fields_as_c.h"

#include "../inc_defs.h"

#define DIM DIM_YZ
#define CALC_J CALC_J_1VB_VAR1
#define INTERPOLATE_1ST INTERPOLATE_1ST_EC

#define NOT_STATIC

#define psc_push_particles_push_a_yz psc_push_particles_1vbec3d_double_push_a_yz

#include "1vb_yz.c"

